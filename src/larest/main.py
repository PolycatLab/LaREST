"""Top-level pipeline orchestration for LaREST.

This module contains the main functions that drive a full LaREST run:

* :func:`run_pipeline` — executes the four stages (RDKit, CREST confgen,
  CENSO, CREST entropy) for a single molecule and writes ``results.json``.
* :func:`compile_results` — combines monomer, polymer, and (for ROR)
  initiator results into per-section ``delta_E / delta_H / delta_S / delta_G``
  CSV summaries inside ``<output>/Monomer/<dir>/summary/``.
* :func:`write_energy_matrices` — computes the energy matrices (see
  :mod:`larest.energy_matrix`) and writes them to ``summary/energy_matrix/``.
* :func:`main` — checks each configured monomer, chooses its initiator, runs
  the pipeline for every molecule, compiles the results and writes the run
  summary (see :mod:`larest.report`).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from larest.censo import check_censo_config, run_censo
from larest.checkpoint import (
    PipelineStage,
    combine_entropy,
    restore_results,
    write_stage_settings,
)
from larest.chem import (
    AUTO_INITIATOR,
    analyse_monomer,
    build_polymer,
    get_isodesmic_initiator,
)
from larest.constants import ENERGY_MATRIX_PARAMS, THERMODYNAMIC_PARAMS
from larest.crest import run_crest_confgen, run_crest_entropy
from larest.data import Initiator, MolResults, Monomer, MonomerInfo, Polymer
from larest.energy_matrix import (
    compute_energy_matrices,
    extrapolate_to_infinite_length,
    low_r2_warning,
)
from larest.output import create_dir, mol_dir_name, remove_dir
from larest.rdkit import run_rdkit
from larest.report import MonomerRecord, collect_results, write_summary
from larest.standard_state import apply_standard_state

logger = logging.getLogger(__name__)

_INDENT = "       "  # 7 spaces — matches LAREST_HEADER box width

_SECTION_LABELS: dict[str, str] = {
    "rdkit": "RDKit",
    "crest": "CREST",
    "censo_prescreening": "CENSO prescreening (Boltzmann average)",
    "censo_screening": "CENSO screening (Boltzmann average)",
    "censo_optimization": "CENSO optimization (Boltzmann average)",
    "censo_refinement": "CENSO refinement (Boltzmann average)",
    "censo_prescreening_min": "CENSO prescreening (lowest conformer)",
    "censo_screening_min": "CENSO screening (lowest conformer)",
    "censo_optimization_min": "CENSO optimization (lowest conformer)",
    "censo_refinement_min": "CENSO refinement (lowest conformer)",
    "censo_corrected": "CENSO corrected (lowest conformer + CREST Sconf)",
}


def get_mol_dir(mol: Monomer | Polymer | Initiator, output_dir: Path) -> Path:
    """Return the output directory of a molecule.

    Parameters
    ----------
    mol : Monomer | Polymer | Initiator
        The molecule.
    output_dir : Path
        Root output directory for the run.

    Returns
    -------
    Path
        ``output_dir/Polymer/<monomer slug>_<length>_<hash>/`` for polymers or
        ``output_dir/<MolType>/<slug>_<hash>/`` otherwise (see
        :func:`~larest.output.mol_dir_name`).
    """
    match mol:
        case Polymer():
            return (
                output_dir
                / "Polymer"
                / mol_dir_name(
                    mol.smiles,
                    suffix=str(mol.length),
                    slug_smiles=mol.monomer_smiles,
                )
            )
        case _:
            return output_dir / type(mol).__name__ / mol_dir_name(mol.smiles)


def run_pipeline(
    mol: Monomer | Polymer | Initiator,
    output_dir: Path,
    config: dict[str, Any],
) -> MolResults:
    """Run the full four-stage pipeline for a single molecule.

    Restores any existing checkpoint results, then conditionally executes
    each enabled stage (controlled by the ``[steps]`` config section) before
    writing a final ``results.json`` file.

    Parameters
    ----------
    mol : Monomer | Polymer | Initiator
        The molecule to process.  Output is placed in the directory given by
        :func:`get_mol_dir`.
    output_dir : Path
        Root output directory for the run.
    config : dict[str, Any]
        Full pipeline configuration dict loaded from ``config.toml``.

    Returns
    -------
    MolResults
        Dataclass containing the molecule's SMILES and a mapping of pipeline
        section names to thermodynamic parameter dicts.
    """
    dir_path = get_mol_dir(mol, output_dir)

    create_dir(output_dir / "temp")
    temperature: float = (
        config.get("censo", {}).get("general", {}).get("temperature", 298.15)
    )
    steps = config["steps"]
    results, stage = restore_results(
        dir_path,
        temperature=temperature,
        steps=steps,
        config=config,
        smiles=mol.smiles,
    )

    if steps["rdkit"] and stage <= PipelineStage.RDKIT:
        results |= run_rdkit(mol.smiles, dir_path, config)
        write_stage_settings(dir_path, PipelineStage.RDKIT, config, mol.smiles)
    if steps["crest_confgen"] and stage <= PipelineStage.CREST_CONFGEN:
        results |= run_crest_confgen(dir_path, config)
        write_stage_settings(dir_path, PipelineStage.CREST_CONFGEN, config, mol.smiles)
    if steps["censo"] and stage <= PipelineStage.CENSO:
        results |= run_censo(dir_path, output_dir, config)
        write_stage_settings(dir_path, PipelineStage.CENSO, config, mol.smiles)
    if steps["crest_entropy"] and stage <= PipelineStage.CREST_ENTROPY:
        crest_entropy_results = run_crest_entropy(dir_path, config)
        write_stage_settings(dir_path, PipelineStage.CREST_ENTROPY, config, mol.smiles)
        results |= combine_entropy(
            results,
            crest_entropy_results,
            temperature,
            censo_enabled=steps["censo"],
        )

    # stage results are stored as computed (1 atm thermostatistics); the
    # standard-state conversion is applied to the reported results only
    results = apply_standard_state(results, config)

    results_file = dir_path / "results.json"
    logger.debug(f"Writing final results to {results_file}")
    with open(results_file, "w") as fstream:
        json.dump(results, fstream, sort_keys=True, indent=4, allow_nan=True)

    remove_dir(output_dir / "temp")
    return MolResults(smiles=mol.smiles, sections=results, dir_path=dir_path)


def format_results_table(
    monomer_smiles: str,
    summary_dfs: dict[str, pd.DataFrame],
) -> str:
    """Format pipeline delta results as a human-readable table string.

    Parameters
    ----------
    monomer_smiles : str
        SMILES string of the monomer (shown in the header).
    summary_dfs : dict[str, pd.DataFrame]
        Mapping of section name to its summary DataFrame, as returned by
        :func:`compile_results`.  Each DataFrame must contain columns
        ``polymer_length``, ``delta_H``, ``delta_S``, and ``delta_G``, and
        optionally ``delta_E``.

    Returns
    -------
    str
        Multi-section formatted table string ready to be printed.
    """
    border = "=" * 59
    lines = [
        "",
        f"{_INDENT}{border}",
        f"{_INDENT}Results: {monomer_smiles}",
        f"{_INDENT}{border}",
    ]
    for section, df in summary_dfs.items():
        deltas = [
            f"delta_{param}"
            for param in THERMODYNAMIC_PARAMS
            if f"delta_{param}" in df.columns
        ]
        display = df[["polymer_length", *deltas]].copy()
        if display[deltas].isna().all(axis=None):
            continue

        for col in ("delta_E", "delta_H", "delta_G"):
            if col in display.columns:
                display[col] = display[col] / 1000
        display["polymer_length"] = display["polymer_length"].apply(
            lambda x: "∞" if np.isinf(x) else str(int(x)),
        )
        display = display.rename(
            columns={
                "polymer_length": "n",
                "delta_E": "ΔE (kJ/mol)",
                "delta_H": "ΔH (kJ/mol)",
                "delta_S": "ΔS (J/mol/K)",
                "delta_G": "ΔG (kJ/mol)",
            },
        )

        label = _SECTION_LABELS.get(section, section)
        lines.append(f"\n{_INDENT}{label}")

        table_lines = display.to_string(
            index=False,
            float_format="{:.4f}".format,
        ).split("\n")
        for i, line in enumerate(table_lines):
            if line.lstrip().startswith("∞"):
                table_lines.insert(i, "-" * len(table_lines[0]))
                break
        lines.append("\n".join(f"{_INDENT}  {line}" for line in table_lines))

    return "\n".join(lines)


def compile_results(
    monomer_smiles: str,
    monomer_results: MolResults,
    polymer_results: list[tuple[int, MolResults]],
    initiator_results: MolResults | None,
    output_dir: Path,
    reaction_type: str,
) -> dict[str, pd.DataFrame]:
    """Compute reaction thermodynamics and write per-section CSV summaries.

    For each pipeline section present in *monomer_results*, constructs a
    DataFrame containing the raw monomer, initiator, and polymer thermodynamic
    values at every polymer length, then appends columns for
    ``delta_E``, ``delta_H``, ``delta_S``, and ``delta_G`` computed as::

        delta_param = (polymer_param - n * monomer_param - initiator_param) / n

    where *n* is the polymer length.  A final row extrapolates each delta to
    infinite chain length by linear regression against 1/n.  One CSV file per
    section is written to ``<monomer dir>/summary/``.

    Parameters
    ----------
    monomer_smiles : str
        SMILES string of the monomer (used to locate the summary directory if
        *monomer_results* has no ``dir_path``).
    monomer_results : MolResults
        Pipeline results for the monomer molecule.
    polymer_results : list[tuple[int, MolResults]]
        Pairs of (polymer_length, MolResults) for each requested chain length.
    initiator_results : MolResults | None
        Pipeline results for the initiator.  Required when *reaction_type* is
        ``"ROR"``; pass ``None`` for ``"RER"`` reactions.
    output_dir : Path
        Root output directory for the run.
    reaction_type : str
        Either ``"ROR"`` or ``"RER"``.  Determines whether initiator
        thermodynamics are included in the delta calculation.

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping of section name to its summary DataFrame (columns:
        ``polymer_length``, ``monomer_*``, ``initiator_*``, ``polymer_*``,
        ``delta_E``, ``delta_H``, ``delta_S``, ``delta_G``).

    Raises
    ------
    ValueError
        If *reaction_type* is ``"ROR"`` but *initiator_results* is ``None``.
    """
    if reaction_type == "ROR" and initiator_results is None:
        raise ValueError(
            f"Initiator results are required for ROR reactions but are missing (monomer: {monomer_smiles})",
        )

    summary_dir = _get_summary_dir(monomer_smiles, monomer_results, output_dir)
    create_dir(summary_dir)

    summary_headings: list[str] = ["polymer_length"]
    for mol_type in ["monomer", "initiator", "polymer"]:
        summary_headings.extend(
            [f"{mol_type}_{param}" for param in THERMODYNAMIC_PARAMS],
        )

    summary_dfs: dict[str, pd.DataFrame] = {}
    initiator_sections = (
        initiator_results.sections
        if reaction_type == "ROR" and initiator_results is not None
        else None
    )

    for section in _thermodynamic_sections(monomer_results):
        summary: dict[str, list[float | None]] = {
            heading: [] for heading in summary_headings
        }

        for length, poly_results in polymer_results:
            summary["polymer_length"].append(length)
            for param in THERMODYNAMIC_PARAMS:
                summary[f"polymer_{param}"].append(
                    poly_results.sections.get(section, {}).get(param),
                )
                summary[f"monomer_{param}"].append(
                    monomer_results.sections[section].get(param),
                )
                summary[f"initiator_{param}"].append(
                    initiator_sections.get(section, {}).get(param)
                    if initiator_sections is not None
                    else 0,
                )

        summary_df = pd.DataFrame(
            summary,
            index=None,
            dtype=np.float64,
        ).sort_values("polymer_length", ascending=True)

        for param in THERMODYNAMIC_PARAMS:
            summary_df[f"delta_{param}"] = (
                summary_df[f"polymer_{param}"]
                - (summary_df["polymer_length"] * summary_df[f"monomer_{param}"])
                - summary_df[f"initiator_{param}"]
            ) / summary_df["polymer_length"]

        inf_row: dict[str, float] = {"polymer_length": np.inf}
        for mol_type in ["monomer", "initiator", "polymer"]:
            for param in THERMODYNAMIC_PARAMS:
                inf_row[f"{mol_type}_{param}"] = np.nan
        for param in THERMODYNAMIC_PARAMS:
            col = f"delta_{param}"
            inf_row[col], _ = extrapolate_to_infinite_length(
                summary_df["polymer_length"].tolist(),
                summary_df[col].tolist(),
            )
        summary_df = pd.concat(
            [summary_df, pd.DataFrame([inf_row])],
            ignore_index=True,
        )

        summary_df.to_csv(summary_dir / f"{section}.csv", header=True, index=False)
        summary_dfs[section] = summary_df

    return summary_dfs


def _get_summary_dir(
    monomer_smiles: str,
    monomer_results: MolResults,
    output_dir: Path,
) -> Path:
    monomer_dir = monomer_results.dir_path or get_mol_dir(
        Monomer(smiles=monomer_smiles),
        output_dir,
    )
    return monomer_dir / "summary"


def _thermodynamic_sections(mol_results: MolResults) -> list[str]:
    """Return the sections of *mol_results* that hold H/S/G-type parameters."""
    return [
        section
        for section, params in mol_results.sections.items()
        if set(params) & set(THERMODYNAMIC_PARAMS)
    ]


def write_energy_matrices(
    monomer_smiles: str,
    monomer_results: MolResults,
    polymer_results: list[tuple[int, MolResults]],
    initiator_results: MolResults | None,
    output_dir: Path,
    reaction_type: str,
) -> dict[str, dict[str, pd.DataFrame]]:
    """Compute energy matrices for every section and write them as CSV files.

    One CSV per section and quantity (E, H, G) is written to
    ``<monomer dir>/summary/energy_matrix/<section>_<param>.csv`` (kJ/mol).
    See :mod:`larest.energy_matrix` for the definitions of the descriptors.

    Parameters
    ----------
    monomer_smiles : str
        SMILES string of the monomer.
    monomer_results : MolResults
        Pipeline results for the monomer molecule.
    polymer_results : list[tuple[int, MolResults]]
        Pairs of (polymer_length, MolResults) for each requested chain length.
    initiator_results : MolResults | None
        Pipeline results for the initiator (ROR only).
    output_dir : Path
        Root output directory for the run.
    reaction_type : str
        Either ``"ROR"`` or ``"RER"``.

    Returns
    -------
    dict[str, dict[str, pd.DataFrame]]
        Energy matrices keyed by section then quantity.
    """
    matrix_dir = (
        _get_summary_dir(monomer_smiles, monomer_results, output_dir) / "energy_matrix"
    )
    create_dir(matrix_dir)

    all_matrices: dict[str, dict[str, pd.DataFrame]] = {}
    for section in _thermodynamic_sections(monomer_results):
        matrices = compute_energy_matrices(
            monomer=monomer_results.sections[section],
            polymers={
                length: results.sections.get(section, {})
                for length, results in polymer_results
            },
            initiator=initiator_results.sections.get(section, {})
            if initiator_results is not None
            else None,
            reaction_type=reaction_type,
            params=ENERGY_MATRIX_PARAMS,
        )
        for param, matrix in matrices.items():
            matrix.to_csv(matrix_dir / f"{section}_{param}.csv")
            if low_r2_warning(matrix):
                logger.warning(
                    f"R² of the 1/n extrapolation of {param} ({section}) for "
                    f"{monomer_smiles} is below 0.8",
                )
        all_matrices[section] = matrices
    return all_matrices


def _choose_initiator(
    monomer_smiles: str,
    group: str | None,
    config: dict[str, Any],
) -> tuple[str, bool]:
    """Return the initiator SMILES for a monomer and whether it was chosen automatically."""
    initiator: str = config["reaction"]["initiator"]
    if initiator.strip().lower() == AUTO_INITIATOR:
        return get_isodesmic_initiator(monomer_smiles, group=group), True
    return initiator, False


def main(output_dir: Path, config: dict[str, Any]) -> None:
    """Run the complete LaREST pipeline for all configured monomers and lengths.

    Iterates over every monomer in ``config["reaction"]["monomers"]``.  Each
    monomer is first checked (:func:`~larest.chem.analyse_monomer`) and skipped
    if it cannot be ring-opened.  For ROR, the initiator is either the
    configured SMILES or, with ``initiator = "auto"``, the isodesmic initiator
    for that monomer.  The pipeline is then run for the monomer, initiator, and
    every requested polymer length, and the results are compiled into summary
    CSVs and energy matrices.  After each monomer, the run summary
    (``larest_summary.txt``) and the collected results (``results/``) are
    updated.  Failures for individual molecules or polymer lengths are logged
    and skipped so that the remaining molecules can still be processed.

    Parameters
    ----------
    output_dir : Path
        Root output directory; sub-directories are created automatically.
    config : dict[str, Any]
        Full pipeline configuration dict loaded from ``config.toml``.
    """
    reaction_type = config["reaction"]["type"]
    result_tables: list[str] = []
    records: list[MonomerRecord] = []

    monomers: list[dict[str, Any]] = [
        {"smiles": m, "group": None} if isinstance(m, str) else m
        for m in config["reaction"]["monomers"]
    ]
    n_monomers = len(monomers)
    lengths = config["reaction"]["lengths"]

    print(f"{_INDENT}Reaction type:   {reaction_type}")  # noqa: T201
    print(f"{_INDENT}Monomers:        {n_monomers}")  # noqa: T201
    if reaction_type == "ROR":
        print(f"{_INDENT}Initiator:       {config['reaction']['initiator']}")  # noqa: T201
    print(f"{_INDENT}Polymer lengths: {lengths}")  # noqa: T201
    print(f"{_INDENT}Output:          {output_dir}")  # noqa: T201

    if config["steps"].get("censo"):
        check_censo_config(config)

    if reaction_type == "ROR" and not config["reaction"].get("initiator"):
        logger.error("ROR reactions require [reaction].initiator (a SMILES or 'auto')")
        return

    logger.info("Running pipeline for monomers")
    for monomer_idx, monomer_entry in enumerate(monomers, 1):
        monomer_smiles: str = monomer_entry["smiles"]
        group: str | None = monomer_entry.get("group")
        print(f"\n{_INDENT}Monomer [{monomer_idx}/{n_monomers}]: {monomer_smiles}")  # noqa: T201

        try:
            info = analyse_monomer(monomer_smiles, group=group)
        except ValueError:
            logger.exception(f"Invalid monomer SMILES {monomer_smiles}")
            info = MonomerInfo(smiles=monomer_smiles, message="invalid SMILES")
        record = MonomerRecord(smiles=monomer_smiles, info=info)
        records.append(record)
        if not info.valid:
            print(f"{_INDENT}  Skipped: {info.message}")  # noqa: T201
            write_summary(output_dir, config, records)
            continue
        if info.unassigned_stereo:
            record.warnings.append(
                "unassigned stereochemistry "
                f"({', '.join(info.unassigned_stereo)}); conformers may mix "
                "stereoisomers",
            )

        monomer = Monomer(smiles=monomer_smiles)
        print(f"{_INDENT}  Monomer...")  # noqa: T201
        try:
            monomer_results = run_pipeline(monomer, output_dir, config)
        except Exception:
            logger.exception(f"Failed to run pipeline for monomer {monomer_smiles}")
            record.molecules["monomer"] = None
            write_summary(output_dir, config, records)
            continue
        record.molecules["monomer"] = monomer_results

        initiator_smiles = None
        initiator_results = None
        if reaction_type == "ROR":
            try:
                initiator_smiles, record.initiator_auto = _choose_initiator(
                    monomer_smiles,
                    info.group,
                    config,
                )
            except Exception:
                logger.exception(f"Failed to choose initiator for {monomer_smiles}")
                write_summary(output_dir, config, records)
                continue
            record.initiator = initiator_smiles
            initiator = Initiator(smiles=initiator_smiles)
            print(f"{_INDENT}  Initiator ({initiator_smiles})...")  # noqa: T201
            try:
                initiator_results = run_pipeline(initiator, output_dir, config)
            except Exception:
                logger.exception(
                    f"Failed to run pipeline for initiator {initiator_smiles}",
                )
                record.molecules["initiator"] = None
                write_summary(output_dir, config, records)
                continue
            record.molecules["initiator"] = initiator_results

        polymers = []
        for length in lengths:
            try:
                polymer_smiles = build_polymer(
                    monomer_smiles=monomer_smiles,
                    polymer_length=length,
                    reaction_type=reaction_type,
                    initiator_smiles=initiator_smiles,
                    group=info.group,
                )
            except Exception:
                logger.exception(
                    f"Failed to build polymer for monomer {monomer_smiles} (length: {length})",
                )
                record.warnings.append(f"failed to build polymer of length {length}")
                continue
            polymers.append(
                Polymer(
                    smiles=polymer_smiles,
                    monomer_smiles=monomer_smiles,
                    length=length,
                    initiator_smiles=initiator_smiles,
                ),
            )

        polymer_results: list[tuple[int, MolResults]] = []
        n_polymers = len(polymers)
        for polymer_idx, polymer in enumerate(polymers, 1):
            print(  # noqa: T201
                f"{_INDENT}  Polymer n={polymer.length} [{polymer_idx}/{n_polymers}]...",
            )
            label = f"polymer_n{polymer.length}"
            try:
                result = run_pipeline(polymer, output_dir, config)
                polymer_results.append((polymer.length, result))
                record.molecules[label] = result
            except Exception:
                logger.exception(
                    f"Failed to run pipeline for polymer {monomer_smiles} (length: {polymer.length})",
                )
                record.molecules[label] = None
                record.warnings.append(
                    f"pipeline failed for polymer of length {polymer.length}; "
                    "it is excluded from the results",
                )
                continue

        print(f"{_INDENT}  Compiling results...")  # noqa: T201
        summary_dfs = compile_results(
            monomer_smiles=monomer_smiles,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=initiator_results,
            output_dir=output_dir,
            reaction_type=reaction_type,
        )
        table = format_results_table(monomer_smiles, summary_dfs)
        result_tables.append(table)
        record.results_table = table
        record.matrices = write_energy_matrices(
            monomer_smiles=monomer_smiles,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=initiator_results,
            output_dir=output_dir,
            reaction_type=reaction_type,
        )
        try:
            collect_results(
                output_dir=output_dir,
                record=record,
                monomer_dir_name=get_mol_dir(monomer, output_dir).name,
                summary_dir=_get_summary_dir(
                    monomer_smiles,
                    monomer_results,
                    output_dir,
                ),
                use_censo=config["steps"].get("censo", True),
            )
        except Exception:
            logger.exception(f"Failed to collect results for {monomer_smiles}")
        write_summary(output_dir, config, records)

    print()  # noqa: T201
    for table in result_tables:
        print(table)  # noqa: T201
    if records:
        print(f"\n{_INDENT}Summary written to {output_dir / 'larest_summary.txt'}")  # noqa: T201
