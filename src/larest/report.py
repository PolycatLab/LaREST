"""Run summary and collected results for the LaREST pipeline.

After each monomer is processed, LaREST writes:

* ``<output>/larest_summary.txt`` — a human-readable overview of the run: the
  key settings, the monomer sanity check, the initiator used, which
  calculations completed, the reaction thermodynamics and the energy matrices.
* ``<output>/results/<monomer>/`` — the files most users need, collected in one
  place: the lowest free energy structure of the monomer, initiator and each
  polymer (``.xyz``), the summary CSVs and energy matrices, and a
  machine-readable ``results.json``.
"""

from __future__ import annotations

import json
import logging
import math
import shutil
from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from typing import TYPE_CHECKING, Any

from larest.constants import CENSO_SECTIONS, ENERGY_MATRIX_PARAMS, MIN_SUFFIX
from larest.crest import get_best_structure
from larest.energy_matrix import format_energy_matrix, low_r2_warning
from larest.output import create_dir

if TYPE_CHECKING:
    from pathlib import Path

    import pandas as pd

    from larest.data import MolResults, MonomerInfo

logger = logging.getLogger(__name__)

SUMMARY_FILE = "larest_summary.txt"
RESULTS_DIR = "results"

_INDENT = "  "
# Sections in order of preference when choosing the "final" results to report.
_FINAL_SECTIONS: list[str] = [
    f"censo_refinement{MIN_SUFFIX}",
    f"censo_optimization{MIN_SUFFIX}",
    f"censo_screening{MIN_SUFFIX}",
    "crest",
    "rdkit",
]
_STAGE_SECTIONS: dict[str, str] = {
    "rdkit": "rdkit",
    "crest": "crest",
    "entropy": "censo_corrected",
}


@dataclass
class MonomerRecord:
    """Everything reported about one monomer in the run summary.

    Attributes
    ----------
    smiles : str
        SMILES string of the monomer.
    info : MonomerInfo
        Result of the monomer sanity check.
    initiator : str | None
        SMILES of the initiator used (ROR only).
    initiator_auto : bool
        Whether the initiator was chosen automatically (isodesmic).
    molecules : dict[str, MolResults | None]
        Results for each molecule, keyed by label (``"monomer"``,
        ``"initiator"``, ``"polymer_n<n>"``); ``None`` if the pipeline failed.
    results_table : str
        Formatted reaction thermodynamics table.
    matrices : dict[str, dict[str, pd.DataFrame]]
        Energy matrices keyed by section then quantity.
    warnings : list[str]
        Warnings to highlight in the summary.
    """

    smiles: str
    info: MonomerInfo
    initiator: str | None = None
    initiator_auto: bool = False
    molecules: dict[str, MolResults | None] = field(default_factory=dict)
    results_table: str = ""
    matrices: dict[str, dict[str, pd.DataFrame]] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)


def get_final_section(results: MolResults | None) -> str | None:
    """Return the most refined section with a free energy for a molecule."""
    if results is None:
        return None
    for section in _FINAL_SECTIONS:
        if results.sections.get(section, {}).get("G") is not None:
            return section
    return None


def collect_results(
    output_dir: Path,
    record: MonomerRecord,
    monomer_dir_name: str,
    summary_dir: Path,
) -> Path:
    """Collect the key outputs for one monomer into ``<output>/results/<monomer>/``.

    Parameters
    ----------
    output_dir : Path
        Root output directory for the run.
    record : MonomerRecord
        Summary record of the monomer.
    monomer_dir_name : str
        Name of the monomer's output directory.
    summary_dir : Path
        Directory containing the summary CSVs and energy matrices.

    Returns
    -------
    Path
        The results directory for the monomer.
    """
    results_dir = output_dir / RESULTS_DIR / monomer_dir_name
    create_dir(results_dir)

    for label, results in record.molecules.items():
        if results is None or results.dir_path is None:
            continue
        structure = get_best_structure(results.dir_path)
        if structure is None:
            logger.warning(f"No final structure found for {label} ({results.smiles})")
            continue
        shutil.copy(structure, results_dir / f"{label}.xyz")

    if summary_dir.exists():
        shutil.copytree(summary_dir, results_dir / "summary", dirs_exist_ok=True)

    data = {
        "monomer": asdict(record.info),
        "initiator": record.initiator,
        "initiator_auto": record.initiator_auto,
        "molecules": {
            label: None
            if results is None
            else {"smiles": results.smiles, "sections": results.sections}
            for label, results in record.molecules.items()
        },
        "energy_matrices": {
            section: {
                param: matrix.astype(object).to_dict(orient="index")
                for param, matrix in matrices.items()
            }
            for section, matrices in record.matrices.items()
        },
        "units": {"E/H/G": "J/mol", "S": "J/mol/K", "energy_matrices": "kJ/mol"},
    }
    with open(results_dir / "results.json", "w") as fstream:
        json.dump(_clean_nan(data), fstream, indent=4)
    return results_dir


def _clean_nan(value: object) -> object:
    """Recursively convert NaN to ``None`` and keys to strings for JSON output."""
    if isinstance(value, dict):
        return {str(k): _clean_nan(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_clean_nan(v) for v in value]
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def write_summary(
    output_dir: Path,
    config: dict[str, Any],
    records: list[MonomerRecord],
) -> Path:
    """Write the human-readable run summary to ``<output>/larest_summary.txt``.

    Parameters
    ----------
    output_dir : Path
        Root output directory for the run.
    config : dict[str, Any]
        Full pipeline configuration dict.
    records : list[MonomerRecord]
        Records of the monomers processed so far.

    Returns
    -------
    Path
        Path to the summary file.
    """
    lines: list[str] = [
        "LaREST run summary",
        "==================",
        f"Written:        {datetime.now(tz=UTC).strftime('%Y-%m-%d %H:%M:%S UTC')}",
        f"LaREST version: {_get_version()}",
        f"Output:         {output_dir}",
        "",
        *_format_settings(config),
        "",
        "Units: energies in kJ/mol (per monomer unless stated), entropies in J/mol/K.",
        "Section names: '<stage>' = Boltzmann average over the conformer ensemble,",
        "'<stage>_min' = lowest free energy conformer, 'censo_corrected' = refinement",
        "lowest free energy conformer with the CREST conformational entropy added.",
    ]
    for idx, record in enumerate(records, 1):
        lines.extend(["", "", f"Monomer [{idx}/{len(records)}]: {record.smiles}"])
        lines.append("-" * len(lines[-1]))
        lines.extend(_format_record(record, config))

    summary_file = output_dir / SUMMARY_FILE
    summary_file.write_text("\n".join(lines) + "\n")
    logger.debug(f"Wrote run summary to {summary_file}")
    return summary_file


def _get_version() -> str:
    try:
        return version("larest")
    except PackageNotFoundError:
        return "unknown"


def _format_settings(config: dict[str, Any]) -> list[str]:
    reaction = config["reaction"]
    steps = config["steps"]
    thermo = config.get("thermo", {})
    lines = [
        "Settings",
        "--------",
        f"Reaction type:  {reaction['type']}",
        f"Lengths:        {reaction['lengths']}",
    ]
    if reaction["type"] == "ROR":
        lines.append(f"Initiator:      {reaction.get('initiator') or '-'}")
    lines.extend(
        [
            f"Temperature:    {thermo.get('temperature', 298.15)} K",
            f"Solvent:        {_format_solvent(config)}",
            "Stages:         "
            + ", ".join(
                f"{stage} {'on' if enabled else 'off'}"
                for stage, enabled in steps.items()
            ),
        ],
    )
    if steps.get("rdkit"):
        rdkit = config.get("rdkit", {})
        lines.append(
            f"RDKit:          {rdkit.get('n_conformers')} conformers, "
            f"{rdkit.get('mmff')}",
        )
    if steps.get("rdkit") or steps.get("xtb"):
        xtb = config.get("xtb", {})
        level = "GFN-FF" if xtb.get("gfnff") else f"GFN{xtb.get('gfn', 2)}-xTB"
        lines.append(f"xTB:            {level}, ohess {xtb.get('ohess', '-')}")
    if steps.get("crest_confgen"):
        lines.append(f"CREST confgen:  {_format_crest(config, 'confgen')}")
    if steps.get("censo"):
        censo = config.get("censo", {})
        general = censo.get("general", {})
        lines.append(
            f"CENSO:          mRRHO {censo.get('screening', {}).get('gfnv', 'gfn2')}, "
            f"rotor cutoff sthr {general.get('sthr', 'default')}",
        )
        lines.extend(
            f"  {part:<14}{_format_censo_part(censo.get(part, {}))}"
            for part in ("prescreening", "screening", "optimization", "refinement")
            if censo.get("cli", {}).get(part)
        )
    if steps.get("crest_entropy"):
        lines.append(f"CREST entropy:  {_format_crest(config, 'entropy')}")
    return lines


def _format_solvent(config: dict[str, Any]) -> str:
    general = config.get("censo", {}).get("general", {})
    if general.get("gas_phase") or general.get("gas-phase"):
        censo = "gas phase"
    else:
        censo = general.get("solvent", "-")
    xtb = config.get("xtb", {}).get("alpb") or "gas phase"
    return f"{config.get('thermo', {}).get('solvent', '-')} (xTB/CREST: {xtb}, CENSO: {censo})"


def _format_crest(config: dict[str, Any], mode: str) -> str:
    crest = config.get("crest", {}).get(mode, {})
    level = "GFN-FF" if crest.get("gfnff") else "GFN2-xTB"
    thresholds = ", ".join(
        f"{key} {crest.get(key, 'default')}" for key in ("ewin", "rthr", "ethr")
    )
    return f"{level}, {thresholds}"


def _format_censo_part(part: dict[str, Any]) -> str:
    func = part.get("func", "default")
    basis = part.get("basis")
    return f"{func}/{basis}" if basis else str(func)


def _format_record(record: MonomerRecord, config: dict[str, Any]) -> list[str]:
    info = record.info
    lines = [
        "Monomer check:",
        f"{_INDENT}Rings:              {info.n_rings} (sizes: "
        f"{', '.join(map(str, info.ring_sizes)) or '-'})",
        f"{_INDENT}Functional groups:  "
        + (
            ", ".join(f"{fg} x{count}" for fg, count in info.functional_groups.items())
            or "none"
        ),
    ]
    if info.valid:
        lines.append(
            f"{_INDENT}Opened:             {info.group} bond "
            f"(ring size {info.ring_size})",
        )
    else:
        lines.append(f"{_INDENT}Status:             SKIPPED - {info.message}")
        return lines
    if info.message:
        lines.append(f"{_INDENT}Note:               {info.message}")
    lines.append(
        f"{_INDENT}Unassigned stereo:  {', '.join(info.unassigned_stereo) or 'none'}",
    )
    if config["reaction"]["type"] == "ROR":
        auto = " (automatic, isodesmic)" if record.initiator_auto else ""
        lines.append(f"Initiator:          {record.initiator}{auto}")

    lines.append("Calculations:")
    for label, results in record.molecules.items():
        lines.append(f"{_INDENT}{label:<18}{_format_status(results, config)}")

    if record.results_table:
        lines.extend(["", "Reaction thermodynamics:", record.results_table])

    if record.matrices:
        lines.extend(
            [
                "",
                "Energy matrices: column m=0 is the reaction energy per monomer, the",
                "diagonal (m=n) the total reaction energy, and 0<m<n the incremental",
                "energy per monomer inserted into a chain of degree m. Matrices for all",
                "sections are in summary/energy_matrix/.",
            ],
        )
    for section, matrices in _report_matrices(record).items():
        for param in ENERGY_MATRIX_PARAMS:
            if param not in matrices or matrices[param].drop(columns="R2").isna().all(
                axis=None,
            ):
                continue
            lines.extend(
                [
                    "",
                    f"Energy matrix {param} ({section}, kJ/mol; rows n, columns m):",
                    format_energy_matrix(matrices[param], indent=_INDENT),
                ],
            )
            if low_r2_warning(matrices[param]):
                lines.append(
                    f"{_INDENT}WARNING: R² of the 1/n extrapolation is below 0.8",
                )

    if record.warnings:
        lines.extend(["", "Warnings:", *(f"{_INDENT}- {w}" for w in record.warnings)])
    return lines


def _report_matrices(record: MonomerRecord) -> dict[str, dict[str, pd.DataFrame]]:
    """Select the energy matrices shown in the summary.

    These are the matrices for the most refined stage completed (lowest free
    energy conformer and, for CENSO, the Boltzmann average) and the
    entropy-corrected results.  All matrices are written to CSV regardless.
    """
    final = get_final_section(record.molecules.get("monomer"))
    if final is None:
        return {}
    sections = [final]
    if final.endswith(MIN_SUFFIX):
        sections.append(final.removesuffix(MIN_SUFFIX))
    sections.append("censo_corrected")
    return {
        section: record.matrices[section]
        for section in sections
        if section in record.matrices
    }


def _format_status(results: MolResults | None, config: dict[str, Any]) -> str:
    if results is None:
        return "FAILED (see larest.log)"
    stages = {
        "rdkit": "rdkit",
        "crest": "crest_confgen",
        "censo": "censo",
        "entropy": "crest_entropy",
    }
    status = []
    for label, step in stages.items():
        if not config["steps"].get(step):
            continue
        if label == "entropy" and not config["steps"].get("censo"):
            done = "crest_entropy" in results.sections
        elif label == "censo":
            done = any(
                results.sections.get(f"{section}{MIN_SUFFIX}", {}).get("G") is not None
                for section in CENSO_SECTIONS
            )
        else:
            done = results.sections.get(_STAGE_SECTIONS[label], {}).get("G") is not None
        status.append(f"{label} {'done' if done else 'missing'}")
    smiles = f"  {results.smiles}"
    return ", ".join(status) + smiles
