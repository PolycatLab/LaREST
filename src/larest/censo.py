"""CENSO DFT refinement stage of the LaREST pipeline.

Wraps the CENSO (Conformer Ensemble Sorting) program, which performs four
sequential sub-stages of increasing DFT accuracy:

0. ``censo_prescreening`` — fast pre-screening of the CREST ensemble
1. ``censo_screening``    — singlepoint DFT screening
2. ``censo_optimization`` — DFT geometry optimisation
3. ``censo_refinement``   — high-accuracy single-point refinement

Thermodynamic parameters (E, H, G) are extracted per conformer from the
``<n>_<PART>.json`` files written by CENSO at each sub-stage, with the enthalpy
correction taken from CENSO's xTB mRRHO calculations.  Both Boltzmann-averaged
and lowest free energy conformer values are reported, and entropy is derived as
S = (H - G) / T.
"""

from __future__ import annotations

import csv
import importlib
import json
import logging
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

from larest.constants import (
    CENSO_PARTS,
    CENSO_SECTIONS,
    GAS_CONSTANT,
    HARTTREE_TO_JMOL,
    MIN_SUFFIX,
    THERMODYNAMIC_PARAMS,
)
from larest.output import create_dir
from larest.setup import parse_command_args

logger = logging.getLogger(__name__)

_CENSO_PATH_EXECUTABLES: dict[str, str] = {
    "orc": "orca",
    "xtb": "xtb",
    "cosmotherm": "cosmotherm",
    "cosmorssetup": "cosmorssetup",
    "tm": "tm",
}


def _resolve_censo_paths(user_paths: dict[str, str]) -> dict[str, str]:
    fallback = shutil.which("ls") or "/bin/ls"
    return {
        key: user_paths.get(key) or shutil.which(exe) or fallback
        for key, exe in _CENSO_PATH_EXECUTABLES.items()
    }


def _get_censo_config_fields() -> dict[str, set[str]] | None:
    """Return the settings CENSO recognises in each ``.censo2rc`` section.

    Read from CENSO's own configuration models, or ``None`` if they cannot be
    imported (e.g. a different CENSO version).
    """
    try:
        # imported dynamically as CENSO is an optional runtime dependency
        parts_config = importlib.import_module("censo.config.parts_config")
    except Exception:
        logger.debug("Could not import CENSO configuration models")
        return None
    return {
        section: set(field.annotation.model_fields)
        for section, field in parts_config.PartsConfig.model_fields.items()
        if hasattr(field.annotation, "model_fields")
    }


def check_censo_config(config: dict[str, Any]) -> list[str]:
    """Warn about ``[censo.*]`` settings that CENSO will silently ignore.

    CENSO ignores unrecognised keys in ``.censo2rc``, so a misspelt or
    outdated setting (e.g. ``bhess``, which CENSO 3 no longer supports) would
    otherwise have no effect without any warning.

    Parameters
    ----------
    config : dict[str, Any]
        Full pipeline configuration dict.

    Returns
    -------
    list[str]
        The unrecognised settings, as ``"[section] key"`` strings.
    """
    fields = _get_censo_config_fields()
    if fields is None:
        return []
    unknown: list[str] = []
    for header, sub_config in config.get("censo", {}).items():
        if header == "cli" or not isinstance(sub_config, dict):
            continue
        if header not in fields:
            unknown.append(f"[{header}]")
            continue
        unknown.extend(
            f"[{header}] {key}"
            for key in sub_config
            if key.replace("-", "_") not in fields[header]
        )
    if unknown:
        logger.warning(
            "CENSO does not recognise these settings, which will be ignored: %s",
            ", ".join(unknown),
        )
    return unknown


def create_censorc(config: dict[str, Any], temp_dir: Path) -> None:
    """Write the CENSO runtime configuration file (``.censo2rc``) to *temp_dir*.

    Iterates over ``config["censo"]`` and writes each sub-section as an INI
    block, e.g.::

        [general]
        temperature = 298.15
        ...

    Parameters
    ----------
    config : dict[str, Any]
        Full pipeline configuration dict.  The ``[censo]`` section is used.
    temp_dir : Path
        Directory in which the ``.censo2rc`` file is written.
    """
    censorc_file = temp_dir / ".censo2rc"
    censo_config: dict[str, Any] = config["censo"]

    with open(censorc_file, "w") as fstream:
        for header, sub_config in censo_config.items():
            if header in ("cli", "paths"):
                continue
            fstream.write(f"[{header}]\n")
            # CENSO's .censo2rc keys are pydantic field names (valid Python
            # identifiers, so always underscores); unknown keys are silently
            # ignored (pydantic extra="ignore"). Normalise any hyphens to
            # underscores so e.g. `gas-phase` reaches CENSO as `gas_phase`
            # instead of being silently dropped (leaving solvation enabled).
            fstream.writelines(
                f"{key.replace('-', '_')} = {value}\n"
                for key, value in sub_config.items()
            )
            fstream.write("\n")

        paths = _resolve_censo_paths(censo_config.get("paths", {}))
        fstream.write("[paths]\n")
        fstream.writelines(f"{key} = {value}\n" for key, value in paths.items())
        fstream.write("\n")
    logger.debug(f"Created censo config file at {censorc_file}")


def build_rdkit_ensemble_xyz(dir_path: Path, output_xyz_file: Path) -> None:
    """Assemble a multi-conformer XYZ ensemble from the xTB-optimised RDKit conformers.

    Used as the CENSO input when the CREST conformer-generation stage is
    skipped.  Reads the xTB ranking in ``xtb/rdkit/results.csv``, then
    concatenates each conformer's xTB-optimised geometry
    (``xtb/rdkit/conformer_<id>/conformer_<id>.xtbopt.xyz``) into a single
    multi-conformer XYZ file, ordered by ascending free energy ``G`` so the
    lowest-energy conformer appears first.

    Parameters
    ----------
    dir_path : Path
        Root molecule output directory (e.g. ``output/Monomer/<slug>``).  Must
        contain a completed RDKit stage.
    output_xyz_file : Path
        Destination path for the assembled multi-conformer XYZ file.

    Raises
    ------
    FileNotFoundError
        If the RDKit xTB results CSV is missing.
    RuntimeError
        If no xTB-optimised conformer geometries could be found.
    """
    results_file = dir_path / "xtb" / "rdkit" / "results.csv"
    if not results_file.exists():
        raise FileNotFoundError(
            f"Cannot build CENSO input: RDKit results not found at {results_file}. "
            "Enable [steps].rdkit or [steps].crest_confgen.",
        )

    with open(results_file, newline="") as fstream:
        rows = [row for row in csv.DictReader(fstream) if row.get("conformer_id")]

    def _sort_key(row: dict[str, str]) -> tuple[bool, float]:
        try:
            g = float(row["G"])
        except (TypeError, ValueError):
            g = math.nan
        # Conformers with an undefined G are placed last.
        return (math.isnan(g), g)

    rows.sort(key=_sort_key)

    n_written = 0
    with open(output_xyz_file, "w") as fout:
        for row in rows:
            conformer_id = int(float(row["conformer_id"]))
            conformer_xyz_file = (
                dir_path
                / "xtb"
                / "rdkit"
                / f"conformer_{conformer_id}"
                / f"conformer_{conformer_id}.xtbopt.xyz"
            )
            if not conformer_xyz_file.exists():
                logger.warning(
                    f"Skipping conformer {conformer_id}: missing {conformer_xyz_file}",
                )
                continue
            block = conformer_xyz_file.read_text()
            fout.write(block if block.endswith("\n") else block + "\n")
            n_written += 1

    if n_written == 0:
        raise RuntimeError(
            f"Failed to build CENSO input: no xTB-optimised RDKit conformers found "
            f"under {dir_path / 'xtb' / 'rdkit'}",
        )

    logger.debug(
        f"Wrote {n_written} RDKit conformers to CENSO input file {output_xyz_file}",
    )


def run_censo(
    dir_path: Path,
    output_dir: Path,
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    """Run CENSO on the CREST conformer ensemble and return thermodynamic results.

    Creates the ``.censo2rc`` config file, invokes the ``censo`` binary with
    the CREST conformers XYZ as input, parses the output for all four
    sub-stages, writes ``<dir_path>/censo/results.json``, and returns the
    parsed results.

    Parameters
    ----------
    dir_path : Path
        Root molecule output directory.  Uses
        ``crest_confgen/crest_conformers.xyz`` from a completed CREST stage when
        present; if the CREST stage was skipped (``[steps].crest_confgen =
        false``), falls back to building the conformer ensemble from the
        xTB-optimised RDKit conformers instead.
    output_dir : Path
        Root run output directory; the ``.censo2rc`` file is written to
        ``<output_dir>/temp/``.
    config : dict[str, Any]
        Full pipeline configuration dict.  Uses ``[censo]`` and
        ``[censo][cli]`` sub-sections.

    Returns
    -------
    dict[str, dict[str, float | None]]
        Mapping of CENSO sub-stage name to thermodynamic parameter dict, e.g.
        ``{"censo_prescreening": {"H": ..., "S": ..., "G": ...}, ...}``.

    Raises
    ------
    subprocess.CalledProcessError
        If the CENSO process exits with a non-zero return code.
    RuntimeError
        If CENSO wrote no per-stage results files.
    """
    temp_dir = output_dir / "temp"
    create_censorc(config=config, temp_dir=temp_dir)

    censo_dir = dir_path / "censo"
    create_dir(censo_dir)

    censo_output_file = censo_dir / "censo.txt"
    censo_config_file = temp_dir / ".censo2rc"

    # CENSO sets its working directory to the parent of --input, so place the
    # conformers file into censo_dir to ensure all CENSO output lands there.
    crest_conformers_src = dir_path / "crest_confgen" / "crest_conformers.xyz"
    censo_input_file = censo_dir / "crest_conformers.xyz"
    if crest_conformers_src.exists():
        shutil.copy(crest_conformers_src, censo_input_file)
    else:
        # CREST stage was skipped; fall back to the RDKit conformer ensemble.
        logger.info(
            f"No CREST ensemble at {crest_conformers_src}; "
            "building CENSO input from RDKit conformers",
        )
        build_rdkit_ensemble_xyz(dir_path=dir_path, output_xyz_file=censo_input_file)

    censo_args: list[str] = [
        "censo",
        "--input",
        str(censo_input_file.absolute()),
        "--inprc",
        str(censo_config_file.absolute()),
        *parse_command_args(sub_config=["censo", "cli"], config=config),
    ]

    with open(censo_output_file, "w") as fstream:
        subprocess.run(
            censo_args,
            stdout=fstream,
            stderr=subprocess.STDOUT,
            cwd=censo_dir,
            check=True,
        )

    if get_last_censo_section(censo_dir) is None:
        raise RuntimeError(
            f"CENSO finished but wrote no <n>_<PART>.json results in {censo_dir}; "
            "check the CENSO version (LaREST expects CENSO 3) and censo.txt",
        )
    censo_results: dict[str, dict[str, float | None]] = parse_censo_results(
        censo_dir=censo_dir,
        temperature=config["censo"]["general"].get("temperature", 298.15),
    )

    censo_results_file = censo_dir / "results.json"
    logger.debug(f"Writing results to {censo_results_file}")
    with open(censo_results_file, "w") as fstream:
        json.dump(censo_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return censo_results


def parse_censo_part(
    censo_dir: Path,
    section: str,
) -> list[dict[str, Any]]:
    """Read the per-conformer energies of one CENSO sub-stage.

    CENSO writes ``<n>_<PART>.json`` with, for every conformer, the DFT
    ``energy``, the solvation free energy ``gsolv``, the mRRHO free energy
    correction ``grrho`` and their sum ``gtot`` (all in Hartree).  The matching
    enthalpy correction is not stored by CENSO, so it is taken from the same
    xTB mRRHO calculation (``<n>_<PART>/<CONF>/xtb_rrho/xtb_rrho.out``) as
    ``TOTAL ENTHALPY - TOTAL ENERGY``.

    For every conformer the following are returned (J/mol):

    * ``E = energy + gsolv`` — the electronic energy, including implicit
      solvation.  In the optimisation stage the geometry is optimised in
      solvent, so ``energy`` already includes solvation and ``gsolv`` is zero.
    * ``H = E + H_mRRHO`` — assuming the solvation enthalpy equals ``gsolv``.
    * ``G = E + G_mRRHO`` (CENSO's ``gtot``).

    ``H`` and ``G`` are ``None`` where no mRRHO contribution was evaluated
    (e.g. the prescreening stage), since ``E + gsolv`` is not a free energy.

    Parameters
    ----------
    censo_dir : Path
        CENSO working directory (``<molecule>/censo``).
    section : str
        CENSO section name, a key of :data:`~larest.constants.CENSO_PARTS`.

    Returns
    -------
    list[dict[str, Any]]
        One dict per conformer with keys ``conformer`` (e.g. ``"CONF3"``),
        ``E``, ``H``, ``G`` (J/mol or ``None``) and ``gtot`` (Hartree, used for
        ranking), sorted by ascending ``gtot``.

    Raises
    ------
    FileNotFoundError
        If the sub-stage JSON file does not exist.
    """
    part = CENSO_PARTS[section]
    part_file = censo_dir / f"{part}.json"
    data: dict[str, dict[str, float]] = json.loads(part_file.read_text())["data"]

    conformers: list[dict[str, Any]] = []
    for conformer, values in data.items():
        e = (values["energy"] + values.get("gsolv", 0.0)) * HARTTREE_TO_JMOL
        h_corr = _parse_xtb_rrho_enthalpy(
            censo_dir / part / conformer / "xtb_rrho" / "xtb_rrho.out",
        )
        has_rrho = values.get("grrho", 0.0) != 0.0
        conformers.append(
            {
                "conformer": conformer,
                "E": e,
                "H": e + h_corr * HARTTREE_TO_JMOL
                if (has_rrho and h_corr is not None)
                else None,
                "G": values["gtot"] * HARTTREE_TO_JMOL if has_rrho else None,
                "gtot": values["gtot"],
            },
        )
    conformers.sort(key=lambda conf: conf["gtot"])
    return conformers


def _parse_xtb_rrho_enthalpy(xtb_rrho_output_file: Path) -> float | None:
    """Return the xTB mRRHO enthalpy correction (Hartree) from an xTB output file.

    The correction is ``TOTAL ENTHALPY - TOTAL ENERGY`` (ZPVE plus thermal
    enthalpy), or ``None`` if the file or either value is missing.
    """
    if not xtb_rrho_output_file.exists():
        return None
    energy = enthalpy = None
    with open(xtb_rrho_output_file) as fstream:
        for line in fstream:
            if "TOTAL ENERGY" in line:
                energy = float(line.split()[3])
            elif "TOTAL ENTHALPY" in line:
                enthalpy = float(line.split()[3])
    if energy is None or enthalpy is None:
        logger.warning(f"Failed to extract mRRHO enthalpy from {xtb_rrho_output_file}")
        return None
    return enthalpy - energy


def parse_censo_results(
    censo_dir: Path,
    temperature: float,
) -> dict[str, dict[str, float | None]]:
    """Extract per-sub-stage thermodynamic parameters from a CENSO run.

    For each CENSO sub-stage two sections are returned:

    * ``<section>`` — Boltzmann-averaged ``E``, ``H``, ``G`` over the
      conformer ensemble, weighted by CENSO's ``gtot`` populations at
      *temperature* (as in CENSO's own ``avE``/``avG``).
    * ``<section>_min`` — values of the lowest free energy (``gtot``)
      conformer.

    Entropy is derived as ``S = (H - G) / T``.  See :func:`parse_censo_part`
    for how the energies are defined.

    Parameters
    ----------
    censo_dir : Path
        CENSO working directory (``<molecule>/censo``).
    temperature : float
        Temperature in Kelvin used for the Boltzmann weights and entropy.

    Returns
    -------
    dict[str, dict[str, float | None]]
        Mapping of section name to ``{"E", "H", "S", "G"}`` (J/mol, J/mol/K).
        Values are ``None`` where they could not be computed (sub-stage not
        run, or no mRRHO contribution).
    """
    censo_output: dict[str, dict[str, float | None]] = {}
    for section in CENSO_SECTIONS:
        average: dict[str, float | None] = dict.fromkeys(THERMODYNAMIC_PARAMS, None)
        minimum: dict[str, float | None] = dict.fromkeys(THERMODYNAMIC_PARAMS, None)
        censo_output[section] = average
        censo_output[f"{section}{MIN_SUFFIX}"] = minimum

        try:
            conformers = parse_censo_part(censo_dir, section)
        except FileNotFoundError:
            logger.debug(f"No CENSO results for {section} in {censo_dir}")
            continue
        except Exception:
            logger.exception(f"Failed to parse CENSO results for {section}")
            continue
        if not conformers:
            continue

        g_min = conformers[0]["gtot"]
        weights = [
            math.exp(
                -(conf["gtot"] - g_min)
                * HARTTREE_TO_JMOL
                / (GAS_CONSTANT * temperature),
            )
            for conf in conformers
        ]
        total_weight = sum(weights)
        for param in ("E", "H", "G"):
            values = [conf[param] for conf in conformers]
            minimum[param] = values[0]
            if all(value is not None for value in values):
                average[param] = (
                    sum(w * v for w, v in zip(weights, values, strict=True))
                    / total_weight
                )
        for params in (average, minimum):
            if params["H"] is not None and params["G"] is not None:
                params["S"] = (params["H"] - params["G"]) / temperature
        logger.debug(f"{section}: average {average}, minimum {minimum}")

    return censo_output


def get_best_censo_conformer(censo_dir: Path, section: str) -> str:
    """Return the ID (e.g. ``"CONF5"``) of the lowest free energy conformer.

    Parameters
    ----------
    censo_dir : Path
        CENSO working directory (``<molecule>/censo``).
    section : str
        CENSO section name, a key of :data:`~larest.constants.CENSO_PARTS`.

    Returns
    -------
    str
        Conformer ID with the lowest ``gtot`` in that sub-stage.
    """
    return parse_censo_part(censo_dir, section)[0]["conformer"]


def get_last_censo_section(censo_dir: Path) -> str | None:
    """Return the last CENSO sub-stage with results in *censo_dir*, if any."""
    for section in reversed(CENSO_SECTIONS):
        if (censo_dir / f"{CENSO_PARTS[section]}.json").exists():
            return section
    return None


def extract_best_conformer_xyz(
    censo_conformers_xyz_file: Path,
    best_conformer_id: str,
    output_xyz_file: Path,
) -> None:
    """Extract a single conformer block from a multi-conformer XYZ file.

    Reads the XYZ file produced by CENSO (which may contain many conformers)
    and writes only the block for *best_conformer_id* to *output_xyz_file*.

    Parameters
    ----------
    censo_conformers_xyz_file : Path
        Path to the multi-conformer XYZ file (e.g. ``3_REFINEMENT.xyz`` as produced by CENSO).
    best_conformer_id : str
        Conformer label to search for (e.g. ``"CONF5"``).  The search looks
        for this string in the comment line of each XYZ block.
    output_xyz_file : Path
        Destination path for the extracted single-conformer XYZ file.
    """
    logger.debug(
        f"Extracting best conformer ({best_conformer_id}) .xyz from {censo_conformers_xyz_file}",
    )
    with open(censo_conformers_xyz_file) as fin:
        conformers_xyz: list[str] = fin.readlines()
        n_atoms: int = int(conformers_xyz[0])
        for i, line in enumerate(conformers_xyz):
            if best_conformer_id in line.split():
                with open(output_xyz_file, "w") as fout:
                    fout.writelines(conformers_xyz[i - 1 : i + n_atoms + 1])
                break
        else:
            raise ValueError(
                f"Conformer {best_conformer_id} not found in {censo_conformers_xyz_file}",
            )

    logger.debug(f"Finished extracting best conformer xyz to {output_xyz_file}")
