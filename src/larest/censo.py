"""CENSO DFT refinement stage of the LaREST pipeline.

Wraps the CENSO (Conformer Ensemble Sorting) program, which performs four
sequential sub-stages of increasing DFT accuracy:

0. ``censo_prescreening`` — fast pre-screening of the CREST ensemble
1. ``censo_screening``    — singlepoint DFT screening
2. ``censo_optimization`` — DFT geometry optimisation
3. ``censo_refinement``   — high-accuracy single-point refinement

Thermodynamic parameters (H, G) are extracted from the CENSO output at each
sub-stage, and entropy is derived as S = (H - G) / T.
"""

from __future__ import annotations

import csv
import json
import logging
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

from larest.constants import CENSO_SECTIONS, HARTTREE_TO_JMOL, THERMODYNAMIC_PARAMS
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

    censo_results: dict[str, dict[str, float | None]] = parse_censo_output(
        censo_output_file=censo_output_file,
        temperature=config["censo"]["general"].get("temperature", 298.15),
    )

    censo_results_file = censo_dir / "results.json"
    logger.debug(f"Writing results to {censo_results_file}")
    with open(censo_results_file, "w") as fstream:
        json.dump(censo_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return censo_results


def parse_censo_output(
    censo_output_file: Path,
    temperature: float,
) -> dict[str, dict[str, float | None]]:
    """Extract per-sub-stage thermodynamic parameters from a CENSO output file.

    Scans for ``"part0"``, ``"part1"``, ``"part2"``, ``"part3"`` markers and
    reads Boltzmann-averaged H and G values from each.  Entropy is derived as
    ``S = (H - G) / T``.  Values are converted from Hartree to J/mol.

    Parameters
    ----------
    censo_output_file : Path
        Path to the plain-text CENSO output file.
    temperature : float
        Temperature in Kelvin used to derive entropy from H and G.

    Returns
    -------
    dict[str, dict[str, float | None]]
        Nested dict mapping each CENSO section (``"censo_prescreening"``,
        ``"censo_screening"``, ``"censo_optimization"``, ``"censo_refinement"``) to a
        thermodynamic parameter dict with keys ``"H"``, ``"S"``, ``"G"``
        (J/mol).  A value is ``None`` if it could not be extracted.
    """
    # WARN: cannot use fromkeys, otherwise they all point to the same mutable dict
    censo_output: dict[str, dict[str, float | None]] = {
        section: dict.fromkeys(THERMODYNAMIC_PARAMS, None) for section in CENSO_SECTIONS
    }

    logger.debug(f"Searching for CENSO results in file {censo_output_file}")
    with open(censo_output_file) as fstream:
        section_no: int = 0
        for i, line in enumerate(fstream):
            if f"part{section_no}" in line:
                try:
                    censo_output[CENSO_SECTIONS[section_no]]["H"] = (
                        float(line.split()[1]) * HARTTREE_TO_JMOL
                    )
                    censo_output[CENSO_SECTIONS[section_no]]["G"] = (
                        float(line.split()[2]) * HARTTREE_TO_JMOL
                    )
                except Exception:
                    logger.exception(
                        f"Failed to extract H and G from line {i}: {line}",
                    )
                else:
                    section_no += 1

    for params in censo_output.values():
        if params["H"] is not None and params["G"] is not None:
            params["S"] = (params["H"] - params["G"]) / temperature
        if not all(v is not None for v in params.values()):
            logger.warning(
                f"Failed to extract necessary data from {censo_output_file}",
            )
            logger.warning("Missing data will be assigned None")
        logger.debug(
            f"Found enthalpy: {params['H']}, free energy: {params['G']}, entropy {params['S']}",
        )

    return censo_output


def parse_best_censo_conformers(
    censo_output_file: Path,
) -> dict[str, str]:
    """Identify the highest-ranked conformer in each CENSO sub-stage.

    Scans the CENSO output for ``"Highest ranked conformer"`` lines and
    records the associated conformer ID (e.g. ``"CONF5"``) for each sub-stage.
    Missing entries fall back to ``"CONF0"`` with a logged warning.

    Parameters
    ----------
    censo_output_file : Path
        Path to the plain-text CENSO output file.

    Returns
    -------
    dict[str, str]
        Mapping of CENSO section name to the ID of the top-ranked conformer
        in that section, e.g. ``{"censo_refinement": "CONF5", ...}``.
    """
    best_censo_conformers: dict[str, str] = {}

    logger.debug(f"Searching for results in file {censo_output_file}")
    section_no: int = 0
    current_conf: str | None = None
    in_results_table: bool = False

    with open(censo_output_file) as fstream:
        for i, line in enumerate(fstream):
            # Each stage's final results table is headed by a dashed line containing
            # "RESULTS".  The rows are sorted by ΔGtot ascending, so the first CONF
            # row is the best conformer.  Reset on each new header so intermediate
            # tables (e.g. "Converged or removed conformers") don't interfere.
            if "--" in line and "RESULTS" in line:
                in_results_table = True
                current_conf = None
            elif in_results_table and current_conf is None:
                token = line.split()[0] if line.split() else ""
                if token.startswith("CONF") and token[4:].isdigit():
                    current_conf = token
            elif f"<<==part{section_no}==" in line and section_no < len(CENSO_SECTIONS):
                section = CENSO_SECTIONS[section_no]
                if current_conf:
                    best_censo_conformers[section] = current_conf
                else:
                    logger.warning(
                        f"Could not find conformer ID for {section} at line {i}",
                    )
                in_results_table = False
                current_conf = None
                section_no += 1

    if len(best_censo_conformers) < len(CENSO_SECTIONS):
        logger.warning(
            f"Failed to extract best conformers from {censo_output_file}",
        )
        logger.warning("Missing sections will be assigned first conformer (CONF1)")
        for section in CENSO_SECTIONS:
            best_censo_conformers.setdefault(section, "CONF1")

    for section in CENSO_SECTIONS:
        logger.debug(
            f"Best conformer in {section}: {best_censo_conformers[section]}",
        )

    return best_censo_conformers


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
