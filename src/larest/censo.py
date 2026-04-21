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

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from larest.constants import CENSO_SECTIONS, HARTTREE_TO_JMOL, THERMODYNAMIC_PARAMS
from larest.output import create_dir
from larest.setup import parse_command_args

logger = logging.getLogger(__name__)


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
            if header == "cli":
                continue
            fstream.write(f"[{header}]\n")
            fstream.writelines(
                f"{key} = {value}\n" for key, value in sub_config.items()
            )
            fstream.write("\n")
    logger.debug(f"Created censo config file at {censorc_file}")


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
        Root molecule output directory.  Must contain
        ``crest_confgen/crest_conformers.xyz`` from a completed CREST stage.
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
    crest_conformers_file = dir_path / "crest_confgen" / "crest_conformers.xyz"

    censo_args: list[str] = [
        "censo",
        "--input",
        str(crest_conformers_file.absolute()),
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
    best_censo_conformers: dict[str, str] = dict.fromkeys(CENSO_SECTIONS, "CONF0")

    logger.debug(f"Searching for results in file {censo_output_file}")
    section_no: int = 0
    with open(censo_output_file) as fstream:
        for i, line in enumerate(fstream):
            if "Highest ranked conformer" in line:
                try:
                    best_censo_conformers[CENSO_SECTIONS[section_no]] = line.split()[-1]
                except Exception:
                    logger.exception(
                        f"Failed to extract best conformer from line {i}: {line}",
                    )
                else:
                    section_no += 1

    if section_no < len(CENSO_SECTIONS):
        logger.warning(
            f"Failed to extract best conformers from {censo_output_file}",
        )
        logger.warning("Missing sections will be assigned first conformer (CONF0)")
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
