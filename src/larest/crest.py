"""CREST conformer generation and entropy calculation stages of the LaREST pipeline.

This module wraps two separate CREST execution modes:

* :func:`run_crest_confgen` — conformer/rotamer ensemble generation starting
  from the best RDKit/xTB conformer, followed by optional xTB re-ranking of
  the CREST best structure.
* :func:`run_crest_entropy` — entropy-mode run using the best CENSO conformer
  to obtain conformational entropy corrections (S_conf, S_rrho, S_total).
"""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from larest.censo import extract_best_conformer_xyz, parse_best_censo_conformers
from larest.constants import CALMOL_TO_JMOL, CREST_ENTROPY_OUTPUT_PARAMS
from larest.output import create_dir
from larest.rdkit import parse_best_rdkit_conformer
from larest.setup import parse_command_args
from larest.xtb import run_xtb

logger = logging.getLogger(__name__)


def run_crest_confgen(
    dir_path: Path,
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    """Run CREST conformer generation from the best RDKit conformer.

    Reads the best RDKit/xTB conformer, runs CREST to explore the conformer
    and rotamer ensemble, then optionally re-ranks the CREST best structure
    with xTB (controlled by ``config["steps"]["xtb"]``).

    Parameters
    ----------
    dir_path : Path
        Root molecule output directory.  Expected to contain
        ``xtb/rdkit/results.csv`` and the corresponding optimised ``.xyz``
        files from the RDKit stage.
    config : dict[str, Any]
        Full pipeline configuration dict.  Uses the ``[crest][confgen]`` and
        ``[xtb]`` sub-sections.

    Returns
    -------
    dict[str, dict[str, float | None]]
        ``{"crest": {"H": ..., "S": ..., "G": ...}}`` if xTB re-ranking is
        enabled, otherwise an empty dict.

    Raises
    ------
    ValueError
        If the best RDKit conformer ID cannot be determined from the results
        CSV.
    subprocess.CalledProcessError
        If the CREST process exits with a non-zero return code.
    """
    crest_dir = dir_path / "crest_confgen"
    create_dir(crest_dir)

    xtb_rdkit_results_file = dir_path / "xtb" / "rdkit" / "results.csv"
    best_rdkit_conformer_results = parse_best_rdkit_conformer(xtb_rdkit_results_file)
    if best_rdkit_conformer_results["conformer_id"] is None:
        raise ValueError("Failed to extract best RDKit conformer for crest_confgen")

    best_rdkit_conformer_id = int(best_rdkit_conformer_results["conformer_id"])
    best_rdkit_conformer_xyz_file = (
        dir_path
        / "xtb"
        / "rdkit"
        / f"conformer_{best_rdkit_conformer_id}"
        / f"conformer_{best_rdkit_conformer_id}.xtbopt.xyz"
    )

    crest_output_file = crest_dir / "crest.txt"
    crest_args: list[str] = [
        "crest",
        str(best_rdkit_conformer_xyz_file.absolute()),
        *parse_command_args(sub_config=["crest", "confgen"], config=config),
    ]

    with open(crest_output_file, "w") as fstream:
        subprocess.run(
            crest_args,
            stdout=fstream,
            stderr=subprocess.STDOUT,
            cwd=crest_dir,
            check=True,
        )

    if config["steps"]["xtb"]:
        best_crest_conformer_xyz_file = crest_dir / "crest_best.xyz"
        xtb_dir = dir_path / "xtb" / "crest"
        create_dir(xtb_dir)
        xtb_results = run_xtb(
            xtb_input_file=best_crest_conformer_xyz_file,
            xtb_dir=xtb_dir,
            config=config,
        )
        return {"crest": xtb_results}

    return {}


def run_crest_entropy(
    dir_path: Path,
    config: dict[str, Any],
) -> dict[str, float | None]:
    """Run CREST in entropy mode using the best CENSO conformer.

    Extracts the highest-ranked CENSO ``censo_refinement`` conformer to an XYZ
    file, runs CREST entropy mode on it, parses the output, writes
    ``results.json``, and returns the conformational entropy values.

    Parameters
    ----------
    dir_path : Path
        Root molecule output directory.  Must contain ``censo/censo.txt`` and
        ``censo/3_REFINEMENT.xyz`` (as produced by CENSO) from a completed CENSO stage.
    config : dict[str, Any]
        Full pipeline configuration dict.  Uses the ``[crest][entropy]``
        sub-section.

    Returns
    -------
    dict[str, float | None]
        Entropy values keyed by ``"S_conf"``, ``"S_rrho"``, and ``"S_total"``
        (all in J/mol/K).  Values are ``None`` if they could not be extracted.

    Raises
    ------
    ValueError
        If the best CENSO conformer cannot be identified or extracted.
    subprocess.CalledProcessError
        If the CREST process exits with a non-zero return code.
    """
    censo_dir = dir_path / "censo"
    best_censo_conformer_xyz_file = censo_dir / "censo_best.xyz"

    try:
        best_censo_conformer = parse_best_censo_conformers(
            censo_output_file=censo_dir / "censo.txt",
        )["censo_refinement"]
        extract_best_conformer_xyz(
            censo_conformers_xyz_file=censo_dir / "3_REFINEMENT.xyz",
            best_conformer_id=best_censo_conformer,
            output_xyz_file=best_censo_conformer_xyz_file,
        )
    except Exception as err:
        raise ValueError(
            "Failed to extract best CENSO conformer for crest_entropy",
        ) from err

    crest_dir = dir_path / "crest_entropy"
    create_dir(crest_dir)

    crest_output_file = crest_dir / "crest.txt"
    crest_args: list[str] = [
        "crest",
        str(best_censo_conformer_xyz_file.absolute()),
        "--entropy",
        *parse_command_args(sub_config=["crest", "entropy"], config=config),
    ]

    with open(crest_output_file, "w") as fstream:
        subprocess.run(
            crest_args,
            stdout=fstream,
            stderr=subprocess.STDOUT,
            cwd=crest_dir,
            check=True,
        )

    crest_results: dict[str, float | None] = parse_crest_entropy_output(
        crest_output_file=crest_output_file,
    )

    crest_results_file = crest_dir / "results.json"
    logger.debug(f"Writing results to {crest_results_file}")
    with open(crest_results_file, "w") as fstream:
        json.dump(crest_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return crest_results


def parse_crest_entropy_output(
    crest_output_file: Path,
) -> dict[str, float | None]:
    """Extract conformational entropy values from a CREST entropy output file.

    Scans for the ``"Sconf"``, ``"δSrrho"``, and ``"S(total)"`` markers and
    converts the extracted values from cal/mol/K to J/mol/K.

    Parameters
    ----------
    crest_output_file : Path
        Path to the plain-text CREST entropy output file.

    Returns
    -------
    dict[str, float | None]
        Dict with keys ``"S_conf"``, ``"S_rrho"``, and ``"S_total"``
        (all in J/mol/K).  Any value that could not be extracted is ``None``;
        a warning is logged in that case.
    """
    crest_output: dict[str, float | None] = dict.fromkeys(
        CREST_ENTROPY_OUTPUT_PARAMS,
        None,
    )

    logger.debug(f"Searching for results in file {crest_output_file}")
    with open(crest_output_file) as fstream:
        for i, line in enumerate(fstream):
            if "Sconf" in line:
                try:
                    crest_output["S_conf"] = float(line.split()[-1]) * CALMOL_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract S_conf from line {i}: {line}",
                    )
            elif ("+" in line) and ("δSrrho" in line) and (len(line.split()) == 4):
                try:
                    crest_output["S_rrho"] = float(line.split()[-1]) * CALMOL_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract S_rrho from line {i}: {line}",
                    )
            elif ("S(total)" in line) and ("cal" in line):
                try:
                    crest_output["S_total"] = float(line.split()[3]) * CALMOL_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract S_total from line {i}: {line}",
                    )

    if not all(value is not None for value in crest_output.values()):
        logger.warning(f"Failed to extract necessary data from {crest_output_file}")
        logger.warning("Missing data will be assigned None")
    logger.debug(
        f"Found S_conf: {crest_output['S_conf']}, S_rrho: {crest_output['S_rrho']}",
    )
    logger.debug(f"S_total {crest_output['S_total']}")

    return crest_output
