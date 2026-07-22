"""xTB execution and output parsing utilities.

Provides :func:`run_xtb` to invoke the ``xtb`` binary on a single geometry
file and :func:`parse_xtb_output` to extract enthalpy (H), entropy (S), and
free energy (G) from the resulting plain-text output.  Entropy is derived from
H and G using the relation S = (H - G) / T.
"""

from __future__ import annotations

import json
import logging
import subprocess
from pathlib import Path
from typing import Any

from larest.constants import HARTTREE_TO_JMOL, THERMODYNAMIC_PARAMS
from larest.setup import parse_command_args

logger = logging.getLogger(__name__)


def write_thermo_input(input_dir: Path, temperature: float) -> Path:
    """Write an xTB ``$thermo`` xcontrol file fixing the thermostatistical temperature.

    xTB evaluates its thermostatistics — and hence the ``TOTAL FREE ENERGY`` it
    reports — at 298.15 K unless told otherwise.  To keep that free energy (and
    the entropy derived from it as ``S = (H - G) / T``) consistent with the
    pipeline's canonical temperature, we pass a ``$thermo`` block via
    ``--input``.  This is distinct from ``--etemp``, the electronic
    (Fermi-smearing) temperature, which does not set the thermostatistical T.

    Parameters
    ----------
    input_dir : Path
        Directory in which the ``xtb.inp`` xcontrol file is written.
    temperature : float
        Thermostatistical temperature in Kelvin.

    Returns
    -------
    Path
        Path to the written xcontrol file.
    """
    thermo_input_file = input_dir / "xtb.inp"
    thermo_input_file.write_text(f"$thermo\n   temp={temperature}\n$end\n")
    logger.debug(f"Wrote xTB thermo input (temp={temperature}) to {thermo_input_file}")
    return thermo_input_file


def run_xtb(
    xtb_input_file: Path,
    xtb_dir: Path,
    config: dict[str, Any],
) -> dict[str, float | None]:
    """Run xTB on a single geometry file and return thermodynamic parameters.

    Constructs the xTB command from *config*, executes it via
    :func:`subprocess.run`, parses the plain-text output, writes a
    ``results.json`` file, and returns the parsed values.

    Parameters
    ----------
    xtb_input_file : Path
        Path to the input geometry file (e.g. an ``.xyz`` file).
    xtb_dir : Path
        Working directory in which xTB is executed and output files are
        written.
    config : dict[str, Any]
        Full pipeline configuration dict.  The ``[xtb]`` section is used to
        build the CLI flags via :func:`~larest.setup.parse_command_args`.

    Returns
    -------
    dict[str, float | None]
        Thermodynamic parameters keyed by ``"H"``, ``"S"``, and ``"G"``
        (all in J/mol).  Values are ``None`` if they could not be extracted.

    Raises
    ------
    subprocess.CalledProcessError
        If the xTB process exits with a non-zero return code.
    """
    temperature: float = config["xtb"].get("temperature", 298.15)
    thermo_input_file = write_thermo_input(xtb_dir, temperature)

    xtb_args: list[str] = [
        "xtb",
        str(xtb_input_file.absolute()),
        "--namespace",
        xtb_input_file.name.split(".")[0],
        "--input",
        str(thermo_input_file.absolute()),
        *parse_command_args(sub_config=["xtb"], config=config),
    ]

    xtb_output_file = xtb_dir / "xtb.txt"
    with open(xtb_output_file, "w") as fstream:
        subprocess.run(
            xtb_args,
            stdout=fstream,
            stderr=subprocess.STDOUT,
            cwd=xtb_dir,
            check=True,
        )

    xtb_results: dict[str, float | None] = parse_xtb_output(
        xtb_output_file=xtb_output_file,
        temperature=temperature,
    )

    xtb_results_file = xtb_dir / "results.json"
    logger.debug(f"Writing results to {xtb_results_file}")
    with open(xtb_results_file, "w") as fstream:
        json.dump(xtb_results, fstream, sort_keys=True, allow_nan=True, indent=4)

    return xtb_results


def parse_xtb_output(
    xtb_output_file: Path,
    temperature: float,
) -> dict[str, float | None]:
    """Extract H, S, and G from an xTB plain-text output file.

    Scans for ``"TOTAL ENTHALPY"`` and ``"TOTAL FREE ENERGY"`` markers.
    Values are converted from Hartree to J/mol.  Entropy is computed as
    ``S = (H - G) / T``.

    Parameters
    ----------
    xtb_output_file : Path
        Path to the xTB output text file to parse.
    temperature : float
        Thermostatistical temperature in Kelvin (the T at which xTB evaluated
        the free energy), used to derive entropy from H and G.

    Returns
    -------
    dict[str, float | None]
        Dict with keys ``"H"``, ``"S"``, ``"G"`` (all in J/mol).  Any value
        that could not be extracted is ``None``; a warning is logged in that
        case.
    """
    xtb_output: dict[str, float | None] = dict.fromkeys(THERMODYNAMIC_PARAMS, None)

    logger.debug(f"Searching for xTB results in file {xtb_output_file}")
    with open(xtb_output_file) as fstream:
        for i, line in enumerate(fstream):
            if "TOTAL ENTHALPY" in line:
                try:
                    xtb_output["H"] = float(line.split()[3]) * HARTTREE_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract total enthalpy from line {i}: {line}",
                    )
            elif "TOTAL FREE ENERGY" in line:
                try:
                    xtb_output["G"] = float(line.split()[4]) * HARTTREE_TO_JMOL
                except Exception:
                    logger.exception(
                        f"Failed to extract total free energy from line {i}: {line}",
                    )

    if xtb_output["H"] is not None and xtb_output["G"] is not None:
        xtb_output["S"] = (xtb_output["H"] - xtb_output["G"]) / temperature
    if not all(v is not None for v in xtb_output.values()):
        logger.warning(f"Failed to extract necessary data from {xtb_output_file}")
        logger.warning("Missing data will be assigned None")
    logger.debug(
        f"Found enthalpy: {xtb_output['H']}, entropy: {xtb_output['S']}, free energy {xtb_output['G']}",
    )

    return xtb_output
