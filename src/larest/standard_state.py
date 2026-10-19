"""Standard-state conversion of free energies and entropies.

The thermostatistical (mRRHO) contributions from xTB, and hence from CENSO,
are those of an ideal gas at 1 atm, whereas implicit solvation free energies
(ALPB in xTB, SMD/CPCM in CENSO) refer to 1 M in the gas phase and 1 M in
solution.  Solution-phase free energies are therefore converted to a 1 M
standard state by adding ``RT ln(c° RT / p°)`` per molecule (7.9 kJ/mol at
298.15 K), with the corresponding entropy change ``-R ln(c° RT / p°)``.  This
matters for reaction energies whenever the number of molecules changes, as in
ring-opening polymerisation, and makes ΔG and ΔS comparable with experimental
values for a 1 M standard state.  E and H are unaffected.
"""

from __future__ import annotations

import logging
import math
from typing import Any

from larest.constants import GAS_CONSTANT

logger = logging.getLogger(__name__)

STANDARD_PRESSURE = 101325.0  # Pa (1 atm)
STANDARD_CONCENTRATION = 1000.0  # mol/m^3 (1 M)
STANDARD_STATES: frozenset[str] = frozenset({"auto", "1M", "1atm"})


def get_standard_state_correction(temperature: float) -> float:
    """Return the free energy change from a 1 atm to a 1 M standard state.

    Parameters
    ----------
    temperature : float
        Temperature in Kelvin.

    Returns
    -------
    float
        ``RT ln(c° RT / p°)`` in J/mol.
    """
    ratio = STANDARD_CONCENTRATION * GAS_CONSTANT * temperature / STANDARD_PRESSURE
    return GAS_CONSTANT * temperature * math.log(ratio)


def _xtb_solvated(config: dict[str, Any]) -> bool:
    xtb = config.get("xtb", {})
    return any(xtb.get(key) for key in ("alpb", "gbsa", "cpcmx"))


def _censo_solvated(config: dict[str, Any]) -> bool:
    general = config.get("censo", {}).get("general", {})
    return bool(
        general.get("solvent")
        and not (general.get("gas_phase") or general.get("gas-phase")),
    )


def apply_standard_state(
    results: dict[str, dict[str, float | None]],
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    """Convert G and S of every section to the configured standard state.

    ``[thermo].standard_state`` selects the standard state: ``"auto"`` (the
    default) uses 1 M for sections computed with implicit solvation and 1 atm
    for gas-phase sections; ``"1M"`` and ``"1atm"`` force either.  Sections from
    xTB (``rdkit``, ``crest``) follow the xTB solvation settings and CENSO
    sections follow the CENSO ones.

    Parameters
    ----------
    results : dict[str, dict[str, float | None]]
        Results of one molecule, as computed (1 atm ideal-gas thermostatistics).
    config : dict[str, Any]
        Full pipeline configuration dict.

    Returns
    -------
    dict[str, dict[str, float | None]]
        New results dict with converted ``G`` and ``S``.
    """
    thermo = config.get("thermo", {})
    standard_state = thermo.get("standard_state", "auto")
    if standard_state not in STANDARD_STATES:
        raise ValueError(
            f"Unknown [thermo].standard_state {standard_state!r}; "
            f"choose from {', '.join(sorted(STANDARD_STATES))}",
        )
    temperature = thermo.get("temperature", 298.15)
    correction = get_standard_state_correction(temperature)

    converted: dict[str, dict[str, float | None]] = {}
    for section, params in results.items():
        if section.startswith("censo"):
            solvated = _censo_solvated(config)
        elif section in ("rdkit", "crest"):
            solvated = _xtb_solvated(config)
        else:
            converted[section] = params
            continue
        to_1m = standard_state == "1M" or (standard_state == "auto" and solvated)
        new_params = dict(params)
        g, s = new_params.get("G"), new_params.get("S")
        if to_1m and g is not None:
            new_params["G"] = g + correction
        if to_1m and s is not None:
            new_params["S"] = s - correction / temperature
        converted[section] = new_params
    return converted
