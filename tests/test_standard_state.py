"""Tests for larest.standard_state."""

from __future__ import annotations

import math

import pytest

from larest.constants import GAS_CONSTANT
from larest.standard_state import apply_standard_state, get_standard_state_correction

T = 298.15


def _results():
    params = {"E": -100.0, "H": -90.0, "S": 300.0, "G": -1000.0}
    return {
        "rdkit": dict(params),
        "censo_refinement_min": dict(params),
        "crest_entropy": {"S_conf": 1.0, "S_rrho": 0.5, "S_total": 1.5},
    }


def _config(standard_state="auto", alpb="toluene", *, gas_phase=False):
    return {
        "thermo": {"temperature": T, "standard_state": standard_state},
        "xtb": {"alpb": alpb},
        "censo": {"general": {"solvent": "toluene", "gas_phase": gas_phase}},
    }


def test_correction_value():
    # RT ln(24.46) = 7.92 kJ/mol at 298.15 K
    assert get_standard_state_correction(T) == pytest.approx(7924, abs=2)
    assert get_standard_state_correction(T) == pytest.approx(
        GAS_CONSTANT * T * math.log(24.465),
        rel=1e-4,
    )


def test_solution_converted_to_1m():
    corrected = apply_standard_state(_results(), _config())
    dg = get_standard_state_correction(T)
    for section in ("rdkit", "censo_refinement_min"):
        assert corrected[section]["G"] == pytest.approx(-1000.0 + dg)
        assert corrected[section]["S"] == pytest.approx(300.0 - dg / T)
        assert corrected[section]["H"] == -90.0
        assert corrected[section]["E"] == -100.0
    # G = H - TS remains consistent
    g, s = corrected["rdkit"]["G"], corrected["rdkit"]["S"]
    assert g is not None
    assert s is not None
    assert g - (-1000.0) == pytest.approx(-T * (s - 300.0))
    assert corrected["crest_entropy"] == _results()["crest_entropy"]


def test_gas_phase_unchanged():
    corrected = apply_standard_state(_results(), _config(alpb=False, gas_phase=True))
    assert corrected == _results()


def test_per_stage_solvation():
    corrected = apply_standard_state(_results(), _config(gas_phase=True))
    assert corrected["rdkit"]["G"] != -1000.0
    assert corrected["censo_refinement_min"]["G"] == -1000.0


def test_forced_states():
    assert apply_standard_state(_results(), _config("1atm")) == _results()
    forced = apply_standard_state(_results(), _config("1M", alpb=False, gas_phase=True))
    assert forced["rdkit"]["G"] != -1000.0


def test_none_values_kept():
    results = {"censo_prescreening": {"E": -1.0, "H": None, "S": None, "G": None}}
    assert apply_standard_state(results, _config()) == results


def test_invalid_state_raises():
    with pytest.raises(ValueError, match="standard_state"):
        apply_standard_state(_results(), _config("1bar"))
