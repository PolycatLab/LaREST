"""Tests for larest.xtb parsing and execution."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from larest.constants import HARTTREE_TO_JMOL
from larest.xtb import parse_xtb_output, run_xtb

TEMPERATURE = 298.15

# Expected parsed values from tests/data/xtb_output.txt
_H_HARTREE = -0.123456789
_G_HARTREE = -0.123456000
_H_EXPECTED = _H_HARTREE * HARTTREE_TO_JMOL
_G_EXPECTED = _G_HARTREE * HARTTREE_TO_JMOL
_S_EXPECTED = (_H_EXPECTED - _G_EXPECTED) / TEMPERATURE


class TestParseXtbOutput:
    def test_parses_enthalpy(self, xtb_output_file):
        result = parse_xtb_output(xtb_output_file, TEMPERATURE)
        assert result["H"] == pytest.approx(_H_EXPECTED, rel=1e-6)

    def test_parses_free_energy(self, xtb_output_file):
        result = parse_xtb_output(xtb_output_file, TEMPERATURE)
        assert result["G"] == pytest.approx(_G_EXPECTED, rel=1e-6)

    def test_computes_entropy(self, xtb_output_file):
        result = parse_xtb_output(xtb_output_file, TEMPERATURE)
        assert result["S"] == pytest.approx(_S_EXPECTED, rel=1e-6)

    def test_returns_none_for_missing_values(self, tmp_path):
        missing_file = tmp_path / "empty.txt"
        missing_file.write_text("no useful data here\n")
        result = parse_xtb_output(missing_file, TEMPERATURE)
        assert result["H"] is None
        assert result["G"] is None
        assert result["S"] is None

    def test_partial_missing(self, data_dir):
        result = parse_xtb_output(data_dir / "xtb_output_missing.txt", TEMPERATURE)
        assert result["H"] is None
        assert result["G"] is None
        assert result["S"] is None

    def test_keys_present(self, xtb_output_file):
        result = parse_xtb_output(xtb_output_file, TEMPERATURE)
        assert set(result.keys()) == {"H", "S", "G"}


class TestRunXtbMocked:
    """Tests for run_xtb with subprocess mocked out."""

    def _make_xtb_output(
        self,
        tmp_path: Path,
        h_hartree: float,
        g_hartree: float,
    ) -> Path:
        output_file = tmp_path / "xtb.txt"
        output_file.write_text(
            f"         :: TOTAL ENTHALPY        {h_hartree} Eh           ::\n"
            f"         :: TOTAL FREE ENERGY     {g_hartree} Eh           ::\n",
        )
        return output_file

    def test_run_xtb_returns_results(self, tmp_path, minimal_config):
        xyz_file = tmp_path / "mol.xyz"
        xyz_file.write_text("1\ntest\nC 0.0 0.0 0.0\n")

        with patch("larest.xtb.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            # Pre-write the expected output so parse_xtb_output finds it
            self._make_xtb_output(tmp_path, -0.12345678, -0.12345000)

            result = run_xtb(xyz_file, tmp_path, minimal_config)

        assert "H" in result
        assert "S" in result
        assert "G" in result

    def test_run_xtb_writes_results_json(self, tmp_path, minimal_config):
        xyz_file = tmp_path / "mol.xyz"
        xyz_file.write_text("1\ntest\nC 0.0 0.0 0.0\n")

        with patch("larest.xtb.subprocess.run"):
            self._make_xtb_output(tmp_path, -0.12345678, -0.12345000)
            run_xtb(xyz_file, tmp_path, minimal_config)

        results_json = tmp_path / "results.json"
        assert results_json.exists()
        data = json.loads(results_json.read_text())
        assert "H" in data

    def test_run_xtb_passes_config_args(self, tmp_path, minimal_config):
        xyz_file = tmp_path / "mol.xyz"
        xyz_file.write_text("1\ntest\nC 0.0 0.0 0.0\n")

        with patch("larest.xtb.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            self._make_xtb_output(tmp_path, -0.12345678, -0.12345000)
            run_xtb(xyz_file, tmp_path, minimal_config)

        call_args = mock_run.call_args[0][0]
        assert "xtb" in call_args
        assert "--gfn" in call_args

    def test_run_xtb_raises_on_nonzero_exit(self, tmp_path, minimal_config):
        import subprocess

        xyz_file = tmp_path / "mol.xyz"
        xyz_file.write_text("1\ntest\nC 0.0 0.0 0.0\n")

        with patch("larest.xtb.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, "xtb")
            with pytest.raises(subprocess.CalledProcessError):
                run_xtb(xyz_file, tmp_path, minimal_config)


# ---------------------------------------------------------------------------
# Integration tests (require xtb installed)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestRunXtbIntegration:
    """Runs xtb on a tiny molecule to verify end-to-end parsing."""

    BETA_PL_XYZ = """\
12
beta-propiolactone
C   0.000000   0.000000   0.000000
C   1.540000   0.000000   0.000000
C   2.100000   1.200000   0.000000
O   1.300000   2.200000   0.000000
O   0.000000   1.400000   0.000000
H  -0.380000  -1.020000   0.000000
H  -0.380000   0.510000   0.883000
H   1.920000  -0.510000   0.883000
H   1.920000  -0.510000  -0.883000
H   2.880000   1.300000   0.882000
H   2.880000   1.300000  -0.882000
H  -0.380000   0.510000  -0.883000
"""

    def test_xtb_returns_finite_values(self, tmp_path, minimal_config):
        xyz_file = tmp_path / "mol.xyz"
        xyz_file.write_text(self.BETA_PL_XYZ)

        result = run_xtb(xyz_file, tmp_path, minimal_config)

        h = result["H"]
        g = result["G"]
        s = result["S"]
        assert h is not None
        assert g is not None
        assert s is not None
        assert g < 0
