"""Tests for larest.rdkit conformer generation and xTB ranking."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from larest.rdkit import parse_best_rdkit_conformer, run_rdkit


class TestParseBestRdkitConformer:
    def test_returns_lowest_g_conformer(self, rdkit_results_file):
        result = parse_best_rdkit_conformer(rdkit_results_file)
        # conformer_id=2 has G=-350000.00 (lowest)
        assert result["conformer_id"] == pytest.approx(2.0)

    def test_returns_correct_h(self, rdkit_results_file):
        result = parse_best_rdkit_conformer(rdkit_results_file)
        assert result["H"] == pytest.approx(-324567.89, rel=1e-6)

    def test_returns_correct_g(self, rdkit_results_file):
        result = parse_best_rdkit_conformer(rdkit_results_file)
        assert result["G"] == pytest.approx(-350000.00, rel=1e-6)

    def test_keys_present(self, rdkit_results_file):
        result = parse_best_rdkit_conformer(rdkit_results_file)
        assert set(result.keys()) == {"conformer_id", "H", "S", "G"}

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_best_rdkit_conformer(tmp_path / "missing.csv")


class TestRunRdkitMocked:
    """Tests for run_rdkit with xTB subprocess mocked."""

    def _fake_xtb_output(self, xtb_dir, conformer_id):
        """Write a minimal xTB output file to satisfy parse_xtb_output."""
        out = xtb_dir / f"conformer_{conformer_id}.txt"
        out.write_text(
            "         :: TOTAL ENTHALPY        -0.12345678 Eh           ::\n"
            "         :: TOTAL FREE ENERGY     -0.12345000 Eh           ::\n",
        )
        return out

    def test_run_rdkit_returns_rdkit_key(self, tmp_path, minimal_config):
        with patch("larest.rdkit.subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0)

            # Mock parse_xtb_output to avoid needing real xTB output files
            fake_xtb = {"H": -324567.0, "S": -85.0, "G": -350000.0}
            with patch("larest.rdkit.parse_xtb_output", return_value=fake_xtb):
                result = run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)

        assert "rdkit" in result
        assert "H" in result["rdkit"]
        assert "G" in result["rdkit"]

    def test_run_rdkit_creates_sdf(self, tmp_path, minimal_config):
        with patch("larest.rdkit.subprocess.run"):
            fake_xtb = {"H": -324567.0, "S": -85.0, "G": -350000.0}
            with patch("larest.rdkit.parse_xtb_output", return_value=fake_xtb):
                run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)

        assert (tmp_path / "rdkit" / "conformers.sdf").exists()

    def test_run_rdkit_creates_results_csv(self, tmp_path, minimal_config):
        with patch("larest.rdkit.subprocess.run"):
            fake_xtb = {"H": -324567.0, "S": -85.0, "G": -350000.0}
            with patch("larest.rdkit.parse_xtb_output", return_value=fake_xtb):
                run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)

        assert (tmp_path / "xtb" / "rdkit" / "results.csv").exists()

    def test_run_rdkit_invalid_smiles_raises(self, tmp_path, minimal_config):
        with pytest.raises(
            ValueError,
            match="SMILES Parse Error",
        ):
            run_rdkit("not_a_smiles!!!", tmp_path, minimal_config)

    def test_run_rdkit_raises_when_all_conformers_fail(self, tmp_path, minimal_config):
        import subprocess

        with (
            patch(
                "larest.rdkit.subprocess.run",
                side_effect=subprocess.CalledProcessError(1, "xtb"),
            ),
            pytest.raises(RuntimeError, match="All xTB conformer calculations failed"),
        ):
            run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)


# ---------------------------------------------------------------------------
# Integration tests (require xtb installed)
# ---------------------------------------------------------------------------


@pytest.mark.integration
class TestRunRdkitIntegration:
    def test_run_rdkit_beta_propiolactone(self, tmp_path, minimal_config):
        result = run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)
        assert "rdkit" in result
        assert result["rdkit"]["H"] is not None
        assert result["rdkit"]["G"] is not None
        assert result["rdkit"]["S"] is not None

    def test_run_rdkit_creates_xtb_results(self, tmp_path, minimal_config):
        run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)
        assert (tmp_path / "xtb" / "rdkit" / "results.csv").exists()

    def test_run_rdkit_creates_conformer_xyz(self, tmp_path, minimal_config):
        run_rdkit("C1CC(=O)O1", tmp_path, minimal_config)
        # At least one xyz file should exist for the single conformer
        xyz_files = list((tmp_path / "rdkit").glob("conformer_*.xyz"))
        assert len(xyz_files) >= 1
