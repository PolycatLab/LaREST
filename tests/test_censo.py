"""Tests for larest.censo output parsing utilities (no ORCA/CENSO required)."""

from __future__ import annotations

import pytest

from larest.censo import (
    build_rdkit_ensemble_xyz,
    create_censorc,
    extract_best_conformer_xyz,
    parse_best_censo_conformers,
    parse_censo_output,
)
from larest.constants import CENSO_SECTIONS, HARTTREE_TO_JMOL

TEMPERATURE = 298.15

# Expected values from tests/data/censo_output.txt
_SECTION_HARTREES = [
    (-0.12345678, -0.12345000),
    (-0.23456789, -0.23456000),
    (-0.34567890, -0.34567000),
    (-0.45678901, -0.45678000),
]


class TestParseCensoOutput:
    def test_parses_all_sections(self, censo_output_file):
        result = parse_censo_output(censo_output_file, TEMPERATURE)
        assert set(result.keys()) == set(CENSO_SECTIONS)

    def test_parses_enthalpy(self, censo_output_file):
        result = parse_censo_output(censo_output_file, TEMPERATURE)
        for section, (h_ha, _) in zip(CENSO_SECTIONS, _SECTION_HARTREES, strict=False):
            assert result[section]["H"] == pytest.approx(
                h_ha * HARTTREE_TO_JMOL,
                rel=1e-6,
            )

    def test_parses_free_energy(self, censo_output_file):
        result = parse_censo_output(censo_output_file, TEMPERATURE)
        for section, (_, g_ha) in zip(CENSO_SECTIONS, _SECTION_HARTREES, strict=False):
            assert result[section]["G"] == pytest.approx(
                g_ha * HARTTREE_TO_JMOL,
                rel=1e-6,
            )

    def test_computes_entropy(self, censo_output_file):
        result = parse_censo_output(censo_output_file, TEMPERATURE)
        for section, (h_ha, g_ha) in zip(
            CENSO_SECTIONS,
            _SECTION_HARTREES,
            strict=False,
        ):
            expected_s = (h_ha - g_ha) * HARTTREE_TO_JMOL / TEMPERATURE
            assert result[section]["S"] == pytest.approx(expected_s, rel=1e-5)

    def test_missing_data_returns_none(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("nothing here\n")
        result = parse_censo_output(empty, TEMPERATURE)
        for section in CENSO_SECTIONS:
            assert result[section]["H"] is None
            assert result[section]["G"] is None


class TestParseBestCensoConformers:
    def test_parses_all_sections(self, censo_output_file):
        result = parse_best_censo_conformers(censo_output_file)
        assert set(result.keys()) == set(CENSO_SECTIONS)

    def test_correct_conformer_ids(self, censo_output_file):
        result = parse_best_censo_conformers(censo_output_file)
        assert result["censo_prescreening"] == "CONF1"
        assert result["censo_screening"] == "CONF2"
        assert result["censo_optimization"] == "CONF3"
        assert result["censo_refinement"] == "CONF5"

    def test_missing_data_defaults_to_conf1(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("nothing here\n")
        result = parse_best_censo_conformers(empty)
        for section in CENSO_SECTIONS:
            assert result[section] == "CONF1"

    def test_missing_sections_logs_warning(self, tmp_path, caplog):
        import logging

        # File contains fewer "Highest ranked conformer" lines than CENSO_SECTIONS
        partial = tmp_path / "partial.txt"
        partial.write_text(
            "Highest ranked conformer CONF1\nHighest ranked conformer CONF2\n",
        )
        with caplog.at_level(logging.WARNING):
            parse_best_censo_conformers(partial)
        assert "Failed to extract best conformers" in caplog.text


class TestExtractBestConformerXyz:
    def test_extracts_correct_conformer(self, censo_conformers_xyz_file, tmp_path):
        out = tmp_path / "best.xyz"
        extract_best_conformer_xyz(censo_conformers_xyz_file, "CONF5", out)
        assert out.exists()
        content = out.read_text()
        assert "CONF5" in content
        assert "CONF1" not in content

    def test_extracted_xyz_has_correct_atom_count(
        self,
        censo_conformers_xyz_file,
        tmp_path,
    ):
        out = tmp_path / "best.xyz"
        extract_best_conformer_xyz(censo_conformers_xyz_file, "CONF1", out)
        lines = out.read_text().splitlines()
        n_atoms = int(lines[0])
        assert n_atoms == 3
        # Should have header + comment + n_atoms lines
        assert len(lines) == n_atoms + 2

    def test_extracts_first_conformer(self, censo_conformers_xyz_file, tmp_path):
        out = tmp_path / "first.xyz"
        extract_best_conformer_xyz(censo_conformers_xyz_file, "CONF1", out)
        assert "CONF1" in out.read_text()

    def test_missing_conformer_id_raises(self, censo_conformers_xyz_file, tmp_path):
        out = tmp_path / "best.xyz"
        with pytest.raises(ValueError, match="CONF99"):
            extract_best_conformer_xyz(censo_conformers_xyz_file, "CONF99", out)


class TestCreateCensorc:
    def test_creates_censorc_file(self, tmp_path, minimal_config):
        create_censorc(minimal_config, tmp_path)
        censorc = tmp_path / ".censo2rc"
        assert censorc.exists()

    def test_censorc_contains_sections(self, tmp_path, minimal_config):
        create_censorc(minimal_config, tmp_path)
        content = (tmp_path / ".censo2rc").read_text()
        assert "[general]" in content
        assert "[cli]" not in content

    def test_censorc_contains_temperature(self, tmp_path, minimal_config):
        create_censorc(minimal_config, tmp_path)
        content = (tmp_path / ".censo2rc").read_text()
        assert "298.15" in content

    def test_censorc_normalises_hyphenated_keys(self, tmp_path, minimal_config):
        # CENSO field names use underscores; a hyphenated key like `gas-phase`
        # must be written as `gas_phase` or CENSO silently ignores it (leaving
        # solvation enabled).
        minimal_config["censo"]["general"]["gas-phase"] = True
        create_censorc(minimal_config, tmp_path)
        content = (tmp_path / ".censo2rc").read_text()
        assert "gas_phase = True" in content
        assert "gas-phase" not in content


class TestBuildRdkitEnsembleXyz:
    @staticmethod
    def _write_conformer(dir_path, conformer_id, comment):
        conf_dir = dir_path / "xtb" / "rdkit" / f"conformer_{conformer_id}"
        conf_dir.mkdir(parents=True, exist_ok=True)
        xyz = conf_dir / f"conformer_{conformer_id}.xtbopt.xyz"
        xyz.write_text(f"1\n{comment}\nO 0.0 0.0 0.0\n")

    @staticmethod
    def _write_results(dir_path, rows):
        results_dir = dir_path / "xtb" / "rdkit"
        results_dir.mkdir(parents=True, exist_ok=True)
        lines = ["conformer_id,H,S,G"]
        lines += [f"{cid},0.0,0.0,{g}" for cid, g in rows]
        (results_dir / "results.csv").write_text("\n".join(lines) + "\n")

    def test_orders_conformers_by_ascending_g(self, tmp_path):
        # conformer 0 has higher G than conformer 1, so conformer 1 comes first.
        self._write_results(tmp_path, [(0, -1.0), (1, -2.0)])
        self._write_conformer(tmp_path, 0, "conf0")
        self._write_conformer(tmp_path, 1, "conf1")

        out = tmp_path / "ensemble.xyz"
        build_rdkit_ensemble_xyz(dir_path=tmp_path, output_xyz_file=out)

        text = out.read_text()
        assert text.index("conf1") < text.index("conf0")
        assert text.count("O 0.0 0.0 0.0") == 2

    def test_skips_missing_geometries(self, tmp_path):
        self._write_results(tmp_path, [(0, -1.0), (1, -2.0)])
        # Only conformer 0 has an xtbopt geometry on disk.
        self._write_conformer(tmp_path, 0, "conf0")

        out = tmp_path / "ensemble.xyz"
        build_rdkit_ensemble_xyz(dir_path=tmp_path, output_xyz_file=out)

        text = out.read_text()
        assert "conf0" in text
        assert "conf1" not in text

    def test_missing_results_csv_raises(self, tmp_path):
        out = tmp_path / "ensemble.xyz"
        with pytest.raises(FileNotFoundError, match="RDKit results"):
            build_rdkit_ensemble_xyz(dir_path=tmp_path, output_xyz_file=out)

    def test_no_geometries_raises(self, tmp_path):
        self._write_results(tmp_path, [(0, -1.0)])
        out = tmp_path / "ensemble.xyz"
        with pytest.raises(RuntimeError, match="no xTB-optimised"):
            build_rdkit_ensemble_xyz(dir_path=tmp_path, output_xyz_file=out)
