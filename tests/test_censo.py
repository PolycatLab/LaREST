"""Tests for larest.censo output parsing utilities (no ORCA/CENSO required)."""

from __future__ import annotations

import pytest

from larest.censo import (
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

    def test_missing_data_defaults_to_conf0(self, tmp_path):
        empty = tmp_path / "empty.txt"
        empty.write_text("nothing here\n")
        result = parse_best_censo_conformers(empty)
        for section in CENSO_SECTIONS:
            assert result[section] == "CONF0"

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
