"""Tests for larest.censo output parsing utilities (no ORCA/CENSO required)."""

from __future__ import annotations

import json
import math

import pytest

from larest.censo import (
    build_rdkit_ensemble_xyz,
    check_censo_config,
    create_censorc,
    extract_best_conformer_xyz,
    get_best_censo_conformer,
    get_last_censo_section,
    parse_censo_part,
    parse_censo_results,
)
from larest.constants import CENSO_SECTIONS, GAS_CONSTANT, HARTTREE_TO_JMOL

TEMPERATURE = 298.15


def _write_xtb_rrho(path, energy, enthalpy):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "          | TOTAL ENERGY            "
        f"{energy:.12f} Eh   |\n"
        "          | TOTAL ENTHALPY          "
        f"{enthalpy:.12f} Eh   |\n"
        "          | TOTAL FREE ENERGY       -16.0 Eh   |\n",
    )


def _write_censo_part(censo_dir, part, conformers, *, rrho=True):
    """Write a CENSO <part>.json and matching xtb_rrho outputs.

    conformers maps conformer ID to (energy, gsolv, grrho, h_rrho) in Hartree.
    """
    censo_dir.mkdir(parents=True, exist_ok=True)
    data = {}
    for conf, (energy, gsolv, grrho, h_rrho) in conformers.items():
        data[conf] = {
            "energy": energy,
            "gsolv": gsolv,
            "grrho": grrho,
            "gtot": energy + gsolv + grrho,
        }
        if rrho:
            _write_xtb_rrho(
                censo_dir / part / conf / "xtb_rrho" / "xtb_rrho.out",
                -16.5,
                -16.5 + h_rrho,
            )
    (censo_dir / f"{part}.json").write_text(json.dumps({"data": data}))


_REFINEMENT = {
    # energy, gsolv, grrho, h_rrho (Hartree)
    "CONF1": (-267.1776, -0.0098, 0.0394, 0.0717),
    "CONF2": (-267.1760, -0.0100, 0.0390, 0.0712),
}


class TestParseCensoPart:
    def test_energies(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        conformers = parse_censo_part(tmp_path, "censo_refinement")
        energy, gsolv, grrho, h_rrho = _REFINEMENT["CONF1"]
        conf1 = next(c for c in conformers if c["conformer"] == "CONF1")
        assert conf1["E"] == pytest.approx((energy + gsolv) * HARTTREE_TO_JMOL)
        assert conf1["H"] == pytest.approx(
            (energy + gsolv + h_rrho) * HARTTREE_TO_JMOL,
        )
        assert conf1["G"] == pytest.approx(
            (energy + gsolv + grrho) * HARTTREE_TO_JMOL,
        )

    def test_enthalpy_above_free_energy(self, tmp_path):
        """H > G, i.e. a positive entropy (the old parser gave H < G)."""
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        for conf in parse_censo_part(tmp_path, "censo_refinement"):
            assert conf["H"] > conf["G"]

    def test_sorted_by_free_energy(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        conformers = parse_censo_part(tmp_path, "censo_refinement")
        assert [c["conformer"] for c in conformers] == ["CONF1", "CONF2"]

    def test_no_rrho_gives_no_h_or_g(self, tmp_path):
        _write_censo_part(
            tmp_path,
            "0_PRESCREENING",
            {"CONF1": (-266.6, -0.01, 0.0, 0.0)},
            rrho=False,
        )
        (conf,) = parse_censo_part(tmp_path, "censo_prescreening")
        assert conf["E"] == pytest.approx(-266.61 * HARTTREE_TO_JMOL)
        assert conf["H"] is None
        assert conf["G"] is None

    def test_missing_rrho_output_gives_no_h(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT, rrho=False)
        conformers = parse_censo_part(tmp_path, "censo_refinement")
        assert all(c["H"] is None for c in conformers)
        assert all(c["G"] is not None for c in conformers)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            parse_censo_part(tmp_path, "censo_refinement")


class TestParseCensoResults:
    def test_sections(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        result = parse_censo_results(tmp_path, TEMPERATURE)
        assert set(result) == {
            *CENSO_SECTIONS,
            *(f"{section}_min" for section in CENSO_SECTIONS),
        }
        assert result["censo_screening"]["G"] is None

    def test_minimum_is_lowest_free_energy_conformer(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        result = parse_censo_results(tmp_path, TEMPERATURE)
        conf1 = parse_censo_part(tmp_path, "censo_refinement")[0]
        for param in ("E", "H", "G"):
            assert result["censo_refinement_min"][param] == pytest.approx(conf1[param])

    def test_boltzmann_average(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        result = parse_censo_results(tmp_path, TEMPERATURE)
        conformers = parse_censo_part(tmp_path, "censo_refinement")
        rt = GAS_CONSTANT * TEMPERATURE
        weights = [math.exp(-(c["G"] - conformers[0]["G"]) / rt) for c in conformers]
        for param in ("E", "H", "G"):
            expected = sum(
                w * c[param] for w, c in zip(weights, conformers, strict=True)
            ) / sum(weights)
            assert result["censo_refinement"][param] == pytest.approx(expected)

    def test_entropy_is_positive_and_consistent(self, tmp_path):
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        result = parse_censo_results(tmp_path, TEMPERATURE)
        for section in ("censo_refinement", "censo_refinement_min"):
            h, s, g = (result[section][p] for p in ("H", "S", "G"))
            assert h is not None
            assert s is not None
            assert g is not None
            assert s > 0
            assert s == pytest.approx((h - g) / TEMPERATURE)

    def test_empty_dir_returns_none(self, tmp_path):
        result = parse_censo_results(tmp_path, TEMPERATURE)
        for params in result.values():
            assert all(value is None for value in params.values())


class TestBestCensoConformer:
    def test_best_conformer(self, tmp_path):
        _write_censo_part(
            tmp_path,
            "3_REFINEMENT",
            {"CONF1": (-1.0, 0.0, 0.01, 0.02), "CONF5": (-1.1, 0.0, 0.01, 0.02)},
        )
        assert get_best_censo_conformer(tmp_path, "censo_refinement") == "CONF5"

    def test_last_section(self, tmp_path):
        assert get_last_censo_section(tmp_path) is None
        _write_censo_part(tmp_path, "1_SCREENING", _REFINEMENT)
        assert get_last_censo_section(tmp_path) == "censo_screening"
        _write_censo_part(tmp_path, "3_REFINEMENT", _REFINEMENT)
        assert get_last_censo_section(tmp_path) == "censo_refinement"


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


class TestRunCensoMissingOutput:
    def test_raises_when_censo_writes_no_results(self, tmp_path, minimal_config):
        from unittest.mock import patch

        from larest.censo import run_censo

        crest_dir = tmp_path / "mol" / "crest_confgen"
        crest_dir.mkdir(parents=True)
        (crest_dir / "crest_conformers.xyz").write_text("1\nCONF1\nO 0 0 0\n")
        (tmp_path / "temp").mkdir()
        with (
            patch("larest.censo.subprocess.run"),
            pytest.raises(RuntimeError, match="wrote no"),
        ):
            run_censo(tmp_path / "mol", tmp_path, minimal_config)


class TestCheckCensoConfig:
    def test_recognised_keys_pass(self, minimal_config):
        pytest.importorskip("censo")
        minimal_config["censo"]["general"]["gas_phase"] = True
        minimal_config["censo"]["refinement"] = {"func": "b3lyp-d3"}
        minimal_config["censo"]["prescreening"] = {}
        minimal_config["censo"]["screening"] = {}
        minimal_config["censo"]["optimization"] = {}
        assert check_censo_config(minimal_config) == []

    def test_unknown_keys_reported(self, minimal_config, caplog):
        pytest.importorskip("censo")
        minimal_config["censo"]["general"]["bhess"] = True
        minimal_config["censo"]["bogus"] = {"x": 1}
        unknown = check_censo_config(minimal_config)
        assert "[general] bhess" in unknown
        assert "[bogus]" in unknown
        assert "will be ignored" in caplog.text
