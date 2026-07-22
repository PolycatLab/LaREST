"""Tests for larest.checkpoint result restoration and entropy correction."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from larest.checkpoint import PipelineStage, apply_entropy_correction, restore_results
from larest.constants import PIPELINE_SECTIONS, THERMODYNAMIC_PARAMS

# ---------------------------------------------------------------------------
# PipelineStage
# ---------------------------------------------------------------------------


class TestPipelineStage:
    def test_ordering(self):
        assert PipelineStage.RDKIT < PipelineStage.CREST_CONFGEN
        assert PipelineStage.CREST_CONFGEN < PipelineStage.CENSO
        assert PipelineStage.CENSO < PipelineStage.CREST_ENTROPY
        assert PipelineStage.CREST_ENTROPY < PipelineStage.FINISH

    def test_values(self):
        assert PipelineStage.RDKIT == 1
        assert PipelineStage.FINISH == 5


# ---------------------------------------------------------------------------
# apply_entropy_correction
# ---------------------------------------------------------------------------


class TestApplyEntropyCorrection:
    def test_corrects_entropy(self):
        refinement: dict[str, float | None] = {
            "H": -100000.0,
            "S": -50.0,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": 15.0,
        }
        result = apply_entropy_correction(refinement, entropy, 298.15)
        assert result["censo_corrected"]["S"] == pytest.approx(-50.0 + 15.0)

    def test_preserves_h_and_recomputes_g(self):
        refinement: dict[str, float | None] = {
            "H": -100000.0,
            "S": -50.0,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": 15.0,
        }
        result = apply_entropy_correction(refinement, entropy, 298.15)
        # H is unchanged; G is recomputed from the corrected entropy as H - T*S
        corrected_s = -50.0 + 15.0
        assert result["censo_corrected"]["H"] == pytest.approx(-100000.0)
        assert result["censo_corrected"]["G"] == pytest.approx(
            -100000.0 - 298.15 * corrected_s,
        )

    def test_does_not_mutate_input(self):
        refinement: dict[str, float | None] = {
            "H": -100000.0,
            "S": -50.0,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": 15.0,
        }
        apply_entropy_correction(refinement, entropy, 298.15)
        assert refinement["S"] == pytest.approx(-50.0)

    def test_raises_when_s_none(self):
        refinement: dict[str, float | None] = {
            "H": -100000.0,
            "S": None,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": 15.0,
        }
        with pytest.raises(
            ValueError,
            match="Failed to apply CREST entropy correction",
        ):
            apply_entropy_correction(refinement, entropy, 298.15)

    def test_raises_when_s_total_none(self):
        refinement: dict[str, float | None] = {
            "H": -100000.0,
            "S": -50.0,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": None,
        }
        with pytest.raises(
            ValueError,
            match="Failed to apply CREST entropy correction",
        ):
            apply_entropy_correction(refinement, entropy, 298.15)

    def test_raises_when_h_none(self):
        refinement: dict[str, float | None] = {
            "H": None,
            "S": -50.0,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": 15.0,
        }
        with pytest.raises(
            ValueError,
            match="Failed to apply CREST entropy correction",
        ):
            apply_entropy_correction(refinement, entropy, 298.15)

    def test_returns_censo_corrected_key(self):
        refinement: dict[str, float | None] = {
            "H": -100000.0,
            "S": -50.0,
            "G": -115000.0,
        }
        entropy: dict[str, float | None] = {
            "S_conf": 10.0,
            "S_rrho": 5.0,
            "S_total": 15.0,
        }
        result = apply_entropy_correction(refinement, entropy, 298.15)
        assert "censo_corrected" in result


# ---------------------------------------------------------------------------
# restore_results
# ---------------------------------------------------------------------------


def _make_rdkit_checkpoint(dir_path: Path) -> None:
    """Create xtb/rdkit/results.csv with one row."""
    xtb_dir = dir_path / "xtb" / "rdkit"
    xtb_dir.mkdir(parents=True)
    (xtb_dir / "results.csv").write_text(
        "conformer_id,H,S,G\n0,-324567.89,-85.12,-350000.00\n",
    )


def _make_crest_checkpoint(dir_path: Path) -> None:
    crest_dir = dir_path / "xtb" / "crest"
    crest_dir.mkdir(parents=True)
    (crest_dir / "results.json").write_text(
        json.dumps({"H": -324567.0, "S": -85.0, "G": -350000.0}),
    )


def _make_censo_checkpoint(dir_path: Path) -> None:
    censo_dir = dir_path / "censo"
    censo_dir.mkdir(parents=True)
    sections_data = {
        section: {"H": -100000.0, "S": -50.0, "G": -115000.0}
        for section in [
            "censo_prescreening",
            "censo_screening",
            "censo_optimization",
            "censo_refinement",
        ]
    }
    (censo_dir / "results.json").write_text(json.dumps(sections_data))


def _make_crest_entropy_checkpoint(dir_path: Path) -> None:
    entropy_dir = dir_path / "crest_entropy"
    entropy_dir.mkdir(parents=True)
    (entropy_dir / "results.json").write_text(
        json.dumps({"S_conf": 10.0, "S_rrho": 5.0, "S_total": 15.0}),
    )


class TestRestoreResults:
    def test_empty_dir_returns_rdkit_stage(self, tmp_path):
        results, stage = restore_results(tmp_path)
        assert stage == PipelineStage.RDKIT

    def test_empty_results_have_none_values(self, tmp_path):
        results, _ = restore_results(tmp_path)
        for section in PIPELINE_SECTIONS:
            for param in THERMODYNAMIC_PARAMS:
                assert results[section][param] is None

    def test_after_rdkit_returns_crest_stage(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        results, stage = restore_results(tmp_path)
        assert stage == PipelineStage.CREST_CONFGEN

    def test_rdkit_results_loaded(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        results, _ = restore_results(tmp_path)
        assert results["rdkit"]["H"] == pytest.approx(-324567.89)
        assert results["rdkit"]["G"] == pytest.approx(-350000.00)

    def test_after_crest_returns_censo_stage(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        results, stage = restore_results(tmp_path)
        assert stage == PipelineStage.CENSO

    def test_crest_results_loaded(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        results, _ = restore_results(tmp_path)
        assert results["crest"]["H"] == pytest.approx(-324567.0)

    def test_after_censo_returns_crest_entropy_stage(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        _, stage = restore_results(tmp_path)
        assert stage == PipelineStage.CREST_ENTROPY

    def test_all_stages_present_returns_finish(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        _make_crest_entropy_checkpoint(tmp_path)
        _, stage = restore_results(tmp_path)
        assert stage == PipelineStage.FINISH

    def test_corrupt_checkpoint_falls_back(self, tmp_path):
        """A corrupt rdkit results.csv should restart from RDKIT."""
        xtb_dir = tmp_path / "xtb" / "rdkit"
        xtb_dir.mkdir(parents=True)
        (xtb_dir / "results.csv").write_text("not,a,valid,csv\n!!!\n")
        _, stage = restore_results(tmp_path)
        assert stage == PipelineStage.RDKIT
