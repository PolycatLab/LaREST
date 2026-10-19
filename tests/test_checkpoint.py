"""Tests for larest.checkpoint result restoration and entropy correction."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from larest.checkpoint import (
    PipelineStage,
    apply_entropy_correction,
    combine_entropy,
    get_entropy_reference_section,
    restore_results,
    write_stage_settings,
)
from larest.constants import HARTTREE_TO_JMOL, PIPELINE_SECTIONS, THERMODYNAMIC_PARAMS

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
        f"{section}{suffix}": {"E": -95000.0, "H": -100000.0, "S": 50.0, "G": -115000.0}
        for section in [
            "censo_prescreening",
            "censo_screening",
            "censo_optimization",
            "censo_refinement",
        ]
        for suffix in ("", "_min")
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

    def test_disabled_crest_does_not_block_censo_checkpoint(self, tmp_path):
        """With crest_confgen off, a completed CENSO stage must still resume.

        The missing CREST checkpoint should be skipped rather than treated as
        the next stage to run, so the pipeline resumes past the (finished)
        expensive CENSO stage.
        """
        _make_rdkit_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        steps = {"crest_confgen": False}
        _, stage = restore_results(tmp_path, steps=steps)
        assert stage == PipelineStage.CREST_ENTROPY

    def test_disabled_crest_reruns_missing_censo(self, tmp_path):
        """With crest_confgen off and CENSO not yet done, resume at CENSO."""
        _make_rdkit_checkpoint(tmp_path)
        steps = {"crest_confgen": False}
        _, stage = restore_results(tmp_path, steps=steps)
        assert stage == PipelineStage.CENSO

    def test_disabled_crest_entropy_returns_finish(self, tmp_path):
        """A disabled trailing stage should not prevent reaching FINISH."""
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        steps = {"crest_entropy": False}
        _, stage = restore_results(tmp_path, steps=steps)
        assert stage == PipelineStage.FINISH


# ---------------------------------------------------------------------------
# Settings stamps and legacy results
# ---------------------------------------------------------------------------


class TestStageSettings:
    def _config(self):
        return {
            "steps": {"xtb": True},
            "rdkit": {"n_conformers": 10, "n_cores": 4},
            "xtb": {"gfn": 2, "alpb": "toluene", "parallel": 4},
            "crest": {"confgen": {"ewin": 6.0, "T": 4}, "entropy": {"T": 4}},
            "censo": {
                "cli": {"maxcores": 4, "refinement": True},
                "refinement": {"func": "wb97x-v"},
            },
        }

    def test_unchanged_settings_restore(self, tmp_path):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        write_stage_settings(tmp_path, PipelineStage.RDKIT, config, "C1CC(=O)O1")
        _, stage = restore_results(
            tmp_path,
            config=config,
            smiles="C1CC(=O)O1",
        )
        assert stage == PipelineStage.CREST_CONFGEN

    def test_changed_settings_rerun(self, tmp_path, caplog):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        write_stage_settings(tmp_path, PipelineStage.RDKIT, config, "C1CC(=O)O1")
        config["xtb"]["alpb"] = False
        _, stage = restore_results(tmp_path, config=config, smiles="C1CC(=O)O1")
        assert stage == PipelineStage.RDKIT
        assert "have changed" in caplog.text

    def test_core_counts_do_not_trigger_rerun(self, tmp_path):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        write_stage_settings(tmp_path, PipelineStage.RDKIT, config, "C1CC(=O)O1")
        config["xtb"]["parallel"] = 64
        config["rdkit"]["n_cores"] = 64
        _, stage = restore_results(tmp_path, config=config, smiles="C1CC(=O)O1")
        assert stage == PipelineStage.CREST_CONFGEN

    def test_different_molecule_reruns(self, tmp_path):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        write_stage_settings(tmp_path, PipelineStage.RDKIT, config, "C1CC(=O)O1")
        _, stage = restore_results(tmp_path, config=config, smiles="C1CCC(=O)O1")
        assert stage == PipelineStage.RDKIT

    def test_censo_settings_change_reruns_censo(self, tmp_path):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        for stage in (
            PipelineStage.RDKIT,
            PipelineStage.CREST_CONFGEN,
            PipelineStage.CENSO,
        ):
            write_stage_settings(tmp_path, stage, config, "C")
        config["censo"]["refinement"]["func"] = "b3lyp-d3"
        _, stage = restore_results(tmp_path, config=config, smiles="C")
        assert stage == PipelineStage.CENSO

    def test_legacy_stage_without_stamp_is_accepted_and_stamped(self, tmp_path):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        (tmp_path / "rdkit").mkdir()
        _, stage = restore_results(tmp_path, config=config, smiles="C")
        assert stage == PipelineStage.CREST_CONFGEN
        assert (tmp_path / "rdkit" / "larest_settings.json").exists()


class TestLegacyResults:
    def test_rdkit_energy_recovered_from_xtb_output(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        conf_dir = tmp_path / "xtb" / "rdkit" / "conformer_0"
        conf_dir.mkdir()
        (conf_dir / "conformer_0.txt").write_text(
            "| TOTAL ENERGY  -0.5 Eh |\n",
        )
        results, _ = restore_results(tmp_path)
        assert results["rdkit"]["E"] == pytest.approx(-0.5 * HARTTREE_TO_JMOL)

    def test_crest_energy_recovered_from_xtb_output(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        (tmp_path / "xtb" / "crest" / "xtb.txt").write_text(
            "| TOTAL ENERGY  -0.25 Eh |\n",
        )
        results, _ = restore_results(tmp_path)
        assert results["crest"]["E"] == pytest.approx(-0.25 * HARTTREE_TO_JMOL)

    def test_censo_rederived_from_raw_output(self, tmp_path):
        """Stored CENSO results.json (old format: H = E) is ignored if raw output exists."""
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        censo_dir = tmp_path / "censo"
        (censo_dir / "3_REFINEMENT.json").write_text(
            json.dumps(
                {
                    "data": {
                        "CONF1": {
                            "energy": -1.0,
                            "gsolv": 0.0,
                            "grrho": 0.04,
                            "gtot": -0.96,
                        },
                    },
                },
            ),
        )
        rrho = censo_dir / "3_REFINEMENT" / "CONF1" / "xtb_rrho" / "xtb_rrho.out"
        rrho.parent.mkdir(parents=True)
        rrho.write_text("| TOTAL ENERGY -5.0 Eh |\n| TOTAL ENTHALPY -4.93 Eh |\n")
        results, _ = restore_results(tmp_path)
        assert results["censo_refinement_min"]["H"] == pytest.approx(
            -0.93 * HARTTREE_TO_JMOL,
        )
        entropy = results["censo_refinement_min"]["S"]
        assert entropy is not None
        assert entropy > 0
        assert results["censo_screening"]["G"] is None

    def test_entropy_correction_uses_lowest_conformer(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        _make_crest_entropy_checkpoint(tmp_path)
        results, stage = restore_results(tmp_path)
        assert stage == PipelineStage.FINISH
        assert results["censo_corrected"]["S"] == pytest.approx(50.0 + 15.0)

    def test_entropy_without_censo_stored_raw(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_crest_entropy_checkpoint(tmp_path)
        results, stage = restore_results(tmp_path, steps={"censo": False})
        assert stage == PipelineStage.FINISH
        assert results["crest_entropy"]["S_total"] == pytest.approx(15.0)


class TestEntropyReferenceSection:
    def test_refinement_preferred(self):
        results = {
            "censo_optimization_min": {"H": -1.0, "S": 1.0},
            "censo_refinement_min": {"H": -2.0, "S": 2.0},
        }
        assert get_entropy_reference_section(results) == "censo_refinement_min"

    def test_falls_back_to_last_completed_stage(self):
        results = {
            "censo_optimization_min": {"H": -1.0, "S": 1.0},
            "censo_refinement_min": {"H": None, "S": None},
        }
        assert get_entropy_reference_section(results) == "censo_optimization_min"


class TestChainedSettings:
    def _config(self):
        return {
            "steps": {
                "rdkit": True,
                "crest_confgen": True,
                "censo": True,
                "crest_entropy": True,
                "xtb": True,
            },
            "rdkit": {"n_conformers": 10},
            "xtb": {"gfn": 2},
            "crest": {"confgen": {"ewin": 6.0}, "entropy": {}},
            "censo": {"refinement": {"func": "wb97x-v"}},
        }

    def test_upstream_change_invalidates_later_stage_after_interrupted_rerun(
        self,
        tmp_path,
    ):
        """An RDKit settings change must invalidate CREST even if the RDKit
        rerun completed (new stamp) but the CREST rerun was interrupted."""
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        for stage in (PipelineStage.RDKIT, PipelineStage.CREST_CONFGEN):
            write_stage_settings(tmp_path, stage, config, "C")
        config["rdkit"]["n_conformers"] = 50
        # RDKit was re-run with the new settings, then the job died
        write_stage_settings(tmp_path, PipelineStage.RDKIT, config, "C")
        _, stage = restore_results(
            tmp_path,
            steps={"censo": False, "crest_entropy": False},
            config=config,
            smiles="C",
        )
        assert stage == PipelineStage.CREST_CONFGEN

    def test_disabling_upstream_stage_invalidates_censo(self, tmp_path):
        config = self._config()
        _make_rdkit_checkpoint(tmp_path)
        _make_crest_checkpoint(tmp_path)
        _make_censo_checkpoint(tmp_path)
        for stage in (
            PipelineStage.RDKIT,
            PipelineStage.CREST_CONFGEN,
            PipelineStage.CENSO,
        ):
            write_stage_settings(tmp_path, stage, config, "C")
        config["steps"]["crest_confgen"] = False
        _, stage = restore_results(
            tmp_path,
            steps=config["steps"],
            config=config,
            smiles="C",
        )
        assert stage == PipelineStage.CENSO


class TestCrestCheckpointWithoutXtb:
    def test_crest_ensemble_is_checkpoint(self, tmp_path):
        _make_rdkit_checkpoint(tmp_path)
        crest_dir = tmp_path / "crest_confgen"
        crest_dir.mkdir()
        (crest_dir / "crest_conformers.xyz").write_text("1\n\nO 0 0 0\n")
        _, stage = restore_results(
            tmp_path,
            steps={"xtb": False, "censo": False, "crest_entropy": False},
        )
        assert stage == PipelineStage.FINISH


class TestCombineEntropy:
    def test_without_censo_enthalpy_stores_raw(self, caplog):
        results = {"censo_refinement_min": {"H": None, "S": None, "G": None}}
        entropy = {"S_conf": 1.0, "S_rrho": 0.5, "S_total": 1.5}
        combined = combine_entropy(results, entropy, 298.15)
        assert combined == {"crest_entropy": entropy}
        assert "No CENSO enthalpy" in caplog.text

    def test_with_censo_applies_correction(self):
        results = {"censo_refinement_min": {"E": -1.0, "H": -2.0, "S": 3.0, "G": -4.0}}
        entropy = {"S_conf": 1.0, "S_rrho": 0.5, "S_total": 1.5}
        combined = combine_entropy(results, entropy, 298.15)
        assert combined["censo_corrected"]["S"] == pytest.approx(4.5)

    def test_censo_disabled_stores_raw(self):
        results = {"censo_refinement_min": {"E": -1.0, "H": -2.0, "S": 3.0, "G": -4.0}}
        entropy = {"S_conf": 1.0, "S_rrho": 0.5, "S_total": 1.5}
        assert combine_entropy(results, entropy, 298.15, censo_enabled=False) == {
            "crest_entropy": entropy,
        }
