"""Tests for larest.main pipeline orchestration."""

from __future__ import annotations

import subprocess
import textwrap
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from larest.data import MolResults, Monomer, Polymer
from larest.main import compile_results, format_results_table, run_pipeline

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_MONOMER_SMILES = "C1CC(=O)O1"
_INITIATOR_SMILES = "CCO"

_SECTIONS = [
    "rdkit",
    "crest",
    "censo_prescreening",
    "censo_screening",
    "censo_optimization",
    "censo_refinement",
]


def _make_mol_results(
    smiles: str,
    h: float = -100000.0,
    s: float = -50.0,
    g: float = -115000.0,
):
    sections: dict[str, dict[str, float | None]] = {
        section: {"H": h, "S": s, "G": g} for section in _SECTIONS
    }
    return MolResults(smiles=smiles, sections=sections)


# ---------------------------------------------------------------------------
# compile_results
# ---------------------------------------------------------------------------


class TestFormatResultsTable:
    def _make_summary_dfs(self):
        return {
            "rdkit": pd.DataFrame(
                {
                    "polymer_length": [2.0, 3.0],
                    "delta_H": [-5000.0, -4500.0],
                    "delta_S": [-10.0, -9.5],
                    "delta_G": [-2000.0, -1800.0],
                },
            ),
        }

    def test_contains_monomer_smiles(self):
        smiles = "C1CC(=O)O1"
        table = format_results_table(smiles, self._make_summary_dfs())
        assert smiles in table

    def test_contains_section_name(self):
        table = format_results_table("C1CC(=O)O1", self._make_summary_dfs())
        assert "RDKit" in table

    def test_contains_delta_columns_with_units(self):
        table = format_results_table("C1CC(=O)O1", self._make_summary_dfs())
        assert "ΔH (kJ/mol)" in table
        assert "ΔS (J/mol/K)" in table
        assert "ΔG (kJ/mol)" in table

    def test_contains_polymer_lengths(self):
        table = format_results_table("C1CC(=O)O1", self._make_summary_dfs())
        assert "2" in table
        assert "3" in table

    def test_contains_delta_values(self):
        table = format_results_table("C1CC(=O)O1", self._make_summary_dfs())
        assert "-5.0000" in table
        assert "-4.5000" in table

    def test_multiple_sections(self):
        dfs = {
            "rdkit": pd.DataFrame(
                {
                    "polymer_length": [2.0],
                    "delta_H": [1.0],
                    "delta_S": [2.0],
                    "delta_G": [3.0],
                },
            ),
            "crest": pd.DataFrame(
                {
                    "polymer_length": [2.0],
                    "delta_H": [4.0],
                    "delta_S": [5.0],
                    "delta_G": [6.0],
                },
            ),
        }
        table = format_results_table("C1CC(=O)O1", dfs)
        assert "RDKit" in table
        assert "CREST" in table

    def test_returns_string(self):
        table = format_results_table("C1CC(=O)O1", self._make_summary_dfs())
        assert isinstance(table, str)

    def test_displays_inf_row_as_unicode(self):
        import numpy as np

        dfs = {
            "rdkit": pd.DataFrame(
                {
                    "polymer_length": [2.0, np.inf],
                    "delta_H": [-5000.0, -4000.0],
                    "delta_S": [-10.0, -8.0],
                    "delta_G": [-2000.0, -1600.0],
                },
            ),
        }
        table = format_results_table("C1CC(=O)O1", dfs)
        assert "∞" in table

    def test_skips_all_nan_sections(self):
        dfs = {
            "rdkit": pd.DataFrame(
                {
                    "polymer_length": [2.0],
                    "delta_H": [1.0],
                    "delta_S": [2.0],
                    "delta_G": [3.0],
                },
            ),
            "crest": pd.DataFrame(
                {
                    "polymer_length": [2.0],
                    "delta_H": [float("nan")],
                    "delta_S": [float("nan")],
                    "delta_G": [float("nan")],
                },
            ),
        }
        table = format_results_table("C1CC(=O)O1", dfs)
        assert "RDKit" in table
        assert "CREST" not in table


class TestCompileResults:
    def test_returns_dict_of_dataframes(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [
            (2, _make_mol_results("P", h=-210000.0, s=-105.0, g=-241500.0)),
        ]

        result = compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        import numpy as np

        assert isinstance(result, dict)
        assert set(result.keys()) == set(_SECTIONS)
        for df in result.values():
            assert isinstance(df, pd.DataFrame)
            assert "delta_H" in df.columns
            assert "delta_S" in df.columns
            assert "delta_G" in df.columns
            assert np.isinf(df["polymer_length"]).any()

    def test_creates_summary_directory(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [
            (
                2,
                _make_mol_results(
                    "OCC(=O)OCC(=O)O",
                    h=-200000.0,
                    s=-100.0,
                    g=-230000.0,
                ),
            ),
        ]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        assert summary_dir.exists()

    def test_creates_csv_per_section(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [
            (
                2,
                _make_mol_results(
                    "OCC(=O)OCC(=O)O",
                    h=-200000.0,
                    s=-100.0,
                    g=-230000.0,
                ),
            ),
        ]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        for section in _SECTIONS:
            assert (summary_dir / f"{section}.csv").exists()

    def test_delta_computation_rer(self, tmp_path):
        """delta_param = (poly - n * mono - 0) / n"""
        monomer_results = _make_mol_results(
            _MONOMER_SMILES,
            h=-100000.0,
            s=-50.0,
            g=-115000.0,
        )
        n = 2
        polymer_results = [
            (n, _make_mol_results("P", h=-210000.0, s=-105.0, g=-241500.0)),
        ]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        import pandas as pd

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        row = df[df["polymer_length"] == n].iloc[0]

        expected_delta_h = (-210000.0 - n * -100000.0 - 0) / n
        assert row["delta_H"] == pytest.approx(expected_delta_h, rel=1e-6)

    def test_ror_raises_without_initiator(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [(1, _make_mol_results("P"))]

        with pytest.raises(ValueError, match="Initiator results are required"):
            compile_results(
                monomer_smiles=_MONOMER_SMILES,
                monomer_results=monomer_results,
                polymer_results=polymer_results,
                initiator_results=None,
                output_dir=tmp_path,
                reaction_type="ROR",
            )

    def test_ror_includes_initiator_in_delta(self, tmp_path):
        monomer_results = _make_mol_results(
            _MONOMER_SMILES,
            h=-100000.0,
            s=-50.0,
            g=-115000.0,
        )
        initiator_results = _make_mol_results(
            _INITIATOR_SMILES,
            h=-50000.0,
            s=-30.0,
            g=-59000.0,
        )
        n = 1
        polymer_results = [
            (n, _make_mol_results("P", h=-160000.0, s=-85.0, g=-180000.0)),
        ]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=initiator_results,
            output_dir=tmp_path,
            reaction_type="ROR",
        )

        import pandas as pd

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        row = df[df["polymer_length"] == n].iloc[0]

        expected_delta_h = (-160000.0 - n * -100000.0 - -50000.0) / n
        assert row["delta_H"] == pytest.approx(expected_delta_h, rel=1e-6)

    def test_inf_row_extrapolated_value(self, tmp_path):
        """Intercept of ΔH vs 1/n regression is written as the inf row."""
        import numpy as np

        # delta_H(n=2) = (-210000 - 2*-100000) / 2 = -5000
        # delta_H(n=4) = (-418000 - 4*-100000) / 4 = -4500
        # polyfit(1/n=[0.5, 0.25], y=[-5000, -4500], 1) → intercept = -4000
        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=_make_mol_results(
                _MONOMER_SMILES,
                h=-100000.0,
                s=-50.0,
                g=-115000.0,
            ),
            polymer_results=[
                (2, _make_mol_results("P2", h=-210000.0, s=-105.0, g=-241500.0)),
                (4, _make_mol_results("P4", h=-418000.0, s=-209.0, g=-480200.0)),
            ],
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        inf_row = df[np.isinf(df["polymer_length"])].iloc[0]

        # delta_H(n=2) = (-210000 - 2*-100000) / 2 = -5000
        # delta_H(n=4) = (-418000 - 4*-100000) / 4 = -4500
        # polyfit(1/n=[0.5,0.25], y=[-5000,-4500], 1) → intercept = -4000
        assert inf_row["delta_H"] == pytest.approx(-4000.0, rel=1e-6)

    def test_multiple_polymer_lengths(self, tmp_path):
        monomer_results = _make_mol_results(_MONOMER_SMILES)
        polymer_results = [
            (2, _make_mol_results("P2")),
            (3, _make_mol_results("P3")),
            (4, _make_mol_results("P4")),
        ]

        compile_results(
            monomer_smiles=_MONOMER_SMILES,
            monomer_results=monomer_results,
            polymer_results=polymer_results,
            initiator_results=None,
            output_dir=tmp_path,
            reaction_type="RER",
        )

        import numpy as np

        from larest.output import slugify

        summary_dir = tmp_path / "Monomer" / slugify(_MONOMER_SMILES) / "summary"
        df = pd.read_csv(summary_dir / "rdkit.csv")
        assert len(df) == 4  # 3 finite lengths + 1 inf row
        finite = df[np.isfinite(df["polymer_length"])]
        assert set(finite["polymer_length"].astype(int).tolist()) == {2, 3, 4}
        assert df["polymer_length"].apply(np.isinf).any()


# ---------------------------------------------------------------------------
# run_pipeline (mocked)
# ---------------------------------------------------------------------------


class TestRunPipelineMocked:
    def _make_config(self):
        return {
            "steps": {
                "rdkit": True,
                "crest_confgen": False,
                "censo": False,
                "crest_entropy": False,
                "xtb": False,
            },
            "rdkit": {
                "n_conformers": 1,
                "random_seed": 42,
                "conformer_box_size": 2.0,
                "mmff": "MMFF94",
                "mmff_iters": 10,
                "align_iters": 5,
                "precision": 6,
                "n_cores": 1,
            },
            "xtb": {
                "parallel": 1,
                "gfn": 2,
                "ohess": "vtight",
                "alpb": "toluene",
                "etemp": 298.15,
            },
            "crest": {"confgen": {}, "entropy": {}},
            "censo": {"cli": {}, "general": {"temperature": 298.15}},
            "reaction": {
                "type": "RER",
                "lengths": [2],
                "initiator": "",
                "monomers": [_MONOMER_SMILES],
            },
        }

    def _mol_dir(self, tmp_path, mol):
        """Return and pre-create the molecule output directory.

        run_rdkit normally creates this via create_dir(rdkit_dir, parents=True).
        When run_rdkit is mocked we must create it ourselves so results.json
        can be written.
        """
        from larest.output import create_dir, slugify

        match mol:
            case Polymer():
                d = tmp_path / "Polymer" / f"{slugify(mol.monomer_smiles)}_{mol.length}"
            case _:
                d = tmp_path / type(mol).__name__ / slugify(mol.smiles)
        create_dir(d)
        return d

    def test_run_pipeline_monomer(self, tmp_path):
        config = self._make_config()
        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        mol = Monomer(_MONOMER_SMILES)
        self._mol_dir(tmp_path, mol)

        with patch("larest.main.run_rdkit", return_value=fake_rdkit):
            result = run_pipeline(mol, tmp_path, config)

        assert isinstance(result, MolResults)
        assert result.smiles == _MONOMER_SMILES
        assert result.sections["rdkit"]["H"] == pytest.approx(-100000.0)

    def test_run_pipeline_writes_results_json(self, tmp_path):
        config = self._make_config()
        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        mol = Monomer(_MONOMER_SMILES)
        mol_dir = self._mol_dir(tmp_path, mol)

        with patch("larest.main.run_rdkit", return_value=fake_rdkit):
            run_pipeline(mol, tmp_path, config)

        assert (mol_dir / "results.json").exists()

    def test_run_pipeline_polymer_dir_name(self, tmp_path):
        config = self._make_config()
        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        polymer = Polymer(
            smiles="OCC(=O)OCC(=O)O",
            monomer_smiles=_MONOMER_SMILES,
            length=2,
        )
        expected_dir = self._mol_dir(tmp_path, polymer)

        with patch("larest.main.run_rdkit", return_value=fake_rdkit):
            run_pipeline(polymer, tmp_path, config)

        assert expected_dir.exists()

    def test_run_pipeline_skips_disabled_steps(self, tmp_path):
        config = self._make_config()
        config["steps"]["rdkit"] = False
        mol = Monomer(_MONOMER_SMILES)
        self._mol_dir(tmp_path, mol)

        with (
            patch("larest.main.run_rdkit") as mock_rdkit,
            patch("larest.main.restore_results", return_value=({}, MagicMock(value=5))),
        ):
            run_pipeline(mol, tmp_path, config)

        mock_rdkit.assert_not_called()

    def test_run_pipeline_stores_crest_entropy_when_censo_disabled(self, tmp_path):
        from larest.checkpoint import PipelineStage

        config = self._make_config()
        config["steps"]["crest_entropy"] = True
        config["steps"]["censo"] = False
        mol = Monomer(_MONOMER_SMILES)
        self._mol_dir(tmp_path, mol)

        fake_rdkit = {"rdkit": {"H": -100000.0, "S": -50.0, "G": -115000.0}}
        fake_entropy = {"S_conf": 10.0, "S_rrho": 5.0, "S_total": 15.0}

        with (
            patch(
                "larest.main.restore_results",
                return_value=({}, PipelineStage.RDKIT),
            ),
            patch("larest.main.run_rdkit", return_value=fake_rdkit),
            patch("larest.main.run_crest_entropy", return_value=fake_entropy),
        ):
            result = run_pipeline(mol, tmp_path, config)

        assert "crest_entropy" in result.sections
        assert result.sections["crest_entropy"] == fake_entropy


# ---------------------------------------------------------------------------
# Full CLI integration test
# ---------------------------------------------------------------------------

_MINIMAL_CONFIG_TOML = textwrap.dedent("""\
    [reaction]
    monomers = ["C1CC(=O)O1"]
    lengths = [2]
    type = "RER"

    [steps]
    rdkit = true
    crest_confgen = false
    censo = false
    crest_entropy = false
    xtb = true

    [rdkit]
    n_conformers = 2
""")


@pytest.mark.integration
class TestCLIIntegration:
    """End-to-end test: invokes the `larest` CLI on a minimal config."""

    def test_cli_exits_zero(self, tmp_path):
        config_file = tmp_path / "config.toml"
        config_file.write_text(_MINIMAL_CONFIG_TOML)
        output_dir = tmp_path / "output"

        result = subprocess.run(
            ["larest", str(config_file), "-o", str(output_dir)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr

    def test_cli_creates_monomer_results(self, tmp_path):
        config_file = tmp_path / "config.toml"
        config_file.write_text(_MINIMAL_CONFIG_TOML)
        output_dir = tmp_path / "output"

        subprocess.run(
            ["larest", str(config_file), "-o", str(output_dir)],
            capture_output=True,
            check=True,
        )

        monomer_results = list((output_dir / "Monomer").glob("*/results.json"))
        assert len(monomer_results) == 1

    def test_cli_creates_polymer_results(self, tmp_path):
        config_file = tmp_path / "config.toml"
        config_file.write_text(_MINIMAL_CONFIG_TOML)
        output_dir = tmp_path / "output"

        subprocess.run(
            ["larest", str(config_file), "-o", str(output_dir)],
            capture_output=True,
            check=True,
        )

        polymer_results = list((output_dir / "Polymer").glob("*_2/results.json"))
        assert len(polymer_results) == 1

    def test_cli_creates_summary_csv(self, tmp_path):
        config_file = tmp_path / "config.toml"
        config_file.write_text(_MINIMAL_CONFIG_TOML)
        output_dir = tmp_path / "output"

        subprocess.run(
            ["larest", str(config_file), "-o", str(output_dir)],
            capture_output=True,
            check=True,
        )

        summary_csvs = list((output_dir / "Monomer").glob("*/summary/rdkit.csv"))
        assert len(summary_csvs) == 1

    def test_cli_summary_has_correct_polymer_length(self, tmp_path):
        import pandas as pd

        config_file = tmp_path / "config.toml"
        config_file.write_text(_MINIMAL_CONFIG_TOML)
        output_dir = tmp_path / "output"

        subprocess.run(
            ["larest", str(config_file), "-o", str(output_dir)],
            capture_output=True,
            check=True,
        )

        import numpy as np

        csv_path = next((output_dir / "Monomer").glob("*/summary/rdkit.csv"))
        df = pd.read_csv(csv_path)
        finite = df[np.isfinite(df["polymer_length"])]
        assert set(finite["polymer_length"].astype(int).tolist()) == {2}
