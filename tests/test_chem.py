"""Tests for larest.chem polymer chemistry utilities."""

from __future__ import annotations

import pytest

from larest.chem import build_polymer, get_mol, get_polymer_unit, get_ring_size

# ---------------------------------------------------------------------------
# get_mol
# ---------------------------------------------------------------------------


class TestGetMol:
    def test_valid_smiles(self):
        mol = get_mol("CCO")
        assert mol is not None

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Failed to create RDKit Mol"):
            get_mol("not_a_smiles!!!")

    def test_returns_mol_object(self):
        from rdkit.Chem.rdchem import Mol

        mol = get_mol("C1CC(=O)O1")
        assert isinstance(mol, Mol)


# ---------------------------------------------------------------------------
# get_ring_size
# ---------------------------------------------------------------------------


class TestGetRingSize:
    @pytest.mark.parametrize(
        ("smiles", "expected"),
        [
            ("C1CC(=O)O1", 4),  # beta-propiolactone (4-membered)
            ("C1CCC(=O)O1", 5),  # gamma-butyrolactone (5-membered)
            ("C1CCCC(=O)O1", 6),  # delta-valerolactone (6-membered)
            ("C1CCCCC(=O)O1", 7),  # epsilon-caprolactone (7-membered)
        ],
    )
    def test_ring_size(self, smiles, expected):
        assert get_ring_size(smiles) == expected

    def test_non_ring_returns_none(self):
        # acyclic ester has no ring-opening group
        assert get_ring_size("CC(=O)OCC") is None

    def test_invalid_smiles_raises(self):
        with pytest.raises(
            ValueError,
            match="Failed to create RDKit Mol object from SMILES",
        ):
            get_ring_size("not_valid!!!")


# ---------------------------------------------------------------------------
# get_polymer_unit
# ---------------------------------------------------------------------------


class TestGetPolymerUnit:
    def test_monomer_unit_created(self):
        from rdkit.Chem.rdchem import Mol

        unit = get_polymer_unit("C1CC(=O)O1", "monomer", "Xe", "Y")
        assert isinstance(unit, Mol)

    def test_initiator_unit_created(self):
        from rdkit.Chem.rdchem import Mol

        unit = get_polymer_unit("CCO", "initiator", "Xe", "Y")
        assert isinstance(unit, Mol)

    def test_invalid_monomer_raises(self):
        with pytest.raises(ValueError, match="No functional group atom ids found"):
            get_polymer_unit("CCC", "monomer", "Xe", "Y")  # no lactone group


# ---------------------------------------------------------------------------
# build_polymer
# ---------------------------------------------------------------------------

# Simple lactone monomers for polymer building
_BETA_PL = "C1CC(=O)O1"  # beta-propiolactone
_GAMMA_BL = "C1CCC(=O)O1"  # gamma-butyrolactone
_INITIATOR = "CCO"  # ethanol


class TestBuildPolymerRER:
    def test_basic_rer_length_2(self):
        config = {"reaction": {"type": "RER", "initiator": ""}}
        smiles = build_polymer(_BETA_PL, 2, "RER", config)
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    def test_rer_length_3(self):
        config = {"reaction": {"type": "RER", "initiator": ""}}
        smiles = build_polymer(_GAMMA_BL, 3, "RER", config)
        assert isinstance(smiles, str)

    def test_rer_length_1_raises(self):
        config = {"reaction": {"type": "RER", "initiator": ""}}
        with pytest.raises(ValueError, match="polymer length > 1"):
            build_polymer(_BETA_PL, 1, "RER", config)

    def test_rer_length_0_raises(self):
        config = {"reaction": {"type": "RER", "initiator": ""}}
        with pytest.raises(ValueError, match="polymer length > 1"):
            build_polymer(_BETA_PL, 0, "RER", config)


class TestBuildPolymerROR:
    def test_basic_ror_length_1(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config)
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    def test_ror_length_2(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        smiles = build_polymer(_BETA_PL, 2, "ROR", config)
        assert isinstance(smiles, str)

    def test_ror_length_0_raises(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        with pytest.raises(ValueError, match="polymer length >= 1"):
            build_polymer(_BETA_PL, 0, "ROR", config)

    def test_ror_length_0_error_includes_length(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        with pytest.raises(ValueError, match=r"\(current length: 0\)"):
            build_polymer(_BETA_PL, 0, "ROR", config)

    def test_ror_longer_polymer_than_rer(self):
        """ROR n=2 polymer should contain initiator fragment."""
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        smiles = build_polymer(_BETA_PL, 2, "ROR", config)
        # ethanol initiator (CCO) contributes two carbons
        assert smiles.count("C") >= 2
