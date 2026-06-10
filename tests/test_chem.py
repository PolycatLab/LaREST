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


# ---------------------------------------------------------------------------
# Methyl formate initiator (non-symmetric, formate ester)
# ---------------------------------------------------------------------------

_METHYL_FORMATE = "COC=O"
# 5-membered ring with methyl substituent adjacent to ring O
_METHYL_GAMMA_BL = "CC1CCC(=O)O1"


def _has_substructure(smiles: str, smarts: str) -> bool:
    from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles

    mol = MolFromSmiles(smiles)
    pattern = MolFromSmarts(smarts)
    assert mol is not None, f"Invalid product SMILES: {smiles}"
    assert pattern is not None, f"Invalid SMARTS: {smarts}"
    return mol.HasSubstructMatch(pattern)


class TestBuildPolymerRORMethylFormate:
    """Methyl formate (HCOOCH3) breaks its acyl-C-O bond so both chain ends
    are esters: HC(=O)-O- at the ring-O end and -C(=O)-OCH3 at the acyl end.
    """

    def test_beta_pl_n1_valid_smiles(self):
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config)
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    def test_beta_pl_n1_has_formate_ester_end(self):
        """Product must have HC(=O)-O- (formate ester) at the ring-O end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config)
        # formyl C: exactly 1 H, double bond to O, single bond to O
        assert _has_substructure(smiles, "[C;H1](=O)[O]")

    def test_beta_pl_n1_has_methyl_ester_end(self):
        """Product must have -C(=O)-OCH3 (methyl ester) at the acyl end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config)
        # methyl ester terminus: ester carbonyl bonded to O bonded to CH3
        assert _has_substructure(smiles, "C(=O)O[CH3]")

    def test_beta_pl_n2_both_ester_ends(self):
        """Same terminal groups survive for a length-2 chain."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 2, "ROR", config)
        assert _has_substructure(smiles, "[C;H1](=O)[O]")
        assert _has_substructure(smiles, "C(=O)O[CH3]")

    def test_substituted_lactone_n1_has_formate_ester_end(self):
        """Substituted lactone (methyl at ring position) still gets formate end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_METHYL_GAMMA_BL, 1, "ROR", config)
        assert _has_substructure(smiles, "[C;H1](=O)[O]")

    def test_substituted_lactone_n1_has_methyl_ester_end(self):
        """Substituted lactone (methyl at ring position) still gets methyl ester end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_METHYL_GAMMA_BL, 1, "ROR", config)
        assert _has_substructure(smiles, "C(=O)O[CH3]")

    def test_substituted_lactone_retains_methyl_substituent(self):
        """The ring methyl group of the monomer must be present in the product."""
        from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles

        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_METHYL_GAMMA_BL, 1, "ROR", config)
        mol = MolFromSmiles(smiles)
        # product has at least two CH3 groups: one from initiator -OCH3,
        # one substituent from the monomer ring
        assert len(mol.GetSubstructMatches(MolFromSmarts("[CH3]"))) >= 2
