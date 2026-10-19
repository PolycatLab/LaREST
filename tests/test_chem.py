"""Tests for larest.chem polymer chemistry utilities."""

from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles, MolToSmiles

from larest.chem import (
    analyse_monomer,
    build_polymer,
    get_isodesmic_initiator,
    get_mol,
    get_polymer_unit,
    get_ring_size,
)
from larest.constants import MONOMER_GROUPS

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
        smiles = build_polymer(_BETA_PL, 2, "RER")
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    def test_rer_length_3(self):
        smiles = build_polymer(_GAMMA_BL, 3, "RER")
        assert isinstance(smiles, str)

    def test_rer_length_1_raises(self):
        with pytest.raises(ValueError, match="polymer length > 1"):
            build_polymer(_BETA_PL, 1, "RER")

    def test_rer_length_0_raises(self):
        with pytest.raises(ValueError, match="polymer length > 1"):
            build_polymer(_BETA_PL, 0, "RER")


class TestBuildPolymerROR:
    def test_basic_ror_length_1(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config["reaction"]["initiator"])
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    def test_ror_length_2(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        smiles = build_polymer(_BETA_PL, 2, "ROR", config["reaction"]["initiator"])
        assert isinstance(smiles, str)

    def test_ror_length_0_raises(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        with pytest.raises(ValueError, match="polymer length >= 1"):
            build_polymer(_BETA_PL, 0, "ROR", config["reaction"]["initiator"])

    def test_ror_length_0_error_includes_length(self):
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        with pytest.raises(ValueError, match=r"\(current length: 0\)"):
            build_polymer(_BETA_PL, 0, "ROR", config["reaction"]["initiator"])

    def test_ror_longer_polymer_than_rer(self):
        """ROR n=2 polymer should contain initiator fragment."""
        config = {"reaction": {"type": "ROR", "initiator": _INITIATOR}}
        smiles = build_polymer(_BETA_PL, 2, "ROR", config["reaction"]["initiator"])
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
        smiles = build_polymer(_BETA_PL, 1, "ROR", config["reaction"]["initiator"])
        assert isinstance(smiles, str)
        assert len(smiles) > 0

    def test_beta_pl_n1_has_formate_ester_end(self):
        """Product must have HC(=O)-O- (formate ester) at the ring-O end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config["reaction"]["initiator"])
        # formyl C: exactly 1 H, double bond to O, single bond to O
        assert _has_substructure(smiles, "[C;H1](=O)[O]")

    def test_beta_pl_n1_has_methyl_ester_end(self):
        """Product must have -C(=O)-OCH3 (methyl ester) at the acyl end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 1, "ROR", config["reaction"]["initiator"])
        # methyl ester terminus: ester carbonyl bonded to O bonded to CH3
        assert _has_substructure(smiles, "C(=O)O[CH3]")

    def test_beta_pl_n2_both_ester_ends(self):
        """Same terminal groups survive for a length-2 chain."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(_BETA_PL, 2, "ROR", config["reaction"]["initiator"])
        assert _has_substructure(smiles, "[C;H1](=O)[O]")
        assert _has_substructure(smiles, "C(=O)O[CH3]")

    def test_substituted_lactone_n1_has_formate_ester_end(self):
        """Substituted lactone (methyl at ring position) still gets formate end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(
            _METHYL_GAMMA_BL,
            1,
            "ROR",
            config["reaction"]["initiator"],
        )
        assert _has_substructure(smiles, "[C;H1](=O)[O]")

    def test_substituted_lactone_n1_has_methyl_ester_end(self):
        """Substituted lactone (methyl at ring position) still gets methyl ester end."""
        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(
            _METHYL_GAMMA_BL,
            1,
            "ROR",
            config["reaction"]["initiator"],
        )
        assert _has_substructure(smiles, "C(=O)O[CH3]")

    def test_substituted_lactone_retains_methyl_substituent(self):
        """The ring methyl group of the monomer must be present in the product."""
        from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles

        config = {"reaction": {"type": "ROR", "initiator": _METHYL_FORMATE}}
        smiles = build_polymer(
            _METHYL_GAMMA_BL,
            1,
            "ROR",
            config["reaction"]["initiator"],
        )
        mol = MolFromSmiles(smiles)
        # product has at least two CH3 groups: one from initiator -OCH3,
        # one substituent from the monomer ring
        assert len(mol.GetSubstructMatches(MolFromSmarts("[CH3]"))) >= 2


# ---------------------------------------------------------------------------
# Bond selection and orientation
# ---------------------------------------------------------------------------


def _same_molecule(smiles_a: str, smiles_b: str) -> bool:
    return MolToSmiles(MolFromSmiles(smiles_a)) == MolToSmiles(MolFromSmiles(smiles_b))


class TestBuildPolymerBondSelection:
    @pytest.mark.parametrize(
        ("monomer", "initiator", "expected"),
        [
            # gamma-butyrolactone + methanol -> methyl 4-hydroxybutanoate
            ("O=C1CCCO1", "CO", "COC(=O)CCCO"),
            # alpha-methylene-gamma-butyrolactone: the acyl-O bond is broken,
            # not the bond between the carbonyl and the alkene carbon
            ("C=C1CCOC1=O", "CO", "C=C(CCO)C(=O)OC"),
            # carbonate: the initiator O attacks the carbonyl C (no peroxide)
            ("O=C1OCCCO1", "CO", "COC(=O)OCCCO"),
            # cyclic carbamate
            ("O=C1OCCN1", "CO", "COC(=O)NCCO"),
            # lactam
            ("O=C1CCCN1", "CO", "COC(=O)CCCN"),
            # water (hydrolysis)
            ("O=C1CCCO1", "O", "O=C(O)CCCO"),
            # methyl acetate (isodesmic, Houk/Hall)
            ("O=C1CCCO1", "CC(=O)OC", "COC(=O)CCCOC(C)=O"),
            # phthalide (aromatic alpha carbon)
            ("O=C1OCc2ccccc21", "CO", "COC(=O)c1ccccc1CO"),
        ],
    )
    def test_ror_products(self, monomer, initiator, expected):
        assert _same_molecule(build_polymer(monomer, 1, "ROR", initiator), expected)

    def test_ror_without_initiator_raises(self):
        with pytest.raises(ValueError, match="initiator"):
            build_polymer(_BETA_PL, 1, "ROR")

    def test_group_override(self):
        # glutarimide-like ring with a lactam and a cyclic carbamate bond
        smiles = "O=C1OCCC(=O)N1"
        default = build_polymer(smiles, 1, "ROR", "CO")
        override = build_polymer(smiles, 1, "ROR", "CO", group="oA")
        assert default != override
        assert _same_molecule(override, "COC(=O)NC(=O)CCO")


# ---------------------------------------------------------------------------
# analyse_monomer
# ---------------------------------------------------------------------------


class TestAnalyseMonomer:
    def test_lactone(self):
        info = analyse_monomer("O=C1CCCO1")
        assert info.valid
        assert info.group == "L"
        assert info.ring_size == 5
        assert info.n_rings == 1
        assert info.functional_groups == {"L": 1}

    def test_acyclic_is_invalid(self):
        info = analyse_monomer("CC(=O)OCC")
        assert not info.valid
        assert "not cyclic" in info.message

    def test_cyclic_without_group_is_invalid(self):
        info = analyse_monomer("C1CCCCC1")
        assert not info.valid
        assert "no supported" in info.message

    def test_priority_lactone_over_lactam(self):
        # ring containing both an ester and an amide bond
        info = analyse_monomer("O=C1CNC(=O)CO1")
        assert info.functional_groups == {"L": 1, "Lm": 1}
        assert info.group == "L"

    def test_priority_order_matches_paper(self):
        assert list(MONOMER_GROUPS)[:6] == ["L", "tnL", "tL", "dtL", "Lm", "tLm"]

    def test_thiolactam(self):
        assert analyse_monomer("S=C1CCCN1").group == "tLm"

    def test_requested_group_missing_is_invalid(self):
        info = analyse_monomer("O=C1CCCO1", group="Lm")
        assert not info.valid

    def test_fused_bicyclic_ring_size(self):
        info = analyse_monomer("O=C1OC[C@@H]2CCCC[C@H]12")
        assert info.n_rings == 2
        assert info.ring_sizes == [5, 6]
        assert info.ring_size == 5
        assert info.unassigned_stereo == []

    def test_unassigned_stereo_reported(self):
        info = analyse_monomer("O=C1OCC2CCCCC12")
        assert len(info.unassigned_stereo) == 2

    def test_bridged_lactone_ring_size(self):
        # 2-oxabicyclo[2.2.2]octan-3-one: the ester bond is in 6-membered rings
        info = analyse_monomer("O=C1OC2CCC1CC2")
        assert info.valid
        assert info.ring_size == 6


# ---------------------------------------------------------------------------
# get_isodesmic_initiator
# ---------------------------------------------------------------------------


class TestIsodesmicInitiator:
    @pytest.mark.parametrize(
        ("monomer", "expected"),
        [
            ("O=C1CCCO1", "CC(=O)OC"),  # lactone -> methyl acetate
            ("S=C1CCCO1", "CC(=S)OC"),  # thionolactone -> O-methyl thioacetate
            ("O=C1CCCS1", "CC(=O)SC"),  # thiolactone -> S-methyl thioacetate
            ("S=C1CCCS1", "CC(=S)SC"),  # dithiolactone -> methyl dithioacetate
            ("O=C1CCCN1", "CC(=O)NC"),  # lactam -> N-methylacetamide
            ("O=C1CCCN1C", "CC(=O)N(C)C"),  # N-methyl lactam
            ("O=C1OCCCO1", "COC(=O)OC"),  # carbonate -> dimethyl carbonate
            ("O=C1OCCN1", "CNC(=O)OC"),  # carbamate -> methyl N-methylcarbamate
        ],
    )
    def test_initiator(self, monomer, expected):
        assert _same_molecule(get_isodesmic_initiator(monomer), expected)

    @pytest.mark.parametrize(
        "monomer",
        ["O=C1CCCO1", "O=C1CCCS1", "S=C1CCCO1", "O=C1CCCN1", "O=C1OCCCO1"],
    )
    def test_reaction_is_isodesmic(self, monomer):
        """The ring-opened product contains two copies of the opened group."""
        info = analyse_monomer(monomer)
        initiator = get_isodesmic_initiator(monomer)
        product = MolFromSmiles(build_polymer(monomer, 1, "ROR", initiator))
        assert info.group is not None
        acyclic = MONOMER_GROUPS[info.group].replace(";R", "").replace("@", "-")
        assert len(product.GetSubstructMatches(MolFromSmarts(acyclic))) == 2

    def test_invalid_monomer_raises(self):
        with pytest.raises(ValueError, match="Cannot choose initiator"):
            get_isodesmic_initiator("CCCC")
