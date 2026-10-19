"""RDKit-based chemistry utilities for monomer analysis and polymer construction.

This module provides functions to parse SMILES strings into RDKit molecule
objects, check that monomers are suitable for ring-opening (the monomer
"sanity check"), identify the bond broken on ring opening, choose an isodesmic
initiator, and assemble polymer chains by iteratively zipping monomer (and
optionally initiator) units using RDKit's ``molzip`` functionality.

Supported reaction types
------------------------
ROR (Ring-Opening polymerization Reaction)
    An initiator caps the chain ends: its nucleophilic atom bonds to the acyl
    end of the chain and its remaining fragment (H for water and alcohols, or
    the acyl group for esters) caps the other end.
RER (Ring Expansion Reaction)
    No initiator; the chain is closed into a cyclic oligomer.
"""

from __future__ import annotations

import logging
from typing import Literal

from rdkit.Chem import FindPotentialStereo, StereoSpecified, StereoType
from rdkit.Chem.rdchem import Atom, BondType, EditableMol, Mol, RingInfo
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles, MolToSmiles
from rdkit.Chem.rdmolops import (
    AddHs,
    MolzipLabel,
    MolzipParams,
    RemoveHs,
    molzip,
)

from larest.constants import INITIATOR_GROUPS, MONOMER_GROUPS
from larest.data import MonomerInfo

logger = logging.getLogger(__name__)

AUTO_INITIATOR = "auto"


def get_mol(smiles: str) -> Mol:
    """Get an RDKit Mol object from a SMILES string

    Parameters
    ----------
    smiles : str
        Input SMILES string of the molecule

    Returns
    -------
    Mol
        RDKit Mol object with the specified SMILES string

    """
    mol: Mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(
            f"Failed to create RDKit Mol object from SMILES: {smiles}",
        )
    logger.debug(f"Created RDKit Mol object from SMILES: {smiles}")
    return mol


def _match_mapped_bonds(mol: Mol, smarts: str) -> list[tuple[int, int]]:
    """Return the (map 1, map 2) atom index pairs of all matches of *smarts*."""
    pattern: Mol = MolFromSmarts(smarts)
    map_positions: dict[int, int] = {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in pattern.GetAtoms()
        if atom.GetAtomMapNum()
    }
    pairs: list[tuple[int, int]] = []
    for match in mol.GetSubstructMatches(pattern, uniquify=False):
        pair = (match[map_positions[1]], match[map_positions[2]])
        if pair not in pairs:
            pairs.append(pair)
    return pairs


def find_breakable_bonds(mol: Mol) -> list[tuple[str, int, int]]:
    """Find all ring bonds that can be cleaved on ring opening.

    Parameters
    ----------
    mol : Mol
        RDKit molecule of the monomer.

    Returns
    -------
    list[tuple[str, int, int]]
        ``(group, acyl_atom_idx, heteroatom_idx)`` for every breakable bond, in
        order of decreasing priority (see
        :data:`~larest.constants.MONOMER_GROUPS`).  Each bond is listed once,
        under its highest priority group.
    """
    bonds: list[tuple[str, int, int]] = []
    seen: set[frozenset[int]] = set()
    for group, smarts in MONOMER_GROUPS.items():
        for acyl_idx, hetero_idx in _match_mapped_bonds(mol, smarts):
            key = frozenset((acyl_idx, hetero_idx))
            if key not in seen:
                seen.add(key)
                bonds.append((group, acyl_idx, hetero_idx))
    return bonds


def _unassigned_stereo(mol: Mol) -> list[str]:
    """Describe stereocentres and stereo double bonds with no assigned configuration."""
    unassigned: list[str] = []
    for info in FindPotentialStereo(mol):
        if info.specified != StereoSpecified.Unspecified:
            continue
        if info.type == StereoType.Atom_Tetrahedral:
            atom = mol.GetAtomWithIdx(info.centeredOn)
            unassigned.append(f"atom {info.centeredOn} ({atom.GetSymbol()})")
        elif info.type == StereoType.Bond_Double:
            bond = mol.GetBondWithIdx(info.centeredOn)
            unassigned.append(
                f"double bond {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}",
            )
    return unassigned


def analyse_monomer(smiles: str, group: str | None = None) -> MonomerInfo:
    """Check that a monomer is suitable for ring-opening and find the bond to break.

    The monomer is checked to be cyclic and to contain at least one supported
    ring-opening functional group.  If several breakable bonds are present,
    the one belonging to the highest priority group is selected (lactone >
    thionolactone > thiolactone > dithiolactone > lactam > thiolactam >
    carbonates and derivatives), unless *group* is given.  Stereocentres and
    stereo double bonds without an assigned configuration are also reported,
    since conformer generation would otherwise pick a random configuration.

    Parameters
    ----------
    smiles : str
        SMILES string of the monomer.
    group : str | None
        Functional group to open (a key of
        :data:`~larest.constants.MONOMER_GROUPS`), overriding the priority
        ranking.

    Returns
    -------
    MonomerInfo
        Summary of the monomer analysis.  ``valid`` is ``False`` (with the
        reason in ``message``) if the monomer cannot be ring-opened.

    Raises
    ------
    ValueError
        If the SMILES string is invalid.
    """
    info = _analyse_monomer(smiles, group=group)
    if info.valid:
        logger.info(
            f"Monomer {smiles}: {info.n_rings} ring(s) {info.ring_sizes}, "
            f"functional groups {info.functional_groups}, opening {info.group} "
            f"(ring size {info.ring_size})",
        )
    else:
        logger.warning(f"Monomer {smiles} is unsuitable: {info.message}")
    if info.unassigned_stereo:
        logger.warning(
            f"Monomer {smiles} has unassigned stereochemistry "
            f"({', '.join(info.unassigned_stereo)}); conformers may mix stereoisomers. "
            "Specify the stereochemistry in the SMILES to avoid this.",
        )
    return info


def _analyse_monomer(smiles: str, group: str | None = None) -> MonomerInfo:
    """Implementation of :func:`analyse_monomer`, without logging."""
    mol = get_mol(smiles)
    ring_info: RingInfo = mol.GetRingInfo()
    breakable_bonds = find_breakable_bonds(mol)

    functional_groups: dict[str, int] = {}
    for fg, _, _ in breakable_bonds:
        functional_groups[fg] = functional_groups.get(fg, 0) + 1

    info = MonomerInfo(
        smiles=smiles,
        n_rings=ring_info.NumRings(),
        ring_sizes=sorted(len(ring) for ring in ring_info.AtomRings()),
        functional_groups=functional_groups,
        unassigned_stereo=_unassigned_stereo(mol),
    )

    if info.n_rings == 0:
        info.message = "monomer is not cyclic"
    elif not breakable_bonds:
        info.message = "no supported ring-opening functional group found"
    elif group is not None and group not in functional_groups:
        info.message = (
            f"requested functional group {group} not found "
            f"(found: {', '.join(functional_groups)})"
        )
    else:
        selected = next(
            bond for bond in breakable_bonds if group is None or bond[0] == group
        )
        info.group = selected[0]
        info.bond = (selected[1], selected[2])
        info.ring_size = ring_info.MinBondRingSize(
            mol.GetBondBetweenAtoms(*info.bond).GetIdx(),
        )
        info.valid = True
        if len(breakable_bonds) > 1:
            info.message = (
                f"{len(breakable_bonds)} breakable bonds found; "
                f"opening the highest priority {info.group} bond"
            )
    return info


def get_ring_size(smiles: str) -> int | None:
    """Get the size of the ring opened for the monomer with the specified SMILES

    Parameters
    ----------
    smiles : str
        Input SMILES string of the molecule

    Returns
    -------
    int | None
        Size of the smallest ring containing the breakable bond, or None if no
        breakable bond is found

    """
    try:
        info = _analyse_monomer(smiles)
    except ValueError:
        logger.exception(f"Failed to get ring size for SMILES: {smiles}")
        raise
    if info.ring_size is None:
        logger.warning(f"Failed to detect ring size for SMILES: {smiles}")
    return info.ring_size


def get_isodesmic_initiator(smiles: str, group: str | None = None) -> str:
    """Get the initiator giving an isodesmic ring-opening reaction for a monomer.

    The initiator is the acyclic analogue of the ring-opening functional group,
    with each ring neighbour of the acyl carbon replaced by a methyl group, e.g.
    methyl acetate for lactones, S-methyl thioacetate for thiolactones,
    N-methylacetamide for lactams and dimethyl carbonate for carbonates.  The
    types of bonds broken and formed are then conserved on ring opening.

    Parameters
    ----------
    smiles : str
        SMILES string of the monomer.
    group : str | None
        Functional group to open, overriding the priority ranking.

    Returns
    -------
    str
        Canonical SMILES string of the isodesmic initiator.

    Raises
    ------
    ValueError
        If no breakable bond is found in the monomer.
    """
    info = _analyse_monomer(smiles, group=group)
    if not info.valid or info.bond is None:
        raise ValueError(f"Cannot choose initiator for {smiles}: {info.message}")

    mol = get_mol(smiles)
    acyl_idx, hetero_idx = info.bond
    acyl: Atom = mol.GetAtomWithIdx(acyl_idx)
    hetero: Atom = mol.GetAtomWithIdx(hetero_idx)

    double_bonded = ""
    alpha = ""
    for bond in acyl.GetBonds():
        neighbour = bond.GetOtherAtom(acyl)
        if neighbour.GetIdx() == hetero_idx:
            continue
        if bond.GetBondType() == BondType.DOUBLE:
            double_bonded = neighbour.GetSymbol()
        else:
            alpha = _methyl_capped(neighbour, leading=True)

    initiator = f"{alpha}C(={double_bonded}){_methyl_capped(hetero, leading=False)}"
    initiator = MolToSmiles(get_mol(initiator))
    logger.info(f"Isodesmic initiator for {smiles} ({info.group}): {initiator}")
    return initiator


def _methyl_capped(atom: Atom, *, leading: bool) -> str:
    """SMILES fragment for a ring atom with its ring substituents replaced by methyl."""
    if atom.GetAtomicNum() == 6:
        return "C"
    symbol = atom.GetSymbol()
    # keep N-H for secondary amides/carbamates, otherwise N,N-dimethyl
    if symbol == "N" and atom.GetTotalNumHs() == 0:
        return "CN(C)" if leading else "N(C)C"
    return f"C{symbol}" if leading else f"{symbol}C"


def get_polymer_unit(
    smiles: str,
    mol_type: Literal["monomer", "initiator"],
    front_dummy: str,
    back_dummy: str,
    group: str | None = None,
) -> Mol:
    """Get singular polymer build unit, constructed from the monomer/initiator

    For a monomer, the breakable ring bond is cleaved, the front dummy is
    attached to the acyl carbon and the back dummy to the ring heteroatom.  For
    an initiator, the front dummy is attached to the nucleophilic atom (which
    bonds to the acyl end of the chain) and the back dummy to the fragment that
    caps the heteroatom end of the chain.

    Parameters
    ----------
    smiles : str
        Input SMILES string of the monomer/initiator
    mol_type : Literal["monomer", "initiator"]
        Whether the build unit is from the monomer or initiator
    front_dummy : str
        Front dummy atom used to attach to the polymer chain
    back_dummy : str
        Back dummy atom used to attach to the polymer chain
    group : str | None
        Monomer functional group to open, overriding the priority ranking

    Returns
    -------
    Mol
        RDKIT Mol object for the polymer build unit

    """
    logger.debug(f"Getting {mol_type} unit for polymer construction (SMILES: {smiles})")
    logger.debug(f"Front dummy: {front_dummy}, Back dummy {back_dummy}")
    match mol_type:
        case "monomer":
            mol: Mol = get_mol(smiles)
            info = _analyse_monomer(smiles, group=group)
            if not info.valid or info.bond is None:
                raise ValueError(
                    f"No functional group atom ids found for {mol_type} SMILES: "
                    f"{smiles} ({info.message})",
                )
            fg_atom_idx: tuple[int, int] = info.bond
        case "initiator":
            # need to add Hs to break O-H bond
            mol: Mol = AddHs(mol=get_mol(smiles), explicitOnly=False)
            for fg, fg_smarts in INITIATOR_GROUPS.items():
                pairs = _match_mapped_bonds(mol, fg_smarts)
                if pairs:
                    fg_atom_idx = pairs[0]
                    logger.debug(f"Functional group detected: {fg}")
                    break
            else:
                raise ValueError(
                    f"No functional group atom ids found for {mol_type} SMILES: {smiles}",
                )
    logger.debug(f"Functional group atom ids: {fg_atom_idx}")

    front_dummy_atom: Atom = MolFromSmiles(f"[{front_dummy}]").GetAtomWithIdx(0)
    back_dummy_atom: Atom = MolFromSmiles(f"[{back_dummy}]").GetAtomWithIdx(0)

    # creating an editable Mol
    emol: EditableMol = EditableMol(mol)
    # breaking functional group bond
    emol.RemoveBond(*fg_atom_idx)
    # adding front dummy to first fg group atom
    emol.AddBond(
        beginAtomIdx=fg_atom_idx[0],
        endAtomIdx=emol.AddAtom(front_dummy_atom),
        order=BondType.SINGLE,
    )
    # adding back dummy to second fg group atom
    emol.AddBond(
        beginAtomIdx=fg_atom_idx[1],
        endAtomIdx=emol.AddAtom(back_dummy_atom),
        order=BondType.SINGLE,
    )

    logger.debug(f"Created polymer unit {MolToSmiles(emol.GetMol())}")
    return emol.GetMol()


def build_polymer(
    monomer_smiles: str,
    polymer_length: int,
    reaction_type: Literal["ROR", "RER"],
    initiator_smiles: str | None = None,
    group: str | None = None,
) -> str:
    """Build a polymer chain SMILES string from a monomer and reaction type.

    Iteratively assembles the chain by zipping monomer units together using
    RDKit's ``molzip``, then caps the chain ends with an initiator unit (ROR)
    or closes the chain into a ring (RER).

    Parameters
    ----------
    monomer_smiles : str
        SMILES string of the monomer.
    polymer_length : int
        Number of monomer repeat units.  Must be >= 2 for RER and >= 1 for ROR.
    reaction_type : {"ROR", "RER"}
        ``"ROR"`` caps the chain with *initiator_smiles*; ``"RER"`` builds a
        cyclic oligomer.
    initiator_smiles : str | None
        SMILES string of the initiator.  Required when *reaction_type* is
        ``"ROR"``.
    group : str | None
        Monomer functional group to open, overriding the priority ranking.

    Returns
    -------
    str
        Canonical SMILES string of the assembled polymer chain.

    Raises
    ------
    ValueError
        If *polymer_length* is out of range for the given *reaction_type*, if
        no initiator is given for ROR, or if a monomer/initiator unit cannot be
        constructed from the provided SMILES.
    """
    if polymer_length <= 1 and reaction_type == "RER":
        raise ValueError(
            f"Please specify a polymer length > 1 for RER reaction (current length: {polymer_length})",
        )
    if polymer_length < 1 and reaction_type == "ROR":
        raise ValueError(
            f"Please specify a polymer length >= 1 for ROR reaction (current length: {polymer_length})",
        )
    if reaction_type == "ROR" and not initiator_smiles:
        raise ValueError("An initiator SMILES is required for ROR reactions")

    logger.debug(f"Building {reaction_type} polymer")
    logger.debug(f"Monomer smiles: {monomer_smiles}, length: {polymer_length}")
    if reaction_type == "ROR":
        logger.debug(f"Initiator smiles: {initiator_smiles}")

    mol_zip_params = MolzipParams()
    mol_zip_params.label = MolzipLabel.AtomType

    try:
        polymer = get_polymer_unit(
            smiles=monomer_smiles,
            mol_type="monomer",
            front_dummy="Xe",
            back_dummy="Y",
            group=group,
        )
    except Exception:
        logger.exception(f"Failed to create monomer unit from SMILES: {monomer_smiles}")
        raise

    monomer_units: int = 1
    front_dummy: str = "Y"
    back_dummy: str = "Zr"

    if reaction_type == "RER":
        polymer_length -= 1
    while monomer_units < polymer_length:
        repeating_unit = get_polymer_unit(
            smiles=monomer_smiles,
            mol_type="monomer",
            front_dummy=front_dummy,
            back_dummy=back_dummy,
            group=group,
        )
        mol_zip_params.setAtomSymbols([front_dummy])
        try:
            polymer = molzip(polymer, repeating_unit, mol_zip_params)
        except Exception:
            logger.exception(
                f"Failed to zip polymer units together: {MolToSmiles(polymer)} {MolToSmiles(repeating_unit)}",
            )
            raise
        else:
            front_dummy, back_dummy = back_dummy, front_dummy
            monomer_units += 1
            logger.debug(
                f"Polymer chain with length={monomer_units}): {MolToSmiles(polymer)}",
            )

    # the initiator caps the chain for ROR; a final monomer closes the ring for RER
    terminal_smiles = (
        initiator_smiles
        if reaction_type == "ROR" and initiator_smiles
        else monomer_smiles
    )
    try:
        match reaction_type:
            case "ROR":
                terminal_unit = get_polymer_unit(
                    smiles=terminal_smiles,
                    mol_type="initiator",
                    front_dummy="Xe",
                    back_dummy=front_dummy,
                )
            case "RER":
                terminal_unit = get_polymer_unit(
                    smiles=monomer_smiles,
                    mol_type="monomer",
                    front_dummy=front_dummy,
                    back_dummy="Xe",
                    group=group,
                )
    except Exception:
        logger.exception(
            f"Failed to create terminal group for {reaction_type} reaction",
        )
        raise

    mol_zip_params.setAtomSymbols([front_dummy, "Xe"])
    try:
        match reaction_type:
            case "ROR":
                # Need to remove Hs following initiator addition
                polymer = RemoveHs(
                    mol=molzip(polymer, terminal_unit, mol_zip_params),
                    implicitOnly=False,
                    sanitize=True,
                )
            case "RER":
                polymer = molzip(polymer, terminal_unit, mol_zip_params)
    except Exception:
        logger.exception(
            "Failed to zip polymer and terminal unit together: "
            f"{MolToSmiles(polymer)} {MolToSmiles(terminal_unit)}",
        )
        raise

    logger.info(
        f"Finished building {reaction_type} polymer: {MolToSmiles(polymer)}",
    )
    return MolToSmiles(polymer)
