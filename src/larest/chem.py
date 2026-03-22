"""RDKit-based chemistry utilities for polymer construction and ring analysis.

This module provides functions to parse SMILES strings into RDKit molecule
objects, identify ring-opening functional groups, and assemble polymer chains
by iteratively zipping monomer (and optionally initiator) units using RDKit's
``molzip`` functionality.

Supported reaction types
------------------------
ROR (Ring-Opening polymerization Reaction)
    An initiator alcohol caps one chain end; the other end is the growing chain.
RER (Ring Equilibrium Reaction)
    No initiator; the chain is terminated by an additional monomer unit on each
    end.
"""

from __future__ import annotations

import logging
from typing import Any, Literal

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

logger = logging.getLogger(__name__)


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


def get_ring_size(smiles: str) -> int | None:
    """Get the ring size of the monomer with the specified SMILES

    Parameters
    ----------
    smiles : str
        Input SMILES string of the molecule

    Returns
    -------
    int | None
        Funtional group ring size of the monomer, or None is ring size cannot be determined

    """
    try:
        mol: Mol = get_mol(smiles)
    except ValueError:
        logger.exception(f"Failed to get ring size for SMILES: {smiles}")
        raise

    ring_info: RingInfo = mol.GetRingInfo()
    functional_groups: list[Mol] = [
        MolFromSmarts(fg_smarts) for fg_smarts in MONOMER_GROUPS.values()
    ]
    for fg_mol in functional_groups:
        substruct_atom_idx: tuple[int] = mol.GetSubstructMatch(fg_mol)
        if substruct_atom_idx:
            for atom_id in substruct_atom_idx:
                substruct_atom: Atom = mol.GetAtomWithIdx(atom_id)
                if (substruct_atom.GetAtomicNum() != 6) and (substruct_atom.IsInRing()):
                    ring_size: int = ring_info.MinAtomRingSize(atom_id)
                    logger.debug(f"Found ring size {ring_size} for SMILES: {smiles}")
                    return ring_size
    logger.warning(f"Failed to detect ring size for SMILES: {smiles}")
    return None


def get_polymer_unit(
    smiles: str,
    mol_type: Literal["monomer", "initiator"],
    front_dummy: str,
    back_dummy: str,
) -> Mol:
    """Get singular polymer build unit, constructed from the monomer/initiator

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

    Returns
    -------
    Mol
        RDKIT Mol object for the polymer build unit

    """
    # NOTE: output for molecules with >1 ring-opening functional group
    # is deterministic but not yet customisable

    # convert monomer/initiator smiles and front/back dummies to RDKit Mol/Atoms
    logger.debug(f"Getting {mol_type} unit for polymer construction (SMILES: {smiles})")
    logger.debug(f"Front dummy: {front_dummy}, Back dummy {back_dummy}")
    match mol_type:
        case "monomer":
            functional_groups: list[Mol] = [
                MolFromSmarts(fg_smarts) for fg_smarts in MONOMER_GROUPS.values()
            ]
            mol: Mol = MolFromSmiles(smiles)
        case "initiator":
            functional_groups: list[Mol] = [
                MolFromSmarts(fg_smarts) for fg_smarts in INITIATOR_GROUPS.values()
            ]
            # need to add Hs to break O-H bond
            mol: Mol = AddHs(
                mol=MolFromSmiles(smiles),
                explicitOnly=False,
            )
    front_dummy_atom: Atom = MolFromSmiles(f"[{front_dummy}]").GetAtomWithIdx(0)
    back_dummy_atom: Atom = MolFromSmiles(f"[{back_dummy}]").GetAtomWithIdx(0)

    # obtaining the atom ids of the functional group in the monomer/initiator
    fg_atom_idx = []
    for fg_mol in functional_groups:
        substruct_atom_idx: tuple[int] = mol.GetSubstructMatch(fg_mol)
        if substruct_atom_idx:
            for atom_id in substruct_atom_idx:
                substruct_atom: Atom = mol.GetAtomWithIdx(atom_id)
                match mol_type:
                    case "monomer":
                        # break bond between carbonyl C and electrophilic neighbour
                        if substruct_atom.IsInRing() and (
                            (
                                BondType.DOUBLE
                                in [
                                    bond.GetBondType()
                                    for bond in substruct_atom.GetBonds()
                                ]
                            )
                            or (substruct_atom.GetAtomicNum() != 6)
                        ):
                            fg_atom_idx.append(atom_id)
                    case "initiator":
                        # break O-H bond
                        fg_atom_idx.append(atom_id)
                if len(fg_atom_idx) == 2:
                    break
        if len(fg_atom_idx) != 0:
            if len(fg_atom_idx) != 2:
                raise ValueError(
                    f"Trying to break bond between >2 atoms for functional group {MolToSmiles(fg_mol)}",
                )
            logger.debug(f"Functional group detected: {MolToSmiles(fg_mol)}")
            break

    # no functional group detected in monomer/initiator
    if len(fg_atom_idx) == 0:
        raise ValueError(
            f"No functional group atom ids found for {mol_type} SMILES: {smiles}",
        )
    logger.debug(f"Functional group atom ids: {fg_atom_idx}")

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
    config: dict[str, Any],
) -> str:
    """Build a polymer chain SMILES string from a monomer and reaction type.

    Iteratively assembles the chain by zipping monomer units together using
    RDKit's ``molzip``, then caps the terminal end with either an initiator
    unit (ROR) or an additional monomer unit (RER).

    Parameters
    ----------
    monomer_smiles : str
        SMILES string of the lactone monomer.
    polymer_length : int
        Number of monomer repeat units.  Must be >= 2 for RER and >= 1 for ROR.
    reaction_type : {"ROR", "RER"}
        ``"ROR"`` uses the initiator SMILES from ``config["reaction"]["initiator"]``
        to cap the chain; ``"RER"`` uses an additional monomer unit instead.
    config : dict[str, Any]
        Full pipeline configuration dict.  The ``config["reaction"]["initiator"]``
        key is required when *reaction_type* is ``"ROR"``.

    Returns
    -------
    str
        Canonical SMILES string of the assembled polymer chain.

    Raises
    ------
    ValueError
        If *polymer_length* is out of range for the given *reaction_type*, or
        if a monomer/initiator unit cannot be constructed from the provided
        SMILES.
    """
    if polymer_length <= 1 and reaction_type == "RER":
        raise ValueError(
            f"Please specify a polymer length > 1 for RER reaction (current length: {polymer_length})",
        )
    if polymer_length < 1 and reaction_type == "ROR":
        raise ValueError(
            f"Please specify a polymer length >= 1 for ROR reaction (current length: {polymer_length})",
        )

    logger.debug(f"Building {reaction_type} polymer")
    logger.debug(f"Monomer smiles: {monomer_smiles}, length: {polymer_length}")
    if reaction_type == "ROR":
        logger.debug(f"Initiator smiles: {config['reaction']['initiator']}")

    mol_zip_params = MolzipParams()
    mol_zip_params.label = MolzipLabel.AtomType

    try:
        polymer = get_polymer_unit(
            smiles=monomer_smiles,
            mol_type="monomer",
            front_dummy="Xe",
            back_dummy="Y",
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

    try:
        match reaction_type:
            case "ROR":
                terminal_unit = get_polymer_unit(
                    smiles=config["reaction"]["initiator"],
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
