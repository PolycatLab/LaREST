"""Dataclasses representing molecules and their computed results in the LaREST pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class Monomer:
    """A cyclic monomer molecule.

    Attributes
    ----------
    smiles : str
        SMILES string of the monomer.
    """

    smiles: str


@dataclass
class MonomerInfo:
    """Result of the monomer sanity check (see :func:`larest.chem.analyse_monomer`).

    Attributes
    ----------
    smiles : str
        SMILES string of the monomer.
    n_rings : int
        Number of rings (smallest set of smallest rings).
    ring_sizes : list[int]
        Sizes of all rings.
    functional_groups : dict[str, int]
        Number of breakable bonds found for each ring-opening functional group.
    unassigned_stereo : list[str]
        Stereocentres / stereo double bonds without an assigned configuration.
    group : str | None
        Functional group selected for ring opening.
    bond : tuple[int, int] | None
        Atom indices (acyl carbon, ring heteroatom) of the bond broken.
    ring_size : int | None
        Size of the smallest ring containing the broken bond.
    valid : bool
        Whether the monomer can be ring-opened.
    message : str
        Reason the monomer is invalid, or other notes.
    """

    smiles: str
    n_rings: int = 0
    ring_sizes: list[int] = field(default_factory=list)
    functional_groups: dict[str, int] = field(default_factory=dict)
    unassigned_stereo: list[str] = field(default_factory=list)
    group: str | None = None
    bond: tuple[int, int] | None = None
    ring_size: int | None = None
    valid: bool = False
    message: str = ""


@dataclass
class Initiator:
    """An initiator (water, alcohol, ester, ...) for ring-opening reactions.

    Attributes
    ----------
    smiles : str
        SMILES string of the initiator.
    """

    smiles: str


@dataclass
class Polymer:
    """A polymer chain at a given repeat unit length.

    Attributes
    ----------
    smiles : str
        SMILES string of the polymer chain.
    monomer_smiles : str
        SMILES string of the monomer used to build the chain.
    length : int
        Number of monomer repeat units in the chain.
    initiator_smiles : str | None
        SMILES string of the initiator capping the chain (ROR only).
    """

    smiles: str
    monomer_smiles: str
    length: int
    initiator_smiles: str | None = None


@dataclass
class MolResults:
    """Thermodynamic results for a single molecule across pipeline sections.

    Attributes
    ----------
    smiles : str
        SMILES string of the molecule.
    sections : dict[str, dict[str, float | None]]
        Mapping of pipeline section name (e.g. ``"rdkit"``, ``"censo_refinement"``)
        to a dict of thermodynamic parameters (``"H"``, ``"S"``, ``"G"``).
        Values are ``None`` where a parameter could not be computed.
    dir_path : Path | None
        Output directory of the molecule.
    """

    smiles: str
    sections: dict[str, dict[str, float | None]] = field(default_factory=dict)
    dir_path: Path | None = None
