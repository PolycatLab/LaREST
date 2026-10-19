"""Migrate output directories from older LaREST versions to the current naming.

Older versions of LaREST named molecule directories with :func:`slugify` only
(``Monomer/<slug>``, ``Initiator/<slug>`` and ``Polymer/<monomer slug>_<n>``).
Because :func:`slugify` drops characters such as ``=`` and ``#``, different
molecules could share a directory.  Current versions append a hash of the
SMILES string (see :func:`~larest.output.mol_dir_name`), so old output
directories are no longer found on resume.

:func:`plan_migration` works out, from the configuration used for the old run,
which old directory belongs to which new directory, and :func:`apply_migration`
renames them.  Old directories claimed by more than one molecule (legacy slug
collisions) are never moved, as their contents are ambiguous.

The ``larest-migrate`` console script (:func:`entry_point`) wraps both and
defaults to a dry run.
"""

from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from larest.chem import (
    AUTO_INITIATOR,
    analyse_monomer,
    build_polymer,
    get_isodesmic_initiator,
)
from larest.data import Initiator, Monomer, Polymer
from larest.main import get_mol_dir
from larest.output import slugify
from larest.setup import get_config

logger = logging.getLogger(__name__)

_MOL_TYPES = ("Monomer", "Polymer", "Initiator")


@dataclass(frozen=True)
class PlannedMove:
    """A directory rename from the legacy to the current naming.

    Attributes
    ----------
    old : Path
        Existing legacy directory.
    new : Path
        Directory name used by the current LaREST version.
    """

    old: Path
    new: Path


@dataclass
class MigrationPlan:
    """Moves and warnings produced by :func:`plan_migration`.

    Attributes
    ----------
    moves : list[PlannedMove]
        Directory renames that can be performed safely.
    warnings : list[str]
        Directories that were not moved, and why.
    """

    moves: list[PlannedMove]
    warnings: list[str]


def _legacy_mol_dir(mol: Monomer | Polymer | Initiator, output_dir: Path) -> Path:
    """Return the directory of a molecule as named by older LaREST versions."""
    match mol:
        case Polymer():
            return (
                output_dir / "Polymer" / f"{slugify(mol.monomer_smiles)}_{mol.length}"
            )
        case _:
            return output_dir / type(mol).__name__ / slugify(mol.smiles)


def _get_molecules(
    config: dict[str, Any],
) -> tuple[list[Monomer | Polymer | Initiator], list[str]]:
    """Return the molecules a run with *config* computes, as :func:`larest.main.main` does."""
    reaction_type: str = config["reaction"]["type"]
    if reaction_type not in ("ROR", "RER"):
        raise ValueError(f"Unknown reaction type {reaction_type!r}")
    rtype: Literal["ROR", "RER"] = "ROR" if reaction_type == "ROR" else "RER"
    initiator_setting: str | None = config["reaction"].get("initiator")
    if rtype == "ROR" and not initiator_setting:
        raise ValueError("ROR reactions require [reaction].initiator")

    molecules: list[Monomer | Polymer | Initiator] = []
    warnings: list[str] = []
    for entry in config["reaction"]["monomers"]:
        monomer_entry: dict[str, Any] = (
            {"smiles": entry, "group": None} if isinstance(entry, str) else entry
        )
        smiles: str = monomer_entry["smiles"]
        try:
            info = analyse_monomer(smiles, group=monomer_entry.get("group"))
        except ValueError:
            warnings.append(f"Skipping monomer {smiles}: invalid SMILES")
            continue
        if not info.valid:
            warnings.append(f"Skipping monomer {smiles}: {info.message}")
            continue
        molecules.append(Monomer(smiles=smiles))

        initiator_smiles = None
        if rtype == "ROR" and initiator_setting is not None:
            if initiator_setting.strip().lower() == AUTO_INITIATOR:
                initiator_smiles = get_isodesmic_initiator(smiles, group=info.group)
            else:
                initiator_smiles = initiator_setting
            molecules.append(Initiator(smiles=initiator_smiles))

        for length in config["reaction"]["lengths"]:
            try:
                polymer_smiles = build_polymer(
                    monomer_smiles=smiles,
                    polymer_length=length,
                    reaction_type=rtype,
                    initiator_smiles=initiator_smiles,
                    group=info.group,
                )
            except ValueError as err:
                warnings.append(
                    f"Skipping polymer of {smiles} (length {length}): {err}",
                )
                continue
            molecules.append(
                Polymer(
                    smiles=polymer_smiles,
                    monomer_smiles=smiles,
                    length=length,
                    initiator_smiles=initiator_smiles,
                ),
            )
    return molecules, warnings


def plan_migration(config: dict[str, Any], output_dir: Path) -> MigrationPlan:
    """Plan the renames needed to resume a legacy output directory.

    Parameters
    ----------
    config : dict[str, Any]
        Configuration of the run that produced *output_dir*.  Legacy polymer
        directories are keyed only by monomer and length, so the reaction type
        and initiator must match the old run.
    output_dir : Path
        Root output directory of the legacy run.

    Returns
    -------
    MigrationPlan
        Safe moves, plus warnings for directories that are skipped because
        they are claimed by several molecules (legacy slug collisions) or
        because the new directory already exists.

    Raises
    ------
    ValueError
        If the reaction type is unknown, or no initiator is set for ROR.
    """
    molecules, warnings = _get_molecules(config)
    if not any((output_dir / sub).is_dir() for sub in _MOL_TYPES):
        warnings.append(
            f"No {'/'.join(_MOL_TYPES)} directories in {output_dir}; "
            "is it a LaREST output directory?",
        )

    claims: dict[Path, dict[Path, str]] = defaultdict(dict)
    for mol in molecules:
        claims[_legacy_mol_dir(mol, output_dir)][get_mol_dir(mol, output_dir)] = (
            mol.smiles
        )

    moves: list[PlannedMove] = []
    for old, targets in claims.items():
        if not old.is_dir():
            logger.debug(f"No legacy directory {old}")
            continue
        if len(targets) > 1:
            smiles = ", ".join(sorted(targets.values()))
            warnings.append(
                f"Not moving {old}: shared by several molecules ({smiles}); its "
                "contents are ambiguous and must be recomputed or assigned manually",
            )
            continue
        (new,) = targets
        if new.exists():
            warnings.append(f"Not moving {old}: {new} already exists")
            continue
        moves.append(PlannedMove(old=old, new=new))
    return MigrationPlan(moves=moves, warnings=warnings)


def apply_migration(plan: MigrationPlan) -> None:
    """Perform the moves of a migration plan, never overwriting a directory.

    Parameters
    ----------
    plan : MigrationPlan
        Plan returned by :func:`plan_migration`.

    Raises
    ------
    FileExistsError
        If a target directory has appeared since the plan was made.
    """
    for move in plan.moves:
        if move.new.exists():
            raise FileExistsError(f"Refusing to overwrite {move.new}")
        move.old.rename(move.new)
        logger.info(f"Moved {move.old} -> {move.new}")


def entry_point(argv: list[str] | None = None) -> None:
    """Parse CLI arguments and migrate a legacy output directory.

    Registered as the ``larest-migrate`` console script.  Without ``--apply``
    the planned moves are only printed.

    Parameters
    ----------
    argv : list[str] | None
        Command-line arguments; defaults to ``sys.argv[1:]``.
    """
    parser = argparse.ArgumentParser(
        description="Rename output directories from older LaREST versions so runs resume",
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to the config.toml of the old run",
    )
    parser.add_argument("output", type=str, help="Output directory of the old run")
    parser.add_argument(
        "--apply",
        action="store_true",
        help="perform the moves (default: dry run)",
    )
    args = parser.parse_args(argv)

    output_dir = Path(args.output).resolve()
    if not output_dir.is_dir():
        parser.error(f"Output directory {output_dir} does not exist")
    try:
        config = get_config(Path(args.config).resolve())
        plan = plan_migration(config, output_dir)
    except Exception as err:
        print(f"Failed to plan migration: {err}")  # noqa: T201
        raise SystemExit(1) from err

    for warning in plan.warnings:
        print(f"WARNING: {warning}")  # noqa: T201
    if not plan.moves:
        print("Nothing to move")  # noqa: T201
        return
    verb = "Moving" if args.apply else "Would move"
    for move in plan.moves:
        print(  # noqa: T201
            f"{verb} {move.old.relative_to(output_dir)} -> {move.new.relative_to(output_dir)}",
        )
    if args.apply:
        apply_migration(plan)
        print(f"Moved {len(plan.moves)} directories")  # noqa: T201
    else:
        print("Dry run: re-run with --apply to perform these moves")  # noqa: T201
