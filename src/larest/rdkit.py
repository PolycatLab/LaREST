"""RDKit conformer generation and xTB thermodynamic ranking for the first pipeline stage.

Generates an MMFF conformer ensemble using RDKit, optimises and aligns the
conformers, then re-ranks them by xTB free energy (G).  The best conformer
(lowest G) is taken forward to the CREST stage.
"""

from __future__ import annotations

import logging
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from rdkit.Chem.AllChem import (
    MMFFGetMoleculeForceField,  # type: ignore  # noqa: PGH003
    MMFFGetMoleculeProperties,  # type: ignore  # noqa: PGH003
)
from rdkit.Chem.MolStandardize.rdMolStandardize import StandardizeSmiles
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.rdForceFieldHelpers import MMFFOptimizeMoleculeConfs
from rdkit.Chem.rdMolAlign import AlignMolConformers
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier, MolToXYZFile, SDWriter
from rdkit.Chem.rdmolops import AddHs

from larest.chem import get_mol
from larest.constants import KCALMOL_TO_JMOL, THERMODYNAMIC_PARAMS
from larest.output import create_dir
from larest.setup import parse_command_args
from larest.xtb import parse_xtb_output

if TYPE_CHECKING:
    from rdkit.Chem.rdchem import Mol
    from rdkit.ForceField.rdForceField import MMFFMolProperties

logger = logging.getLogger(__name__)


def run_rdkit(
    smiles: str,
    dir_path: Path,
    config: dict[str, Any],
) -> dict[str, dict[str, float | None]]:
    """Generate MMFF conformers, rank by xTB free energy, and return the best result.

    Workflow:

    1. Embed ``n_conformers`` conformers using ``EmbedMultipleConfs``.
    2. Optimise all conformers with MMFF.
    3. Align conformers geometrically.
    4. Write conformers to ``<dir_path>/rdkit/conformers.sdf``.
    5. For each conformer, run xTB and parse H, S, G.
    6. Return the thermodynamic parameters of the lowest-G conformer.

    Parameters
    ----------
    smiles : str
        SMILES string of the molecule to process.
    dir_path : Path
        Root molecule output directory (e.g. ``output/Monomer/<slug>``).
    config : dict[str, Any]
        Full pipeline configuration dict.

    Returns
    -------
    dict[str, dict[str, float | None]]
        Single-entry dict ``{"rdkit": {"H": ..., "S": ..., "G": ...}}``
        containing the xTB thermodynamic parameters for the best conformer.
    """
    rdkit_dir = dir_path / "rdkit"
    create_dir(rdkit_dir)

    logger.debug("Generating conformers and computing energies using RDKit")
    rdkit_mol: Mol = AddHs(get_mol(StandardizeSmiles(smiles)))

    n_conformers: int = config["rdkit"]["n_conformers"]
    logger.debug(f"Generating {n_conformers} conformers")
    conformer_ids: Sequence[int] = EmbedMultipleConfs(
        rdkit_mol,
        n_conformers,
        useRandomCoords=True,
        randomSeed=config["rdkit"]["random_seed"],
        boxSizeMult=config["rdkit"]["conformer_box_size"],
        numThreads=config["rdkit"]["n_cores"],
    )

    logger.debug(f"Optimising the {n_conformers} conformers")
    MMFFOptimizeMoleculeConfs(
        rdkit_mol,
        numThreads=config["rdkit"]["n_cores"],
        maxIters=config["rdkit"]["mmff_iters"],
        mmffVariant=config["rdkit"]["mmff"],
    )

    logger.debug("Computing molecular properties for MMFF")
    mp: MMFFMolProperties = MMFFGetMoleculeProperties(
        rdkit_mol,
        mmffVariant=config["rdkit"]["mmff"],
        mmffVerbosity=0,
    )

    logger.debug(f"Computing energies for the {n_conformers} conformers")
    conformer_energies: list[tuple[int, float]] = sorted(
        [
            (
                conformer_id,
                MMFFGetMoleculeForceField(
                    rdkit_mol,
                    mp,
                    confId=conformer_id,
                ).CalcEnergy()
                * KCALMOL_TO_JMOL,
            )
            for conformer_id in conformer_ids
        ],
        key=lambda x: x[1],
    )

    logger.debug(f"Aligning the {n_conformers} conformers by their geometries")
    AlignMolConformers(rdkit_mol, maxIters=config["rdkit"]["align_iters"])

    sdf_file = rdkit_dir / "conformers.sdf"
    logger.debug(f"Writing conformers and their energies to {sdf_file}")
    with open(sdf_file, "w") as fstream:
        writer: SDWriter = SDWriter(fstream)
        for cid, energy in conformer_energies:
            rdkit_mol.SetIntProp("conformer_id", cid)
            rdkit_mol.SetDoubleProp("energy", energy)
            writer.write(rdkit_mol, confId=cid)
        writer.close()

    logger.debug("Computing thermodynamic parameters of conformers using xTB")
    xtb_default_args: list[str] = parse_command_args(sub_config=["xtb"], config=config)

    xtb_base_dir = dir_path / "xtb" / "rdkit"
    xtb_results: dict[str, list[float | None]] = {
        "conformer_id": [],
        **{param: [] for param in THERMODYNAMIC_PARAMS},
    }

    logger.debug(f"Getting conformer coordinates from {sdf_file}")
    with open(sdf_file, "rb") as sdfstream:
        mol_supplier: ForwardSDMolSupplier = ForwardSDMolSupplier(
            fileobj=sdfstream,
            sanitize=False,
            removeHs=False,
        )
        for conformer in mol_supplier:
            conformer_id: int = conformer.GetIntProp("conformer_id")
            conformer_xyz_file = rdkit_dir / f"conformer_{conformer_id}.xyz"

            MolToXYZFile(
                mol=rdkit_mol,
                filename=str(conformer_xyz_file),
                confId=conformer_id,
                precision=config["rdkit"]["precision"],
            )

            xtb_dir = xtb_base_dir / f"conformer_{conformer_id}"
            create_dir(xtb_dir)

            xtb_output_file = xtb_dir / f"conformer_{conformer_id}.txt"
            xtb_args: list[str] = [
                "xtb",
                str(conformer_xyz_file.absolute()),
                "--namespace",
                f"conformer_{conformer_id}",
                *xtb_default_args,
            ]

            try:
                with open(xtb_output_file, "w") as fstream:
                    subprocess.run(
                        xtb_args,
                        stdout=fstream,
                        stderr=subprocess.STDOUT,
                        cwd=xtb_dir,
                        check=True,
                    )
            except Exception:
                logger.exception(f"Failed to run xTB for conformer {conformer_id}")
                continue

            try:
                xtb_output = parse_xtb_output(
                    xtb_output_file=xtb_output_file,
                    temperature=config["xtb"]["etemp"],
                )
            except Exception:
                logger.exception(
                    f"Failed to parse xTB results for conformer {conformer_id}",
                )
                continue

            xtb_results["conformer_id"].append(conformer_id)
            for param in THERMODYNAMIC_PARAMS:
                xtb_results[param].append(xtb_output[param])

    xtb_results_file = xtb_base_dir / "results.csv"
    logger.debug(f"Writing results to {xtb_results_file}")

    xtb_results_df = pd.DataFrame(xtb_results, dtype=np.float64).sort_values("G")
    xtb_results_df.to_csv(xtb_results_file, header=True, index=False)

    if xtb_results_df.empty:
        raise RuntimeError(
            "All xTB conformer calculations failed; no results to return",
        )

    best_rdkit_conformer_results = xtb_results_df.iloc[0].to_dict()
    del best_rdkit_conformer_results["conformer_id"]
    return {"rdkit": best_rdkit_conformer_results}


def parse_best_rdkit_conformer(xtb_rdkit_results_file: Path) -> dict[str, float | None]:
    """Read the xTB results CSV and return the lowest-G conformer's parameters.

    Parameters
    ----------
    xtb_rdkit_results_file : Path
        Path to the ``results.csv`` produced by :func:`run_rdkit`, containing
        columns ``conformer_id``, ``H``, ``S``, and ``G``.

    Returns
    -------
    dict[str, float | None]
        Row dict for the conformer with the smallest ``G`` value, including the
        ``conformer_id`` column.
    """
    xtb_results_df: pd.DataFrame = pd.read_csv(
        xtb_rdkit_results_file,
        header=0,
        index_col=False,
        dtype=np.float64,
    ).sort_values("G", ascending=True)

    return xtb_results_df.iloc[0].to_dict()
