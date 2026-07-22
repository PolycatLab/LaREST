"""Checkpoint and result-restoration utilities for the LaREST pipeline.

At the start of each molecule's run, :func:`restore_results` walks through the
expected output files in pipeline order and returns both the accumulated results
and the first stage that has not yet been completed.  This allows interrupted
runs to resume from exactly where they left off without re-executing earlier
stages.

:func:`apply_entropy_correction` combines CENSO ``censo_refinement`` enthalpy /
free-energy values with the conformational entropy obtained from a CREST
entropy run, producing the ``censo_corrected`` section stored alongside the
other thermodynamic sections.
"""

import json
import logging
from collections.abc import Callable
from enum import IntEnum
from pathlib import Path

from larest.constants import PIPELINE_SECTIONS, THERMODYNAMIC_PARAMS
from larest.rdkit import parse_best_rdkit_conformer

logger = logging.getLogger(__name__)


class PipelineStage(IntEnum):
    """Ordered enumeration of pipeline stages used for checkpointing.

    Stages are assigned integer values so that simple comparisons (``<=``)
    can determine whether a given stage still needs to be executed.
    """

    RDKIT = 1
    CREST_CONFGEN = 2
    CENSO = 3
    CREST_ENTROPY = 4
    FINISH = 5


def restore_results(
    dir_path: Path,
    temperature: float = 298.15,
) -> tuple[dict[str, dict[str, float | None]], PipelineStage]:
    """Load previously completed pipeline results and determine the next stage to run.

    Inspects the expected output files for each pipeline stage in order.  The
    first missing or unreadable file determines the stage that must be (re-)run.
    All successfully loaded results are merged and returned alongside that stage.

    Parameters
    ----------
    dir_path : Path
        Root output directory for the molecule (e.g.
        ``output/Monomer/<slug>``).
    temperature : float
        Temperature in Kelvin used to recompute the corrected free energy when
        restoring a CREST entropy checkpoint (see :func:`apply_entropy_correction`).

    Returns
    -------
    results : dict[str, dict[str, float | None]]
        Mapping of section name to thermodynamic parameters populated from any
        checkpoint files that were found.  Parameters that have not yet been
        computed remain ``None``.
    stage : PipelineStage
        The earliest pipeline stage whose results were not found.  The caller
        should (re-)run this stage and all subsequent ones.
    """
    results: dict[str, dict[str, float | None]] = {
        section: dict.fromkeys(THERMODYNAMIC_PARAMS, None)
        for section in PIPELINE_SECTIONS
    }

    if not _load_stage(
        path=dir_path / "xtb" / "rdkit" / "results.csv",
        load_fn=lambda path: _parse_rdkit(results, path),
        label="rdkit",
    ):
        return results, PipelineStage.RDKIT

    if not _load_stage(
        path=dir_path / "xtb" / "crest" / "results.json",
        load_fn=lambda path: results.update({"crest": json.loads(path.read_text())}),
        label="crest_confgen",
    ):
        return results, PipelineStage.CREST_CONFGEN

    if not _load_stage(
        path=dir_path / "censo" / "results.json",
        load_fn=lambda path: results.update(json.loads(path.read_text())),
        label="censo",
    ):
        return results, PipelineStage.CENSO

    if not _load_stage(
        path=dir_path / "crest_entropy" / "results.json",
        load_fn=lambda path: _parse_crest_entropy(results, path, temperature),
        label="crest_entropy",
    ):
        return results, PipelineStage.CREST_ENTROPY

    _load_stage(
        path=dir_path / "results.json",
        load_fn=lambda path: results.update(json.loads(path.read_text())),
        label="final",
    )
    return results, PipelineStage.FINISH


def _load_stage(path: Path, load_fn: Callable[[Path], None], label: str) -> bool:
    """Attempt to load a single checkpoint stage from disk.

    Parameters
    ----------
    path : Path
        Expected location of the checkpoint file.
    load_fn : Callable[[Path], None]
        Function that reads *path* and mutates the shared results dict in
        the enclosing scope.
    label : str
        Human-readable name of the stage, used only for log messages.

    Returns
    -------
    bool
        ``True`` if the file was found and loaded without error, ``False``
        if the file was absent or ``load_fn`` raised an exception.
    """
    if not path.exists():
        logger.info(f"No pre-existing {label} results at {path}, running pipeline...")
        return False
    try:
        load_fn(path)
    except Exception:
        logger.exception(f"Failed to load pre-existing {label} results from {path}")
        return False
    logger.info(f"Loaded pre-existing {label} results from {path}")
    return True


def _parse_rdkit(results: dict, path: Path) -> None:
    """Load the best RDKit/xTB conformer results into *results*.

    Parameters
    ----------
    results : dict
        Shared results dict to update in-place with the ``"rdkit"`` section.
    path : Path
        Path to the xTB results CSV produced during the RDKit stage.
    """
    best = parse_best_rdkit_conformer(path)
    del best["conformer_id"]
    results.update({"rdkit": best})


def apply_entropy_correction(
    refinement_results: dict[str, float | None],
    crest_entropy_results: dict[str, float | None],
    temperature: float,
) -> dict[str, dict[str, float | None]]:
    """Apply a CREST conformational entropy correction to CENSO refinement results.

    The CENSO ``censo_refinement`` stage yields the single-conformer
    thermostatistical (RRHO / vibrational) entropy.  A CREST entropy run adds
    the conformational entropy contribution — ``S_total`` (``Sconf + δSrrho``,
    where ``δSrrho`` is the small ensemble-averaging correction, *not* a second
    full RRHO term) — which is a distinct contribution and is therefore *added*
    to, not substituted for, the CENSO entropy.  The free energy is then
    recomputed from the corrected entropy as ``G = H - T·S`` so that H, S, and G
    remain mutually consistent.

    Parameters
    ----------
    refinement_results : dict[str, float | None]
        Thermodynamic parameters from the ``censo_refinement`` CENSO stage, keyed
        by ``"H"``, ``"S"``, and ``"G"``.
    crest_entropy_results : dict[str, float | None]
        Entropy parameters from the CREST entropy run, keyed by
        ``"S_conf"``, ``"S_rrho"``, and ``"S_total"``.
    temperature : float
        Temperature in Kelvin at which the CENSO values were evaluated; used to
        recompute the corrected free energy.

    Returns
    -------
    dict[str, dict[str, float | None]]
        A single-entry dict ``{"censo_corrected": {...}}`` where ``"S"`` is the
        sum of the CENSO refinement entropy and the CREST total conformational
        entropy, and ``"G"`` has been recomputed as ``H - T·S``.

    Raises
    ------
    ValueError
        If any of ``refinement_results["H"]``, ``refinement_results["S"]``, or
        ``crest_entropy_results["S_total"]`` is ``None``.
    """
    h = refinement_results["H"]
    s = refinement_results["S"]
    s_total = crest_entropy_results["S_total"]
    if h is None or s is None or s_total is None:
        raise ValueError("Failed to apply CREST entropy correction to CENSO results")
    corrected = refinement_results.copy()
    corrected["S"] = s + s_total
    corrected["G"] = h - temperature * corrected["S"]
    return {"censo_corrected": corrected}


def _parse_crest_entropy(results: dict, path: Path, temperature: float) -> None:
    """Load CREST entropy results and apply the entropy correction to *results*.

    Parameters
    ----------
    results : dict
        Shared results dict; must already contain a ``"censo_refinement"`` section.
        Updated in-place with a ``"censo_corrected"`` section.
    path : Path
        Path to the CREST entropy ``results.json`` file.
    temperature : float
        Temperature in Kelvin used to recompute the corrected free energy.
    """
    crest_entropy_results: dict[str, float | None] = json.loads(path.read_text())
    results.update(
        apply_entropy_correction(
            results["censo_refinement"],
            crest_entropy_results,
            temperature,
        ),
    )
