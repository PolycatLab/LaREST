"""Checkpoint and result-restoration utilities for the LaREST pipeline.

At the start of each molecule's run, :func:`restore_results` walks through the
expected output files in pipeline order and returns both the accumulated results
and the first stage that has not yet been completed.  This allows interrupted
runs to resume from exactly where they left off without re-executing earlier
stages.  Each completed stage also records the settings it was run with
(:func:`write_stage_settings`); if those settings have since changed, the stage
is re-run rather than silently reusing stale results.

:func:`apply_entropy_correction` combines CENSO ``censo_refinement_min`` enthalpy /
free-energy values with the conformational entropy obtained from a CREST
entropy run, producing the ``censo_corrected`` section stored alongside the
other thermodynamic sections.
"""

import json
import logging
import math
from collections.abc import Callable, Mapping
from enum import IntEnum
from pathlib import Path
from typing import Any

from larest.censo import get_last_censo_section, parse_censo_results
from larest.constants import (
    CENSO_SECTIONS,
    MIN_SUFFIX,
    PIPELINE_SECTIONS,
    THERMODYNAMIC_PARAMS,
)
from larest.rdkit import parse_best_rdkit_conformer
from larest.xtb import parse_xtb_output

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


# Output sub-directory of each stage, where its settings stamp is written.
_STAGE_DIRS: dict[PipelineStage, str] = {
    PipelineStage.RDKIT: "rdkit",
    PipelineStage.CREST_CONFGEN: "crest_confgen",
    PipelineStage.CENSO: "censo",
    PipelineStage.CREST_ENTROPY: "crest_entropy",
}
_SETTINGS_FILE = "larest_settings.json"
# [steps] keys of the pipeline stages, in order.
_STEP_KEYS: tuple[str, ...] = ("rdkit", "crest_confgen", "censo", "crest_entropy")
# Keys that only affect performance, not results, and so never invalidate a stage.
_PERFORMANCE_KEYS: frozenset[str] = frozenset(
    {"n_cores", "parallel", "T", "maxcores", "paths", "verbose", "json"},
)


def _strip_performance_keys(settings: dict[str, Any]) -> dict[str, Any]:
    return {
        k: _strip_performance_keys(v) if isinstance(v, dict) else v
        for k, v in settings.items()
        if k not in _PERFORMANCE_KEYS
    }


def get_stage_settings(
    stage: PipelineStage,
    config: dict[str, Any],
    smiles: str,
) -> dict[str, Any]:
    """Return the settings that determine the results of a pipeline stage.

    A stage's results depend on every enabled stage before it, so the settings
    of those stages (and which stages are enabled) are included too.  A change
    upstream therefore invalidates all later stages, even if a re-run was
    interrupted before reaching them.

    Parameters
    ----------
    stage : PipelineStage
        The pipeline stage.
    config : dict[str, Any]
        Full pipeline configuration dict.
    smiles : str
        SMILES string of the molecule.

    Returns
    -------
    dict[str, Any]
        JSON-serialisable settings, excluding performance-only keys such as
        core counts.
    """
    steps = config.get("steps", {})
    crest = config.get("crest", {})
    own_settings: dict[PipelineStage, dict[str, Any]] = {
        PipelineStage.RDKIT: {"rdkit": config.get("rdkit"), "xtb": config.get("xtb")},
        PipelineStage.CREST_CONFGEN: {
            "crest.confgen": crest.get("confgen"),
            "xtb": config.get("xtb") if steps.get("xtb", True) else None,
        },
        PipelineStage.CENSO: {"censo": config.get("censo")},
        PipelineStage.CREST_ENTROPY: {"crest.entropy": crest.get("entropy")},
    }
    # which stages up to and including this one are enabled (and, from CREST
    # on, whether xTB re-ranking is used) determines this stage's input
    step_flags = {key: steps.get(key, True) for key in _STEP_KEYS[:stage]}
    if stage >= PipelineStage.CREST_CONFGEN:
        step_flags["xtb"] = steps.get("xtb", True)
    settings: dict[str, Any] = {"smiles": smiles, "steps": step_flags}
    for upstream in PipelineStage:
        if upstream > stage or upstream == PipelineStage.FINISH:
            break
        if upstream == stage or steps.get(_STEP_KEYS[upstream - 1], True):
            settings[upstream.name.lower()] = own_settings[upstream]
    return _strip_performance_keys(settings)


def write_stage_settings(
    dir_path: Path,
    stage: PipelineStage,
    config: dict[str, Any],
    smiles: str,
) -> None:
    """Record the settings a completed stage was run with in its output directory."""
    settings_file = dir_path / _STAGE_DIRS[stage] / _SETTINGS_FILE
    settings_file.parent.mkdir(parents=True, exist_ok=True)
    settings_file.write_text(
        json.dumps(get_stage_settings(stage, config, smiles), indent=4, sort_keys=True),
    )


def _settings_changed(
    dir_path: Path,
    stage: PipelineStage,
    config: dict[str, Any] | None,
    smiles: str | None,
) -> bool:
    """Return ``True`` if a stage was previously run with different settings.

    Stages without a settings stamp (from older LaREST versions) are assumed
    unchanged, with a warning, and the current settings are recorded for them.
    """
    if config is None or smiles is None:
        return False
    settings_file = dir_path / _STAGE_DIRS[stage] / _SETTINGS_FILE
    if not settings_file.exists():
        if (dir_path / _STAGE_DIRS[stage]).exists():
            logger.warning(
                f"No settings record for {stage.name.lower()} stage in {dir_path}; "
                "assuming results match the current settings",
            )
            write_stage_settings(dir_path, stage, config, smiles)
        return False
    previous = json.loads(settings_file.read_text())
    current = json.loads(json.dumps(get_stage_settings(stage, config, smiles)))
    if previous != current:
        logger.warning(
            f"Settings for the {stage.name.lower()} stage in {dir_path} have changed "
            "since it was run; re-running from this stage",
        )
        return True
    return False


def restore_results(
    dir_path: Path,
    temperature: float = 298.15,
    steps: dict[str, bool] | None = None,
    config: dict[str, Any] | None = None,
    smiles: str | None = None,
) -> tuple[dict[str, dict[str, float | None]], PipelineStage]:
    """Load previously completed pipeline results and determine the next stage to run.

    Inspects the expected output files for each *enabled* pipeline stage in
    order.  The first enabled stage with a missing or unreadable file determines
    the stage that must be (re-)run.  Stages that are disabled in *steps* are
    skipped entirely — they neither trigger a re-run nor short-circuit the
    inspection of later stages, so a completed downstream stage (e.g. CENSO) is
    still recognised as a valid checkpoint even when an upstream stage (e.g.
    ``crest_confgen``) is turned off.  All successfully loaded results are merged
    and returned alongside that stage.  A stage whose recorded settings differ
    from the current *config* is treated as missing.

    CENSO results are re-derived from CENSO's raw output files where present,
    so that runs from older LaREST versions pick up the current parsing.

    Parameters
    ----------
    dir_path : Path
        Root output directory for the molecule (e.g.
        ``output/Monomer/<slug>``).
    temperature : float
        Temperature in Kelvin used to recompute the corrected free energy when
        restoring a CREST entropy checkpoint (see :func:`apply_entropy_correction`).
    steps : dict[str, bool] | None
        The ``[steps]`` config section indicating which stages are enabled.  Keys
        are ``rdkit``, ``crest_confgen``, ``censo``, and ``crest_entropy``.  A
        missing key is treated as enabled.  When ``None`` (the default), every
        stage is treated as enabled.
    config : dict[str, Any] | None
        Full pipeline configuration, used to detect stages whose settings have
        changed.  When ``None``, no settings check is performed.
    smiles : str | None
        SMILES string of the molecule, stored with the stage settings.

    Returns
    -------
    results : dict[str, dict[str, float | None]]
        Mapping of section name to thermodynamic parameters populated from any
        checkpoint files that were found.  Parameters that have not yet been
        computed remain ``None``.
    stage : PipelineStage
        The earliest enabled pipeline stage whose results were not found.  The
        caller should (re-)run this stage and all subsequent ones.
    """
    if steps is None:
        steps = {}

    results: dict[str, dict[str, float | None]] = {
        section: dict.fromkeys(THERMODYNAMIC_PARAMS, None)
        for section in PIPELINE_SECTIONS
    }

    if steps.get("rdkit", True) and (
        _settings_changed(dir_path, PipelineStage.RDKIT, config, smiles)
        or not _load_stage(
            path=dir_path / "xtb" / "rdkit" / "results.csv",
            load_fn=lambda path: _parse_rdkit(results, path),
            label="rdkit",
        )
    ):
        return results, PipelineStage.RDKIT

    # without xTB re-ranking, CREST writes no results; its ensemble is the checkpoint
    crest_checkpoint = (
        dir_path / "xtb" / "crest" / "results.json"
        if steps.get("xtb", True)
        else dir_path / "crest_confgen" / "crest_conformers.xyz"
    )
    if steps.get("crest_confgen", True) and (
        _settings_changed(dir_path, PipelineStage.CREST_CONFGEN, config, smiles)
        or not _load_stage(
            path=crest_checkpoint,
            load_fn=lambda path: (
                _parse_crest(results, path, temperature)
                if steps.get("xtb", True)
                else None
            ),
            label="crest_confgen",
        )
    ):
        return results, PipelineStage.CREST_CONFGEN

    if steps.get("censo", True) and (
        _settings_changed(dir_path, PipelineStage.CENSO, config, smiles)
        or not _load_stage(
            path=dir_path / "censo" / "results.json",
            load_fn=lambda path: _parse_censo(results, path, temperature),
            label="censo",
        )
    ):
        return results, PipelineStage.CENSO

    if steps.get("crest_entropy", True) and (
        _settings_changed(dir_path, PipelineStage.CREST_ENTROPY, config, smiles)
        or not _load_stage(
            path=dir_path / "crest_entropy" / "results.json",
            load_fn=lambda path: results.update(
                combine_entropy(
                    results,
                    json.loads(path.read_text()),
                    temperature,
                    censo_enabled=steps.get("censo", True),
                ),
            ),
            label="crest_entropy",
        )
    ):
        return results, PipelineStage.CREST_ENTROPY

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

    Results from older LaREST versions lack the electronic energy ``E``; it is
    recovered from the best conformer's xTB output where available.

    Parameters
    ----------
    results : dict
        Shared results dict to update in-place with the ``"rdkit"`` section.
    path : Path
        Path to the xTB results CSV produced during the RDKit stage.
    """
    best = parse_best_rdkit_conformer(path)
    conformer_id = int(best.pop("conformer_id") or 0)
    energy = best.get("E")
    if energy is None or math.isnan(energy):
        xtb_output_file = (
            path.parent / f"conformer_{conformer_id}" / f"conformer_{conformer_id}.txt"
        )
        best["E"] = _parse_xtb_energy(xtb_output_file)
    results.update({"rdkit": _complete(best)})


def _parse_crest(results: dict, path: Path, temperature: float) -> None:
    """Load the CREST/xTB results into *results*, recovering ``E`` if missing."""
    crest = json.loads(path.read_text())
    if crest.get("E") is None:
        crest["E"] = _parse_xtb_energy(path.parent / "xtb.txt", temperature)
    results.update({"crest": _complete(crest)})


def _parse_censo(results: dict, path: Path, temperature: float) -> None:
    """Load CENSO results, re-deriving them from CENSO's raw output if present."""
    censo_dir = path.parent
    if get_last_censo_section(censo_dir) is not None:
        results.update(parse_censo_results(censo_dir, temperature))
    else:
        logger.warning(
            f"CENSO output files missing in {censo_dir}; using stored results.json",
        )
        results.update(
            {
                section: _complete(params)
                for section, params in json.loads(path.read_text()).items()
            },
        )


def _parse_xtb_energy(
    xtb_output_file: Path,
    temperature: float = 298.15,
) -> float | None:
    if not xtb_output_file.exists():
        return None
    return parse_xtb_output(xtb_output_file, temperature)["E"]


def _complete(params: dict[str, float | None]) -> dict[str, float | None]:
    """Return *params* with any missing thermodynamic parameters set to ``None``."""
    return {param: params.get(param) for param in THERMODYNAMIC_PARAMS}


def get_entropy_reference_section(
    results: Mapping[str, Mapping[str, float | None]],
) -> str:
    """Return the CENSO section the CREST entropy correction is applied to.

    This is the lowest free energy conformer section (``<stage>_min``) of the
    last CENSO sub-stage with an enthalpy and entropy, normally
    ``censo_refinement_min``, matching the reference structure used for the
    CREST entropy run (see :func:`larest.crest.get_best_structure`).
    """
    for section in reversed(CENSO_SECTIONS):
        params = results.get(f"{section}{MIN_SUFFIX}", {})
        if params.get("H") is not None and params.get("S") is not None:
            return f"{section}{MIN_SUFFIX}"
    return f"censo_refinement{MIN_SUFFIX}"


def apply_entropy_correction(
    refinement_results: dict[str, float | None],
    crest_entropy_results: dict[str, float | None],
    temperature: float,
) -> dict[str, dict[str, float | None]]:
    """Apply a CREST conformational entropy correction to CENSO refinement results.

    Following Grimme and co-workers, the correction is applied to the lowest
    free energy conformer of the CENSO refinement (``censo_refinement_min``),
    which is also the reference structure used for the CREST entropy run.  That
    stage yields the single-conformer
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
        Thermodynamic parameters of the lowest free energy conformer from the
        CENSO refinement stage, keyed by ``"E"``, ``"H"``, ``"S"``, and ``"G"``.
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


def combine_entropy(
    results: Mapping[str, Mapping[str, float | None]],
    crest_entropy_results: Mapping[str, float | None],
    temperature: float,
    *,
    censo_enabled: bool = True,
) -> dict[str, dict[str, float | None]]:
    """Combine CREST entropy results with the CENSO results of a molecule.

    If CENSO results with an enthalpy and entropy are available, the CREST
    conformational entropy correction is applied to them (see
    :func:`get_entropy_reference_section` and :func:`apply_entropy_correction`),
    giving a ``"censo_corrected"`` section.  Otherwise (CENSO disabled, or no
    mRRHO contribution evaluated) the raw entropies are returned as a
    ``"crest_entropy"`` section.

    Parameters
    ----------
    results : dict[str, dict[str, float | None]]
        Results of the molecule, including any CENSO sections.
    crest_entropy_results : dict[str, float | None]
        Entropy parameters from the CREST entropy run.
    temperature : float
        Temperature in Kelvin used to recompute the corrected free energy.
    censo_enabled : bool
        Whether the CENSO stage is enabled.

    Returns
    -------
    dict[str, dict[str, float | None]]
        Either ``{"censo_corrected": {...}}`` or ``{"crest_entropy": {...}}``.
    """
    if censo_enabled:
        reference = results.get(get_entropy_reference_section(results), {})
        if reference.get("H") is not None and reference.get("S") is not None:
            return apply_entropy_correction(
                dict(reference),
                dict(crest_entropy_results),
                temperature,
            )
        logger.warning(
            "No CENSO enthalpy and entropy to apply the CREST entropy correction "
            "to; storing the CREST entropies as crest_entropy",
        )
    return {"crest_entropy": dict(crest_entropy_results)}
