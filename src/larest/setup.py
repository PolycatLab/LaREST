"""Configuration loading and logging setup for the LaREST pipeline.

:func:`get_config` deep-merges the user's ``config.toml`` on top of the
built-in defaults in ``defaults.toml``, then propagates the top-level
``[parallelisation] n_cores``, ``[thermo] temperature`` and ``[thermo] solvent``
to each stage unless overridden.

:func:`get_logger` applies the ``[logging]`` section of that dict via
:func:`logging.config.dictConfig`, substituting the resolved path for the
log file handler.

:func:`parse_command_args` converts a nested sub-config dict into a flat
list of CLI flags suitable for passing directly to ``subprocess.run``.
"""

from __future__ import annotations

import copy
import importlib.resources
import logging
import logging.config
import tomllib
from pathlib import Path
from typing import Any

from rdkit.Chem.rdmolfiles import MolFromSmiles
from rdkit.rdBase import BlockLogs

logger = logging.getLogger(__name__)

# Maps (config path, key) to the parallelisation key used by each stage.
# If the user does not explicitly set the key, it is filled from [parallelisation].n_cores.
_PARALLELISATION_KEYS: dict[tuple[str, ...], str] = {
    ("rdkit",): "n_cores",
    ("xtb",): "parallel",
    ("crest", "confgen"): "T",
    ("crest", "entropy"): "T",
    ("censo", "cli"): "maxcores",
}

# Maps (config path, key) to the temperature key used by each stage.
# If the user does not explicitly set the key, it is filled from [thermo].temperature.
# Note: [xtb].temperature is a code-only key (see parse_command_args); the other
# entries map onto the native temperature flags of CENSO (temperature) and CREST (temp).
_TEMPERATURE_KEYS: dict[tuple[str, ...], str] = {
    ("xtb",): "temperature",
    ("censo", "general"): "temperature",
    ("crest", "confgen"): "temp",
    ("crest", "entropy"): "temp",
}

# Maps config path to the solvent key used by each stage. If the user does not
# set a solvent key for a stage, it is filled from [thermo].solvent.
_SOLVENT_KEYS: dict[tuple[str, ...], str] = {
    ("xtb",): "alpb",
    ("crest", "confgen"): "alpb",
    ("crest", "entropy"): "alpb",
    ("censo", "general"): "solvent",
}
# Keys that, if set by the user in an xTB/CREST stage, mean they have chosen the
# solvation for that stage themselves.
_USER_SOLVATION_KEYS: frozenset[str] = frozenset({"alpb", "gbsa", "cpcmx"})
_CENSO_GAS_PHASE_KEYS: frozenset[str] = frozenset({"gas_phase", "gas-phase"})
# [thermo].solvent values that select gas-phase calculations.
GAS_PHASE_VALUES: frozenset[str] = frozenset({"", "gas", "none", "vacuum"})

# Keys that are used internally by LaREST and must never be emitted as CLI flags
# by parse_command_args (e.g. the code-only thermostatistical temperature, which
# xTB receives via a $thermo xcontrol block rather than a --temperature flag).
_NON_CLI_KEYS: frozenset[str] = frozenset({"temperature"})


def _load_defaults() -> dict[str, Any]:
    data = importlib.resources.files("larest").joinpath("defaults.toml").read_bytes()
    return tomllib.loads(data.decode("utf-8"))


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _propagate_default(
    config: dict[str, Any],
    user_config: dict[str, Any],
    keys: dict[tuple[str, ...], str],
    value: object,
) -> dict[str, Any]:
    """Fill a shared default *value* into each per-stage key in *keys*.

    For every ``(path, key)`` entry, the key is set to *value* in *config*
    unless the user explicitly set it in *user_config* (in which case their
    value is preserved).
    """
    for path, key in keys.items():
        # Walk user_config to check if the user explicitly set this key
        user_section: Any = user_config
        explicitly_set = False
        for p in path:
            if not isinstance(user_section, dict) or p not in user_section:
                break
            user_section = user_section[p]
        else:
            explicitly_set = isinstance(user_section, dict) and key in user_section

        if not explicitly_set:
            section = config
            for p in path:
                section = section.setdefault(p, {})
            section[key] = value

    return config


def _apply_parallelisation(
    config: dict[str, Any],
    user_config: dict[str, Any],
) -> dict[str, Any]:
    n_cores = config.get("parallelisation", {}).get("n_cores")
    if n_cores is None:
        return config
    return _propagate_default(config, user_config, _PARALLELISATION_KEYS, n_cores)


def _apply_temperature(
    config: dict[str, Any],
    user_config: dict[str, Any],
) -> dict[str, Any]:
    temperature = config.get("thermo", {}).get("temperature")
    if temperature is None:
        return config
    return _propagate_default(config, user_config, _TEMPERATURE_KEYS, temperature)


def _get_section(config: dict[str, Any], path: tuple[str, ...]) -> object:
    section: Any = config
    for p in path:
        if not isinstance(section, dict) or p not in section:
            return None
        section = section[p]
    return section


def _apply_solvent(
    config: dict[str, Any],
    user_config: dict[str, Any],
) -> dict[str, Any]:
    """Propagate ``[thermo].solvent`` to the solvation settings of every stage.

    A solvent name sets ``--alpb <solvent>`` for xTB and CREST and the CENSO
    ``solvent``; a gas-phase value (``"gas"``, ``"none"``, ``"vacuum"``, ``""``
    or ``false``) disables implicit solvation in xTB and CREST and sets CENSO's
    ``gas_phase = true``.  Stages for which the user explicitly set a
    solvation key are left untouched, and a warning is logged if the stages
    end up with inconsistent solvation.
    """
    solvent = config.get("thermo", {}).get("solvent")
    if solvent is None:
        return config
    gas_phase = solvent is False or str(solvent).strip().lower() in GAS_PHASE_VALUES

    for path, key in _SOLVENT_KEYS.items():
        user_section = _get_section(user_config, path)
        user_keys = set(user_section) if isinstance(user_section, dict) else set()
        section = config
        for p in path:
            section = section.setdefault(p, {})
        if path == ("censo", "general"):
            # the solvent and gas-phase switch are set independently, so that
            # e.g. a user-set gas_phase = false still gets the solvent
            # an explicit CENSO solvent means the user chose solvation for CENSO
            if not user_keys & (_CENSO_GAS_PHASE_KEYS | {key}):
                section["gas_phase"] = gas_phase
            if key not in user_keys and not gas_phase:
                section[key] = solvent
        elif not user_keys & _USER_SOLVATION_KEYS:
            section[key] = False if gas_phase else solvent

    solvation = {
        "xtb": config.get("xtb", {}).get("alpb") or None,
        "crest.confgen": config.get("crest", {}).get("confgen", {}).get("alpb") or None,
        "crest.entropy": config.get("crest", {}).get("entropy", {}).get("alpb") or None,
        "censo": None
        if config.get("censo", {}).get("general", {}).get("gas_phase")
        or config.get("censo", {}).get("general", {}).get("gas-phase")
        else config.get("censo", {}).get("general", {}).get("solvent"),
    }
    if len(set(solvation.values())) > 1:
        logger.warning(
            "Implicit solvation differs between stages "
            f"({', '.join(f'{k}: {v or "gas phase"}' for k, v in solvation.items())})",
        )
    return config


def get_config(config_file: Path) -> dict[str, Any]:
    """Load and merge the pipeline configuration.

    Deep-merges the user's ``config.toml`` on top of the built-in defaults
    from ``src/larest/defaults.toml``, then propagates
    ``[parallelisation] n_cores`` to every stage that the user has not
    explicitly configured.

    Parameters
    ----------
    config_file : Path
        Path to the ``config.toml`` file.

    Returns
    -------
    dict[str, Any]
        Merged configuration dict.

    Raises
    ------
    Exception
        Re-raises any exception raised during config loading after printing
        a human-readable error message.
    """
    defaults = _load_defaults()

    try:
        with open(config_file, "rb") as fstream:
            user_config = tomllib.load(fstream)
    except Exception:
        print(f"Failed to load config from {config_file}")  # noqa: T201
        raise

    config = _deep_merge(defaults, user_config)
    config = _apply_parallelisation(config, user_config)
    config = _apply_temperature(config, user_config)
    config = _apply_solvent(config, user_config)
    return normalise_monomers(config)


_SMART_QUOTES = "\u201c\u201d\u2018\u2019"


def _check_smiles(smiles: object, label: str) -> str | None:
    """Return an error message if *smiles* is not a valid SMILES string."""
    if not isinstance(smiles, str) or not smiles.strip():
        return f"{label} must be a non-empty SMILES string (got {smiles!r})"
    if any(quote in smiles for quote in _SMART_QUOTES):
        return (
            f'{label} {smiles!r} contains curly quotes; use straight quotes (") '
            "in config.toml (curly quotes are often added by word processors)"
        )
    with BlockLogs():
        if MolFromSmiles(smiles) is None:
            return f"{label} {smiles!r} is not a valid SMILES string"
    return None


def validate_reaction_config(config: dict[str, Any]) -> list[str]:
    """Check the ``[reaction]`` settings before any calculation is started.

    Parameters
    ----------
    config : dict[str, Any]
        Full pipeline configuration dict (after :func:`get_config`).

    Returns
    -------
    list[str]
        Human-readable problems; empty if the settings are valid.
    """
    reaction = config.get("reaction", {})
    errors: list[str] = []

    reaction_type = reaction.get("type")
    if reaction_type not in ("RER", "ROR"):
        errors.append(f'[reaction].type must be "RER" or "ROR" (got {reaction_type!r})')

    lengths = reaction.get("lengths", [])
    min_length = 2 if reaction_type == "RER" else 1
    if not lengths:
        errors.append("[reaction].lengths must list at least one polymer length")
    errors.extend(
        f"[reaction].lengths must be integers >= {min_length} for "
        f"{reaction_type} (got {length!r})"
        for length in lengths
        if not isinstance(length, int)
        or isinstance(length, bool)
        or length < min_length
    )

    monomers = reaction.get("monomers", [])
    if not monomers:
        errors.append("[reaction].monomers must list at least one monomer SMILES")
    for monomer in monomers:
        smiles = monomer["smiles"] if isinstance(monomer, dict) else monomer
        if (error := _check_smiles(smiles, "monomer")) is not None:
            errors.append(error)

    if reaction_type == "ROR":
        initiator = reaction.get("initiator", "")
        if not initiator:
            errors.append(
                '[reaction].initiator is required for ROR (a SMILES string or "auto")',
            )
        elif not (isinstance(initiator, str) and initiator.strip().lower() == "auto"):
            error = _check_smiles(initiator, "[reaction].initiator")
            if error is not None:
                errors.append(f'{error}; use a SMILES string or "auto"')
    return errors


def normalise_monomers(config: dict[str, Any]) -> dict[str, Any]:
    """Normalise ``[reaction].monomers`` entries to dicts.

    Each monomer may be given either as a SMILES string or as a table with a
    ``smiles`` key and an optional ``group`` key selecting the functional group
    to open (overriding the priority ranking), e.g.
    ``{ smiles = "O=C1CCC(=O)N1", group = "Lm" }``.

    Parameters
    ----------
    config : dict[str, Any]
        Full pipeline configuration dict.

    Returns
    -------
    dict[str, Any]
        The configuration, with ``config["reaction"]["monomers"]`` as a list of
        ``{"smiles": str, "group": str | None}`` dicts.

    Raises
    ------
    ValueError
        If a monomer entry is neither a string nor a table with ``smiles``.
    """
    monomers: list[dict[str, Any]] = []
    for monomer in config.get("reaction", {}).get("monomers", []):
        if isinstance(monomer, str):
            monomers.append({"smiles": monomer, "group": None})
        elif isinstance(monomer, dict) and "smiles" in monomer:
            monomers.append(
                {"smiles": monomer["smiles"], "group": monomer.get("group")},
            )
        else:
            raise ValueError(f"Invalid monomer entry in config: {monomer}")
    config.setdefault("reaction", {})["monomers"] = monomers
    return config


def get_logger(
    output_dir: Path,
    config: dict[str, Any],
    *,
    verbose: bool = False,
) -> None:
    """Configure the Python logging system from the ``[logging]`` config section.

    Deep-copies the logging sub-config before patching the file handler's
    ``filename`` with an absolute path inside *output_dir*, then applies it
    with :func:`logging.config.dictConfig`.

    Parameters
    ----------
    output_dir : Path
        Root output directory; the log file is written here using the filename
        specified in ``config["logging"]["handlers"]["file"]["filename"]``.
    config : dict[str, Any]
        Full pipeline configuration dict.  Must contain a ``[logging]`` section
        compatible with :func:`logging.config.dictConfig`.

    Raises
    ------
    Exception
        Re-raises any exception raised during logging setup after printing a
        human-readable error message.
    """
    try:
        log_config: dict[str, Any] = copy.deepcopy(config["logging"])
        log_config["handlers"]["file"]["filename"] = str(
            (output_dir / log_config["handlers"]["file"]["filename"]).resolve(),
        )
        if verbose:
            log_config["handlers"]["stream"]["level"] = "DEBUG"
        logging.config.dictConfig(log_config)
    except Exception:
        print(f"Failed to setup logging config from {log_config}")  # noqa: T201
        raise


def parse_command_args(sub_config: list[str], config: dict[str, Any]) -> list[str]:
    """Convert a nested config section into a flat list of CLI flags.

    Traverses *config* using the keys in *sub_config* to locate the target
    section, then converts each key-value pair to ``--key value`` tokens.
    Boolean values are handled specially: ``True`` emits ``--key`` with no
    value; ``False`` is omitted entirely.  Code-only keys in
    :data:`_NON_CLI_KEYS` (e.g. ``temperature``) are skipped entirely.

    Parameters
    ----------
    sub_config : list[str]
        Sequence of keys used to navigate to the desired section, e.g.
        ``["xtb"]`` or ``["crest", "confgen"]``.
    config : dict[str, Any]
        Full pipeline configuration dict.

    Returns
    -------
    list[str]
        Flat list of CLI argument tokens ready to be appended to a
        ``subprocess.run`` command list.  Returns an empty list if the
        sub-config path cannot be found.
    """
    cfg = config
    try:
        for key in sub_config:
            cfg = cfg[key]
    except KeyError:
        logger.warning(
            f"Failed to find sub-config {sub_config} in config, using no arguments",
        )
        return []

    args = []
    for k, v in cfg.items():
        if k in _NON_CLI_KEYS:
            continue
        if v is False:
            continue
        args.append(f"--{k}")
        if v is not True:
            args.append(str(v))
    logger.debug(f"Returning {sub_config} args: {args}")
    return args
