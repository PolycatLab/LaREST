"""Configuration loading and logging setup for the LaREST pipeline.

:func:`get_config` deep-merges the user's ``config.toml`` on top of the
built-in defaults in ``defaults.toml``, then propagates the top-level
``[parallelisation] n_cores`` to each stage unless overridden.

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
    return _apply_temperature(config, user_config)


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
