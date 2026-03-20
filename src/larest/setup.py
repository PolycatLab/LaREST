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


def _apply_parallelisation(config: dict[str, Any], user_config: dict[str, Any]) -> dict[str, Any]:
    n_cores = config.get("parallelisation", {}).get("n_cores")
    if n_cores is None:
        return config

    for path, key in _PARALLELISATION_KEYS.items():
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
            section[key] = n_cores

    return config


def get_config(config_dir: Path) -> dict[str, Any]:
    """Load and merge the pipeline configuration.

    Deep-merges the user's ``config.toml`` on top of the built-in defaults
    from ``src/larest/defaults.toml``, then propagates
    ``[parallelisation] n_cores`` to every stage that the user has not
    explicitly configured.

    Parameters
    ----------
    config_dir : Path
        Directory containing ``config.toml``.

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

    config_file = config_dir / "config.toml"
    try:
        with open(config_file, "rb") as fstream:
            user_config = tomllib.load(fstream)
    except Exception:
        print(f"Failed to load config from {config_file}")
        raise

    config = _deep_merge(defaults, user_config)
    return _apply_parallelisation(config, user_config)


def get_logger(output_dir: Path, config: dict[str, Any]) -> None:
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
        log_config["handlers"]["file"]["filename"] = (
            output_dir / log_config["handlers"]["file"]["filename"]
        ).resolve()
        logging.config.dictConfig(log_config)
    except Exception:
        print(f"Failed to setup logging config from {config}")
        raise


def parse_command_args(sub_config: list[str], config: dict[str, Any]) -> list[str]:
    """Convert a nested config section into a flat list of CLI flags.

    Traverses *config* using the keys in *sub_config* to locate the target
    section, then converts each key-value pair to ``--key value`` tokens.
    Boolean values are handled specially: ``True`` emits ``--key`` with no
    value; ``False`` is omitted entirely.

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
        logger.warning(f"Failed to find sub-config {sub_config} in config, using no arguments")
        return []

    args = []
    for k, v in cfg.items():
        if v is False:
            continue
        args.append(f"--{k}")
        if v is not True:
            args.append(str(v))
    logger.debug(f"Returning {sub_config} args: {args}")
    return args

