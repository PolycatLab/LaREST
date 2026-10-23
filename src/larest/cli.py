"""Command-line interface entry point for LaREST.

Parses the config file (positional or ``-c``/``--config``, default
``config.toml``), ``-o``/``--output``, and ``-v``/``--verbose`` arguments,
initialises the logger and configuration, then delegates to
:func:`larest.main.main`.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from larest import LAREST_HEADER
from larest.main import main
from larest.setup import get_config, get_logger, validate_reaction_config

logger = logging.getLogger(__name__)

DEFAULT_CONFIG_FILE = "config.toml"


def entry_point() -> None:
    """Parse CLI arguments and launch the LaREST pipeline.

    This function is registered as the ``larest`` console script entry point
    in ``pyproject.toml``.  It exits with code 1 if the config or logger
    cannot be initialised.
    """
    print(LAREST_HEADER)  # noqa: T201
    parser = argparse.ArgumentParser(description="LaREST")
    parser.add_argument(
        "config",
        type=str,
        nargs="?",
        default=None,
        help=f"Path to the config file (default: {DEFAULT_CONFIG_FILE})",
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        dest="config_option",
        metavar="CONFIG",
        default=None,
        help="Path to the config file (alternative to the positional argument)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="./output",
        help="Output directory",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="increase output verbosity",
        action="store_true",
    )
    args = parser.parse_args()

    if args.config and args.config_option and args.config != args.config_option:
        parser.error(
            f"two different config files given ({args.config} and "
            f"{args.config_option}); use one",
        )
    output_dir = Path(args.output).resolve()
    config_file = Path(
        args.config_option or args.config or DEFAULT_CONFIG_FILE,
    ).resolve()
    if not config_file.is_file():
        parser.error(f"config file not found: {config_file}")

    try:
        config = get_config(config_file)
    except Exception as err:
        raise SystemExit(1) from err

    errors = validate_reaction_config(config)
    if errors:
        print(f"Invalid settings in {config_file}:")  # noqa: T201
        for error in errors:
            print(f"  - {error}")  # noqa: T201
        raise SystemExit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        get_logger(output_dir=output_dir, config=config, verbose=args.verbose)
    except Exception as err:
        raise SystemExit(1) from err

    logger.info("LaREST Initialised")
    for config_key, config_value in config.items():
        logger.debug(f"{config_key} config:\n{config_value}")

    main(output_dir=output_dir, config=config)
