"""Shared pytest fixtures and configuration for the LaREST test suite."""

from __future__ import annotations

from pathlib import Path

import pytest

DATA_DIR = Path(__file__).parent / "data"


def pytest_addoption(parser):
    parser.addoption(
        "--integration",
        action="store_true",
        default=False,
        help="Run integration tests that invoke xtb/crest binaries.",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: marks tests that require xtb/crest to be installed",
    )


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--integration"):
        skip = pytest.mark.skip(reason="Pass --integration to run xtb/crest tests")
        for item in items:
            if "integration" in item.keywords:
                item.add_marker(skip)


# ---------------------------------------------------------------------------
# Minimal pipeline config fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def minimal_config():
    """Minimal pipeline config for unit tests (no external tools needed)."""
    return {
        "steps": {
            "rdkit": True,
            "crest_confgen": True,
            "censo": False,
            "crest_entropy": False,
            "xtb": True,
        },
        "parallelisation": {"n_cores": 1},
        "reaction": {
            "type": "RER",
            "lengths": [2],
            "initiator": "CCO",
            "monomers": ["C1CC(=O)O1"],
        },
        "rdkit": {
            "n_conformers": 1,
            "random_seed": 42,
            "conformer_box_size": 2.0,
            "mmff": "MMFF94",
            "mmff_iters": 10,
            "align_iters": 5,
            "precision": 6,
            "n_cores": 1,
        },
        "xtb": {
            "parallel": 1,
            "gfn": 2,
            "ohess": "vtight",
            "alpb": "toluene",
            "etemp": 298.15,
        },
        "crest": {
            "confgen": {
                "T": 1,
                "gfn2": True,
                "alpb": "toluene",
                "optlev": "vtight",
                "ewin": 6.0,
            },
            "entropy": {
                "T": 1,
                "gfnff": True,
                "alpb": "toluene",
                "optlev": "vtight",
                "ewin": 6.0,
            },
        },
        "censo": {
            "cli": {"maxcores": 1},
            "general": {"temperature": 298.15, "solvent": "toluene"},
            "prescreening": {"run": True},
            "screening": {"run": True},
            "optimization": {"run": True},
            "refinement": {"run": True},
        },
        "logging": {
            "version": 1,
            "disable_existing_loggers": False,
            "formatters": {"default": {"format": "%(levelname)s:%(message)s"}},
            "handlers": {
                "file": {
                    "class": "logging.FileHandler",
                    "level": "DEBUG",
                    "encoding": "utf-8",
                    "filename": "larest.log",
                    "formatter": "default",
                    "mode": "a",
                },
                "stream": {
                    "class": "logging.StreamHandler",
                    "level": "WARNING",
                    "formatter": "default",
                },
            },
            "root": {"handlers": ["file", "stream"], "level": "DEBUG"},
        },
    }


@pytest.fixture
def ror_config(minimal_config):
    """Minimal ROR config variant."""
    cfg = {**minimal_config}
    cfg["reaction"] = {**minimal_config["reaction"], "type": "ROR", "lengths": [1]}
    return cfg


# ---------------------------------------------------------------------------
# Path fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def data_dir():
    return DATA_DIR


@pytest.fixture
def xtb_output_file():
    return DATA_DIR / "xtb_output.txt"


@pytest.fixture
def crest_entropy_output_file():
    return DATA_DIR / "crest_entropy_output.txt"


@pytest.fixture
def rdkit_results_file():
    return DATA_DIR / "rdkit_results.csv"


@pytest.fixture
def censo_output_file():
    return DATA_DIR / "censo_output.txt"


@pytest.fixture
def censo_conformers_xyz_file():
    return DATA_DIR / "censo_conformers.xyz"
