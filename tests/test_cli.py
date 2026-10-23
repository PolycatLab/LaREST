"""Tests for the larest command-line interface (argument handling only)."""

from __future__ import annotations

import sys
from unittest.mock import patch

import pytest

from larest.cli import entry_point

_CONFIG = '[reaction]\ntype = "RER"\nmonomers = ["O=C1CCCO1"]\nlengths = [2]\n'


def _run(monkeypatch, tmp_path, *args):
    """Run the CLI in *tmp_path* with the pipeline mocked; return the config used."""
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(sys, "argv", ["larest", *args, "-o", str(tmp_path / "out")])
    with patch("larest.cli.main") as mock_main, patch("larest.cli.get_config") as get:
        get.return_value = {
            "reaction": {
                "type": "RER",
                "lengths": [2],
                "monomers": [{"smiles": "O=C1CCCO1", "group": None}],
            },
            "logging": {"version": 1, "handlers": {"file": {"filename": "l.log"}}},
        }
        with patch("larest.cli.get_logger"):
            entry_point()
    mock_main.assert_called_once()
    return get.call_args.args[0]


def test_default_config_file(monkeypatch, tmp_path):
    (tmp_path / "config.toml").write_text(_CONFIG)
    assert _run(monkeypatch, tmp_path) == tmp_path.resolve() / "config.toml"


def test_positional_config(monkeypatch, tmp_path):
    (tmp_path / "my.toml").write_text(_CONFIG)
    assert _run(monkeypatch, tmp_path, "my.toml") == tmp_path.resolve() / "my.toml"


@pytest.mark.parametrize("flag", ["-c", "--config"])
def test_config_option(monkeypatch, tmp_path, flag):
    (tmp_path / "my.toml").write_text(_CONFIG)
    used = _run(monkeypatch, tmp_path, flag, "my.toml")
    assert used == tmp_path.resolve() / "my.toml"


def test_same_config_twice_allowed(monkeypatch, tmp_path):
    (tmp_path / "my.toml").write_text(_CONFIG)
    used = _run(monkeypatch, tmp_path, "my.toml", "-c", "my.toml")
    assert used == tmp_path.resolve() / "my.toml"


def test_conflicting_configs_rejected(monkeypatch, tmp_path, capsys):
    (tmp_path / "a.toml").write_text(_CONFIG)
    (tmp_path / "b.toml").write_text(_CONFIG)
    with pytest.raises(SystemExit) as exc:
        _run(monkeypatch, tmp_path, "a.toml", "-c", "b.toml")
    assert exc.value.code == 2
    assert "two different config files" in capsys.readouterr().err


def test_missing_config_reported(monkeypatch, tmp_path, capsys):
    with pytest.raises(SystemExit) as exc:
        _run(monkeypatch, tmp_path)
    assert exc.value.code == 2
    assert "config file not found" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()
