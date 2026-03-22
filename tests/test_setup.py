"""Tests for larest.setup configuration utilities."""

from __future__ import annotations

import pytest

from larest.setup import (
    _apply_parallelisation,
    _deep_merge,
    get_config,
    parse_command_args,
)


class TestDeepMerge:
    def test_simple_override(self):
        base = {"a": 1, "b": 2}
        override = {"b": 99}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": 99}

    def test_nested_merge(self):
        base = {"rdkit": {"n_conformers": 50, "random_seed": 42}}
        override = {"rdkit": {"n_conformers": 10}}
        result = _deep_merge(base, override)
        assert result["rdkit"]["n_conformers"] == 10
        assert result["rdkit"]["random_seed"] == 42

    def test_does_not_mutate_base(self):
        base = {"a": {"x": 1}}
        override = {"a": {"x": 2}}
        _deep_merge(base, override)
        assert base["a"]["x"] == 1

    def test_new_key_added(self):
        base = {"a": 1}
        override = {"b": 2}
        result = _deep_merge(base, override)
        assert result == {"a": 1, "b": 2}

    def test_override_scalar_with_dict(self):
        base = {"a": 1}
        override = {"a": {"nested": 2}}
        result = _deep_merge(base, override)
        assert result["a"] == {"nested": 2}


class TestApplyParallelisation:
    def test_propagates_n_cores(self):
        config = {
            "parallelisation": {"n_cores": 4},
            "rdkit": {},
            "xtb": {},
            "crest": {"confgen": {}, "entropy": {}},
            "censo": {"cli": {}},
        }
        result = _apply_parallelisation(config, user_config={})
        assert result["rdkit"]["n_cores"] == 4
        assert result["xtb"]["parallel"] == 4
        assert result["crest"]["confgen"]["T"] == 4
        assert result["crest"]["entropy"]["T"] == 4
        assert result["censo"]["cli"]["maxcores"] == 4

    def test_user_override_not_clobbered(self):
        config = {
            "parallelisation": {"n_cores": 4},
            "xtb": {"parallel": 2},
        }
        user_config = {"xtb": {"parallel": 2}}
        result = _apply_parallelisation(config, user_config=user_config)
        assert result["xtb"]["parallel"] == 2

    def test_no_n_cores_key_is_noop(self):
        config = {"rdkit": {}}
        result = _apply_parallelisation(config, user_config={})
        assert result == {"rdkit": {}}


class TestParseCommandArgs:
    def test_scalar_value(self):
        config = {"xtb": {"gfn": 2, "etemp": 298.15}}
        args = parse_command_args(["xtb"], config)
        assert "--gfn" in args
        assert "2" in args
        assert "--etemp" in args
        assert "298.15" in args

    def test_bool_true_emits_flag_only(self):
        config = {"xtb": {"ohess": True}}
        args = parse_command_args(["xtb"], config)
        assert "--ohess" in args
        # no value token after the flag
        idx = args.index("--ohess")
        assert idx == len(args) - 1 or not args[idx + 1].startswith("True")

    def test_bool_false_omitted(self):
        config = {"xtb": {"ohess": False}}
        args = parse_command_args(["xtb"], config)
        assert "--ohess" not in args

    def test_nested_path(self):
        config = {"crest": {"confgen": {"T": 4, "gfn2": True}}}
        args = parse_command_args(["crest", "confgen"], config)
        assert "--T" in args
        assert "4" in args
        assert "--gfn2" in args

    def test_missing_key_returns_empty(self):
        config = {}
        args = parse_command_args(["nonexistent"], config)
        assert args == []

    def test_empty_section_returns_empty(self):
        config = {"xtb": {}}
        args = parse_command_args(["xtb"], config)
        assert args == []


class TestGetConfig:
    def test_loads_and_merges_config(self, tmp_path):
        config_file = tmp_path / "config.toml"
        config_file.write_text(
            '[reaction]\nmonomers = ["C1CC(=O)O1"]\nlengths = [2]\n',
        )
        config = get_config(config_file)
        assert config["reaction"]["monomers"] == ["C1CC(=O)O1"]
        # defaults should be present
        assert "rdkit" in config
        assert config["rdkit"]["n_conformers"] == 50

    def test_user_overrides_default(self, tmp_path):
        config_file = tmp_path / "config.toml"
        config_file.write_text("[rdkit]\nn_conformers = 5\n")
        config = get_config(config_file)
        assert config["rdkit"]["n_conformers"] == 5
        # other rdkit defaults still present
        assert config["rdkit"]["random_seed"] == 42

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            get_config(tmp_path / "nonexistent.toml")
