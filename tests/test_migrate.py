"""Tests for larest.migrate (legacy output directory migration)."""

from __future__ import annotations

import pytest
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolToMolBlock
from rdkit.Chem.rdmolops import AddHs

from larest.chem import build_polymer, get_isodesmic_initiator
from larest.data import Initiator, Monomer, Polymer
from larest.main import get_mol_dir
from larest.migrate import apply_migration, entry_point, plan_migration
from larest.output import slugify

GBL = "C1CCC(=O)O1"  # gamma-butyrolactone
BL = "C1CC(=O)O1"  # beta-propiolactone
# beta-propiolactone and the unsaturated 2-oxetenone share a legacy slug
ALKENE = "C1=CC(=O)O1"


def _config(monomers, reaction_type="RER", lengths=(2, 3), initiator=None):
    reaction = {"type": reaction_type, "lengths": list(lengths), "monomers": monomers}
    if initiator is not None:
        reaction["initiator"] = initiator
    return {"reaction": reaction}


def _make_legacy(output_dir, sub, name, smiles):
    """Create a legacy molecule directory holding the RDKit conformers of *smiles*."""
    path = output_dir / sub / name
    (path / "rdkit").mkdir(parents=True)
    (path / "results.json").write_text("{}")
    (path / "rdkit" / "conformers.sdf").write_text(
        MolToMolBlock(AddHs(MolFromSmiles(smiles))) + "$$$$\n",
    )
    return path


def test_rer_plan(tmp_path):
    old_mon = _make_legacy(tmp_path, "Monomer", slugify(GBL), GBL)
    old_p2 = _make_legacy(
        tmp_path,
        "Polymer",
        f"{slugify(GBL)}_2",
        build_polymer(GBL, 2, "RER"),
    )
    plan = plan_migration(_config([GBL]), tmp_path)

    p2 = build_polymer(GBL, 2, "RER")
    expected = {
        old_mon: get_mol_dir(Monomer(GBL), tmp_path),
        old_p2: get_mol_dir(Polymer(p2, GBL, 2), tmp_path),
    }
    assert {m.old: m.new for m in plan.moves} == expected
    # missing length-3 polymer directory is ignored silently
    assert plan.warnings == []


def test_ror_plan_with_initiator(tmp_path):
    old_init = _make_legacy(tmp_path, "Initiator", slugify("CCO"), "CCO")
    old_p1 = _make_legacy(
        tmp_path,
        "Polymer",
        f"{slugify(GBL)}_1",
        build_polymer(GBL, 1, "ROR", "CCO"),
    )
    plan = plan_migration(
        _config([GBL], reaction_type="ROR", lengths=[1], initiator="CCO"),
        tmp_path,
    )
    p1 = build_polymer(GBL, 1, "ROR", initiator_smiles="CCO")
    moves = {m.old: m.new for m in plan.moves}
    assert moves[old_init] == get_mol_dir(Initiator("CCO"), tmp_path)
    assert moves[old_p1] == get_mol_dir(Polymer(p1, GBL, 1, "CCO"), tmp_path)


def test_ror_auto_initiator_shared(tmp_path):
    # two lactones share the same isodesmic initiator (methyl acetate)
    init = get_isodesmic_initiator(GBL)
    assert init == get_isodesmic_initiator(BL)
    old_init = _make_legacy(tmp_path, "Initiator", slugify(init), init)
    plan = plan_migration(
        _config([GBL, BL], reaction_type="ROR", lengths=[1], initiator="auto"),
        tmp_path,
    )
    assert [(m.old, m.new) for m in plan.moves] == [
        (old_init, get_mol_dir(Initiator(init), tmp_path)),
    ]
    assert plan.warnings == []


def test_slug_collision_not_moved(tmp_path):
    assert slugify(BL) == slugify(ALKENE)
    old_mon = _make_legacy(tmp_path, "Monomer", slugify(BL), BL)
    old_p2 = _make_legacy(
        tmp_path,
        "Polymer",
        f"{slugify(BL)}_2",
        build_polymer(BL, 2, "RER"),
    )
    old_gbl = _make_legacy(tmp_path, "Monomer", slugify(GBL), GBL)
    plan = plan_migration(_config([BL, ALKENE, GBL], lengths=[2]), tmp_path)

    assert [m.old for m in plan.moves] == [old_gbl]
    assert len(plan.warnings) == 2
    assert any(str(old_mon) in w and "ambiguous" in w for w in plan.warnings)
    assert any(str(old_p2) in w and "ambiguous" in w for w in plan.warnings)


def test_duplicate_monomer_not_collision(tmp_path):
    old_mon = _make_legacy(tmp_path, "Monomer", slugify(GBL), GBL)
    plan = plan_migration(_config([GBL, GBL], lengths=[2]), tmp_path)
    assert [m.old for m in plan.moves] == [old_mon]
    assert plan.warnings == []


def test_existing_target_not_overwritten(tmp_path):
    _make_legacy(tmp_path, "Monomer", slugify(GBL), GBL)
    get_mol_dir(Monomer(GBL), tmp_path).mkdir(parents=True)
    plan = plan_migration(_config([GBL], lengths=[2]), tmp_path)
    assert plan.moves == []
    assert len(plan.warnings) == 1
    assert "already exists" in plan.warnings[0]


def test_invalid_monomer_skipped(tmp_path):
    (tmp_path / "Monomer").mkdir()
    plan = plan_migration(_config(["CCO"], lengths=[2]), tmp_path)
    assert plan.moves == []
    assert len(plan.warnings) == 1
    assert "Skipping monomer CCO" in plan.warnings[0]


def test_not_an_output_dir(tmp_path):
    plan = plan_migration(_config([GBL], lengths=[2]), tmp_path)
    assert plan.moves == []
    assert len(plan.warnings) == 1
    assert "is it a LaREST output directory" in plan.warnings[0]


def test_unknown_reaction_type(tmp_path):
    with pytest.raises(ValueError, match="reaction type"):
        plan_migration(_config([GBL], reaction_type="XYZ"), tmp_path)


def test_apply_migration(tmp_path):
    old = _make_legacy(tmp_path, "Monomer", slugify(GBL), GBL)
    plan = plan_migration(_config([GBL], lengths=[2]), tmp_path)
    apply_migration(plan)
    new = get_mol_dir(Monomer(GBL), tmp_path)
    assert not old.exists()
    assert (new / "results.json").read_text() == "{}"


def test_apply_migration_refuses_overwrite(tmp_path):
    _make_legacy(tmp_path, "Monomer", slugify(GBL), GBL)
    plan = plan_migration(_config([GBL], lengths=[2]), tmp_path)
    plan.moves[0].new.mkdir(parents=True)
    with pytest.raises(FileExistsError):
        apply_migration(plan)


def _write_config(tmp_path):
    config_file = tmp_path / "config.toml"
    config_file.write_text(
        f'[reaction]\ntype = "RER"\nlengths = [2]\nmonomers = ["{GBL}"]\n',
    )
    return config_file


def test_entry_point_dry_run(tmp_path, capsys):
    out = tmp_path / "out"
    old = _make_legacy(out, "Monomer", slugify(GBL), GBL)
    entry_point([str(_write_config(tmp_path)), str(out)])
    assert old.exists()
    assert "Would move" in capsys.readouterr().out


def test_entry_point_apply(tmp_path, capsys):
    out = tmp_path / "out"
    old = _make_legacy(out, "Monomer", slugify(GBL), GBL)
    entry_point([str(_write_config(tmp_path)), str(out), "--apply"])
    assert not old.exists()
    assert get_mol_dir(Monomer(GBL), out.resolve()).is_dir()
    assert "Moved 1 directories" in capsys.readouterr().out


def test_mismatched_polymer_not_moved(tmp_path):
    """A legacy polymer directory holding a different molecule is not moved."""
    _make_legacy(tmp_path, "Polymer", f"{slugify(GBL)}_2", build_polymer(GBL, 3, "RER"))
    plan = plan_migration(_config([GBL]), tmp_path)
    assert plan.moves == []
    assert any("current version builds" in w for w in plan.warnings)


def test_unverifiable_directory_not_moved(tmp_path):
    path = tmp_path / "Monomer" / slugify(GBL)
    path.mkdir(parents=True)
    plan = plan_migration(_config([GBL]), tmp_path)
    assert plan.moves == []
    assert any("cannot verify" in w for w in plan.warnings)
