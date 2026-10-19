"""Tests for larest.output filesystem utilities."""

from larest.output import create_dir, mol_dir_name, remove_dir, slugify


class TestSlugify:
    def test_simple_smiles(self):
        assert slugify("CCO") == "CCO"

    def test_parentheses_become_hyphens(self):
        result = slugify("CC(=O)O")
        assert "(" not in result
        assert ")" not in result
        assert "-" in result

    def test_ring_smiles(self):
        result = slugify("C1CC(=O)O1")
        assert "(" not in result
        assert ")" not in result

    def test_at_sign_preserved(self):
        result = slugify("[C@@H](O)(N)C")
        assert "@" in result

    def test_consecutive_hyphens_collapsed(self):
        result = slugify("C(C)(C)C")
        assert "--" not in result

    def test_no_leading_trailing_hyphens(self):
        result = slugify("(CCO)")
        assert not result.startswith("-")
        assert not result.endswith("-")

    def test_special_chars_removed(self):
        result = slugify("C/C=C/C")
        assert "/" not in result

    def test_square_brackets(self):
        # [] brackets are word chars via \w in regex — [ and ] are not special
        # but they are not in [\w\s@-] so they get stripped
        result = slugify("[NH3+]")
        assert "[" not in result
        assert "]" not in result


class TestCreateDir:
    def test_creates_directory(self, tmp_path):
        new_dir = tmp_path / "subdir"
        assert not new_dir.exists()
        create_dir(new_dir)
        assert new_dir.exists()
        assert new_dir.is_dir()

    def test_creates_nested_directories(self, tmp_path):
        nested = tmp_path / "a" / "b" / "c"
        create_dir(nested)
        assert nested.exists()

    def test_idempotent(self, tmp_path):
        d = tmp_path / "existing"
        create_dir(d)
        create_dir(d)  # should not raise
        assert d.exists()


class TestRemoveDir:
    def test_removes_directory(self, tmp_path):
        d = tmp_path / "to_remove"
        d.mkdir()
        assert d.exists()
        remove_dir(d)
        assert not d.exists()

    def test_removes_non_empty_directory(self, tmp_path):
        d = tmp_path / "non_empty"
        d.mkdir()
        (d / "file.txt").write_text("content")
        remove_dir(d)
        assert not d.exists()

    def test_nonexistent_dir_does_not_raise(self, tmp_path):
        missing = tmp_path / "missing"
        remove_dir(missing)  # should log warning, not raise


class TestMolDirName:
    def test_distinct_for_molecules_with_same_slug(self):
        # slugify drops "=", so an alkene and its alkane share a slug
        alkene = "O=C1OCC2C1C1C=CC2O1"
        alkane = "O=C1OCC2C1C1CCC2O1"
        assert slugify(alkene) == slugify(alkane)
        assert mol_dir_name(alkene) != mol_dir_name(alkane)

    def test_deterministic(self):
        assert mol_dir_name("C1CC(=O)O1") == mol_dir_name("C1CC(=O)O1")

    def test_suffix_and_slug_source(self):
        name = mol_dir_name("O=C1CCOC(=O)CCO1", suffix="2", slug_smiles="O=C1CCO1")
        slug, length, digest = name.split("_")
        assert slug == slugify("O=C1CCO1")
        assert length == "2"
        assert len(digest) == 8

    def test_long_smiles_truncated(self):
        assert len(mol_dir_name("C" * 500)) < 100
