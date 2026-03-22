# Contributing to LaREST

## Development setup

```bash
conda env create -f env.yaml
conda activate larest
pip install -e ".[dev]"
pre-commit install
```

This installs the package in editable mode along with ruff, ty, pre-commit, and pytest.

## Linting and type checking

Pre-commit runs ruff (lint + format) and ty (type checking) automatically on each commit. To run against all files manually:

```bash
pre-commit run --all-files
```

Or run tools directly:

```bash
ruff check src/
ruff format src/
ty check src/
```

## Testing

```bash
# Unit tests only (no external tools required)
pytest tests/

# Include integration tests (requires xtb and crest on PATH)
pytest tests/ --integration
```

## Submitting changes

1. Fork the repository and create a branch from `main`.
2. Make your changes, ensuring all tests pass and pre-commit checks are clean.
3. Open a pull request against `main` with a clear description of what was changed and why.

## Versioning

Versions are derived from git tags of the form `vX.Y.Z`. To cut a release, create an annotated tag:

```bash
git tag -a v1.0.0 -m "Release v1.0.0"
```
