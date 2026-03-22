# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is LaREST

LaREST (Lactone Ring-opening Energetics Sorting Tool) is a computational chemistry pipeline for calculating thermodynamic parameters (H, S, G) of lactone ring-opening polymerization reactions. It chains together RDKit, xTB, CREST, CENSO, and ORCA to generate and refine conformer ensembles, ultimately predicting ring-opening energetics.

## Installation

```bash
conda env create -f env.yaml
conda activate larest
pip install .
```

CENSO and ORCA must be installed separately (see README).

## Running LaREST

```bash
# Local run
larest <config.toml> -o <output_dir>

# HPC (Imperial College PBS)
qsub pipeline.sh
```

The `-c` flag points to a `config.toml` file. See `config/reference.toml` for documentation of all available options.

## Development setup

Install the dev dependencies and set up pre-commit hooks:

```bash
pip install -e ".[dev]"
pre-commit install
```

## Linting

Pre-commit runs ruff and ty automatically on each commit. To run manually:

```bash
pre-commit run --all-files
```

Or run tools directly:

```bash
ruff check src/
ruff format src/
ty check src/
```

Ruff is configured in `pyproject.toml` with `select = ["ALL"]` and several categories ignored (see `[tool.ruff.lint].ignore`). Docstrings (`D`) and pathlib enforcement (`PTH`) are disabled.

## Testing

```bash
# Unit tests only (no external tools required)
pytest tests/

# Include integration tests (requires xtb and crest on PATH)
pytest tests/ --integration
```

Tests live in `tests/` and are structured as follows:

- `conftest.py` — shared fixtures (minimal pipeline config, sample file paths) and the `--integration` CLI flag
- `data/` — sample output files used as fixtures (xTB output, CREST entropy output, rdkit results CSV, CENSO output, multi-conformer XYZ)
- `test_data.py` — dataclass unit tests
- `test_output.py` — filesystem utility tests (`slugify`, `create_dir`, `remove_dir`)
- `test_setup.py` — config loading, deep-merge, parallelisation propagation, `parse_command_args`
- `test_chem.py` — SMILES parsing, ring-size detection, RER/ROR polymer construction
- `test_xtb.py` — `parse_xtb_output`, `run_xtb` (mocked subprocess + integration)
- `test_crest.py` — `parse_crest_entropy_output`, `run_crest_confgen/entropy` (mocked + integration)
- `test_censo.py` — all pure-parsing functions in `censo.py` (no ORCA/CENSO binary required)
- `test_checkpoint.py` — `PipelineStage`, `restore_results`, `apply_entropy_correction`
- `test_rdkit_stage.py` — `parse_best_rdkit_conformer`, `run_rdkit` (mocked + integration)
- `test_main.py` — `compile_results` delta calculations, `run_pipeline` orchestration

Integration tests are marked `@pytest.mark.integration` and skipped by default. Pass `--integration` to enable them. CENSO/ORCA tests are not included as they require a full DFT installation.

## Architecture

### Pipeline stages (in order)

Each molecule passes through up to four stages, controlled by `[steps]` in config:

1. **rdkit** — Generate MMFF conformers, rank by xTB free energy (`G`)
2. **crest_confgen** — Conformer/rotamer ensemble via CREST v3 from best RDKit conformer, xTB re-ranking
3. **censo** — DFT refinement of CREST ensemble via CENSO (4 sub-stages: `censo_prescreening` → `censo_screening` → `censo_optimization` → `censo_refinement`), using ORCA as the QM backend
4. **crest_entropy** — Conformational entropy via CREST entropy mode; adds `censo_corrected` result by applying CREST entropy correction to `censo_refinement` results

### Molecule dataclasses (`src/larest/data.py`)

Molecules are simple dataclasses with no shared base class:
- `Monomer` — the lactone monomer (`smiles: str`)
- `Polymer` — polymer chain at a given length (`smiles`, `monomer_smiles`, `length`); SMILES is built on-the-fly via `build_polymer()` in `chem.py`
- `Initiator` — the initiating alcohol for ROR reactions (`smiles: str`)
- `MolResults` — holds results for one molecule: `smiles` and `sections: dict[str, dict[str, float|None]]` where sections are `rdkit`, `crest`, `censo_prescreening`, ..., `censo_refinement`, `censo_corrected` and params are `H`, `S`, `G`

### Checkpointing (`src/larest/checkpoint.py`)

At the start of each molecule's run, `MolPipeline.run()` calls `restore_results()` which walks through result files in order and returns the first missing `PipelineStage`. This lets interrupted runs resume from where they left off without re-running earlier stages.

### Configuration system

`get_config()` in `setup.py` loads `config.toml` directly (no defaults merging). Config keys are passed directly as CLI flags to external tools via `parse_command_args()` in `parsers.py` — a `true` boolean value becomes `--flag`, a scalar becomes `--key value`, `false` is omitted.

### Reaction types

- **RER** (Ring Equilibrium Reaction) — polymer chain only, no initiator contribution; requires `lengths >= 2`
- **ROR** (Ring-Opening polymerization Reaction) — initiator is consumed; requires `lengths >= 1`

Final thermodynamics are computed in `compile_results()` in `main.py`:
`delta_param = (polymer_param - n * monomer_param - initiator_param) / n`

### Output structure

```
output/
  Monomer/<slugified_smiles>/
    rdkit/            # conformers.sdf, xtb/rdkit/results.csv
    crest_confgen/    # crest_best.xyz, xtb/crest/results.json
    censo/            # censo.txt, results.json, 3_REFINEMENT.xyz
    crest_entropy/    # crest.txt, results.json
    results.json      # final merged results for all sections
    summary/          # <section>.csv with delta_H, delta_S, delta_G per polymer length
  Polymer/<slugified_smiles>_<length>/
    ...               # same structure as Monomer
  Initiator/<slugified_smiles>/
    ...
```

### Tests directory

`tests/` contains the pytest test suite (see **Testing** above) plus sample fixture data in `tests/data/`.
