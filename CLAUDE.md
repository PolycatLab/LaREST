# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What is LaREST

LaREST (Lactone Ring-opening Energetics Sorting Tool) is a computational chemistry pipeline for calculating thermodynamic parameters (H, S, G) of lactone ring-opening polymerization reactions. It chains together RDKit, xTB, CREST, CENSO, and ORCA to generate and refine conformer ensembles, ultimately predicting ring-opening energetics.

## Installation

```bash
conda env create -f environment.yaml
conda activate larest
pip install .
```

CENSO and ORCA must be installed separately (see README). The `orcaversion` in `config/config.toml` must match the installed ORCA version.

## Running LaREST

```bash
# Local run
larest -o <output_dir> -c <config_dir>

# HPC (Imperial College PBS)
qsub pipeline.sh
```

The `-c` flag points to a directory containing `config.toml`. Settings in `config.toml` override `config/defaults.toml`.

## Linting

```bash
ruff check src/
ruff format src/
```

Ruff is configured in `pyproject.toml` with `select = ["ALL"]` and several categories ignored (see `[tool.ruff.lint].ignore`). Docstrings (`D`) and pathlib enforcement (`PTH`) are disabled.

## Architecture

### Pipeline stages (in order)

Each molecule passes through up to four stages, controlled by `[steps]` in config:

1. **rdkit** — Generate MMFF conformers, rank by xTB free energy (`G`)
2. **crest_confgen** — Conformer/rotamer ensemble via CREST v3 from best RDKit conformer, xTB re-ranking
3. **censo** — DFT refinement of CREST ensemble via CENSO (4 sub-stages: `0_PRESCREENING` → `1_SCREENING` → `2_OPTIMIZATION` → `3_REFINEMENT`), using ORCA as the QM backend
4. **crest_entropy** — Conformational entropy via CREST entropy mode; adds `censo_corrected` result by applying CREST entropy correction to `3_REFINEMENT` results

### Molecule class hierarchy (`src/larest/base.py`)

All molecule types inherit from `LarestMol` (abstract base):
- `Monomer` — the lactone monomer; owns `Initiator` and a list of `Polymer` instances
- `Polymer` — polymer chain at a given length; SMILES is built on-the-fly via `build_polymer()` in `chem.py`
- `Initiator` — the initiating alcohol (ROR reactions only)

Each molecule stores results as `dict[section, dict[param, float|None]]` where sections are `rdkit`, `crest`, `0_PRESCREENING`, ..., `3_REFINEMENT`, `censo_corrected` and params are `H`, `S`, `G`.

### Checkpointing (`src/larest/checkpoint.py`)

On startup, each `LarestMol` calls `restore_results()` which walks through result files in order and sets `_pipeline_stage` to the first missing stage. This lets interrupted runs resume from where they left off without re-running earlier stages.

### Configuration system

`get_config()` in `setup.py` deep-merges `defaults.toml` (all options) with the user's `config.toml` (overrides only). Config keys are passed directly as CLI flags to external tools via `parse_command_args()` in `parsers.py` — a `true` boolean value becomes `--flag`, a scalar becomes `--key value`, `false` is omitted.

### Reaction types

- **RER** (Ring Equilibrium Reaction) — polymer chain only, no initiator contribution; requires `lengths >= 2`
- **ROR** (Ring-Opening polymerization Reaction) — initiator is consumed; requires `lengths >= 1`

Final thermodynamics are computed in `Monomer.compile_results()`:
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

`tests/` contains benchmarking and validation scripts (`ring_size/`, `tropic/`, `tu-2023/`) with their own `pipeline.sh` and plotting scripts. These are research validation datasets, not a pytest test suite.
