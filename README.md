# LaREST: Lactone Ring-opening Energetics Sorting Tool

LaREST is a computational chemistry pipeline for predicting the thermodynamics of lactone ring-opening polymerization reactions. Given a set of monomer SMILES strings and polymer chain lengths, it computes Boltzmann-averaged enthalpies (H), entropies (S), and free energies (G) by chaining together four levels of theory:

1. **RDKit** — MMFF conformer generation, ranked by xTB free energy
2. **CREST** — conformer/rotamer ensemble exploration (iMTD-GC)
3. **CENSO** — DFT refinement of the ensemble via ORCA (four sub-stages of increasing accuracy)
4. **CREST entropy** — conformational entropy correction applied to the CENSO results

Two reaction types are supported:

- **RER** (Ring Expansion Reaction) — no initiator; polymer chain only
- **ROR** (Ring-Opening polymerization Reaction) — includes an initiating alcohol



## Installation

```bash
git clone https://github.com/PolycatLab/LaREST.git
cd LaREST
conda env create -f env.yaml
conda activate larest
pip install .
```

This installs Python ≥ 3.12, `xTB` and `CREST` (via conda) and `LaREST` with `CENSO` (via pip) into the `larest` conda environment. Python 3.12 or later is required (CENSO 3 needs it).

**ORCA** must be installed separately and available on your `PATH` (before `/usr/bin`, which on some Linux systems contains an unrelated `orca` screen reader). See the [ORCA website](https://orcaforum.kofo.mpg.de) for download instructions. The shared-library ORCA builds need the matching OpenMPI version for parallel runs (e.g. `conda install -c conda-forge "openmpi>=4.1,<5"` for the `openmpi418` builds).

## Usage

### Running locally

```bash
larest <config.toml> -o <output_dir>
```

- `config` — path to `config.toml` (required)
- `-o` / `--output` — directory where all results are written (default: `./output`)
- `-v` / `--verbose` — increase console log verbosity

### Running on HPC (Imperial College PBS)

A PBS job script is included for Imperial College's HPC service. Edit the resource directives and the `MODIFY ME` block at the top of `pipeline.sh`, then submit from the LaREST root directory:

```bash
qsub pipeline.sh
```

Key variables to set in `pipeline.sh`:

| Variable | Description |
|---|---|
| `CONDA_DIR` | Path to the directory containing your conda binary |
| `CONDA_ENV` | Conda environment name (default: `larest`) |
| `N_CORES` | Number of cores matching the PBS `ncpus` directive |
| `CONFIG_FILE` | Path to the config file (default `config.toml`) |
| `OUTPUT_DIR` | Output directory (default `output`) |

Relative paths are relative to the directory the job is submitted from. `CONFIG_FILE` and `OUTPUT_DIR` can also be set per job without editing the script:

```bash
qsub -v CONFIG_FILE=configs/gbl.toml,OUTPUT_DIR=runs/gbl pipeline.sh
```

### Checkpointing

LaREST automatically checkpoints after each pipeline stage. If a run is interrupted, re-running the same command will resume from the first incomplete stage — completed stages are not re-executed.

## Configuration

Pipeline behaviour is controlled by `config.toml` in the config directory. User settings are **deep-merged** on top of the built-in defaults in [`src/larest/defaults.toml`](./src/larest/defaults.toml), so you only need to set the values that differ from defaults. A minimal starting point is provided in [`config/example.config`](./config/example.config). See [`config/reference.toml`](./config/reference.toml) for documentation of every available option.

### Parallelisation

Set the core count once at the top level and it propagates to every stage:

```toml
[parallelisation]
n_cores = 16
```

This fills in the stage-specific parallelisation keys automatically:

| Stage | Config key | CLI flag |
|---|---|---|
| RDKit | `[rdkit].n_cores` | (internal, `numThreads`) |
| xTB | `[xtb].parallel` | `--parallel N` |
| CREST confgen | `[crest.confgen].T` | `--T N` |
| CREST entropy | `[crest.entropy].T` | `--T N` |
| CENSO | `[censo.cli].maxcores` | `--maxcores N` |

To override for a specific stage, set that key directly in your `config.toml`:

```toml
[parallelisation]
n_cores = 16      # default for all stages

[censo.cli]
maxcores = 8      # CENSO uses fewer cores (e.g. limited by memory per ORCA job)
```

> [!NOTE]
> CENSO's `maxcores` sets the total CPU budget across all concurrent ORCA jobs.
> A common pattern is to set `n_cores` to the total available cores and lower
> `maxcores` if ORCA jobs are memory-constrained.

### `[thermo]`

Sets the temperature once; it propagates to every stage that evaluates
thermostatistical quantities (H, S, G):

```toml
[thermo]
temperature = 298.15    # K; sets both the T each tool evaluates G at and the
                        # T used to derive S = (H - G) / T
```

| Stage | Config key | How it is applied |
|---|---|---|
| xTB | `[xtb].temperature` | injected via a `$thermo` xcontrol block |
| CENSO | `[censo.general].temperature` | written to `.censo2rc` |
| CREST confgen | `[crest.confgen].temp` | `--temp N` |
| CREST entropy | `[crest.entropy].temp` | `--temp N` |

Override a single stage by setting its key directly (`temperature` for xTB/CENSO,
`temp` for CREST). This is distinct from `[xtb].etemp`, the electronic
(Fermi-smearing) temperature, which does **not** set the thermostatistical T.

### `[reaction]`

The most important section — defines what to compute.

```toml
[reaction]
type = "ROR"                  # "RER" (no initiator) or "ROR" (with initiator)
monomers = ["O=C1OCC1"]       # list of monomer SMILES
lengths = [1, 2, 3]           # polymer chain lengths to evaluate
initiator = "auto"            # initiator SMILES, or "auto" (required for ROR only)
```

Each monomer is first checked (cyclic, supported functional group, ring size,
unassigned stereochemistry) and skipped if it cannot be ring-opened. When a
monomer has several breakable bonds, the highest priority group is opened
(lactone > thionolactone > thiolactone > dithiolactone > lactam > thiolactam >
carbonates).

Supported ROR initiators are water (`"O"`), alcohols (e.g. `"CO"`), and
esters, thioesters and amides (e.g. methyl acetate, `"CC(=O)OC"`).
`initiator = "auto"` chooses the isodesmic initiator for each monomer (e.g.
methyl acetate for lactones, N-methylacetamide for lactams). See
[docs/configuration.md](docs/configuration.md) for details.

### `[thermo]`

```toml
[thermo]
temperature = 298.15  # K, used by every stage
solvent = "toluene"   # implicit solvent for every stage, or "gas"
```

### `[steps]`

Toggle individual pipeline stages on or off.

```toml
[steps]
rdkit = true          # RDKit conformer generation + xTB ranking
crest_confgen = true  # CREST conformer ensemble
censo = true          # CENSO DFT refinement
crest_entropy = true  # CREST conformational entropy correction
xtb = true            # xTB re-ranking after rdkit and crest_confgen stages
```

### `[rdkit]`

Controls MMFF conformer generation. Key options:

| Key | Description |
|---|---|
| `n_conformers` | Number of conformers to generate |
| `n_cores` | CPU cores for parallel embedding and optimisation |
| `mmff` | MMFF variant: `"MMFF94"` or `"MMFF94s"` |
| `random_seed` | Seed for reproducibility |

### `[xtb]`

Options passed directly as CLI flags to the `xtb` binary, except `temperature`. The `temperature` key (inherited from `[thermo]`) is injected via a `$thermo` xcontrol block so xTB evaluates its free energy at that T, and is used to derive entropy via S = (H − G) / T. `etemp` is the electronic (Fermi-smearing) temperature only and does not affect the thermostatistics.

| Key | Description |
|---|---|
| `temperature` | Thermostatistical temperature in K for H/S/G (inherited from `[thermo].temperature`) |
| `etemp` | Electronic (Fermi-smearing) temperature in K — a numerical SCF parameter, not the thermostatistical T |
| `gfn` | GFN-xTB method level (0, 1, or 2) |
| `alpb` | Implicit solvent (default from `[thermo].solvent`) |
| `ohess` | Hessian level after optimisation (e.g. `"vtight"`) |
| `parallel` | Number of CPU threads |

### `[crest.confgen]` and `[crest.entropy]`

Options passed directly as CLI flags to `crest`. Separate sections control the conformer generation and entropy calculation runs. Key options:

| Key | Description |
|---|---|
| `T` | Number of CPU threads |
| `gfn2` / `gfnff` | Level of theory for CREST |
| `alpb` | Implicit solvent (default from `[thermo].solvent`) |
| `ewin` | Energy window for conformer ensemble (kcal/mol) |
| `rthr` | RMSD threshold for conformer deduplication (Å) |

### `[censo.*]`

Controls the four CENSO DFT sub-stages. Each sub-stage (`prescreening`, `screening`, `optimization`, `refinement`) has its own section with `func`, `basis`, `sm`, and `threshold` keys. The `[censo.general]` section sets global CENSO settings including `temperature`, `solvent`, and `gas_phase` (all inherited from `[thermo]` unless set here).

> [!NOTE]
> Keys under `[censo.*]` are written verbatim into CENSO's `.censo2rc`, and CENSO
> **silently ignores** unrecognised keys. CENSO field names use underscores — use
> `gas_phase`, not `gas-phase`. (LaREST normalises hyphens to underscores when
> writing the file, but the underscore form is preferred.)

## Output

All results are written under the output directory. The most useful files are
collected in `larest_summary.txt` (a human-readable overview of the run) and
`results/<monomer>/` (final structures, summary CSVs, energy matrices and a
machine-readable `results.json`). See [docs/output.md](docs/output.md) for the
full structure and definitions.

## Development

Install the dev dependencies (ruff, ty, pre-commit, pytest) and register the git hooks:

```bash
pip install -e ".[dev]"
pre-commit install
```

Pre-commit will run ruff (lint + format) and ty (type checking) automatically on each commit. To run against all files manually:

```bash
pre-commit run --all-files
```

## Testing

The test suite lives in `tests/` and uses `pytest`. Tests are split into fast unit tests (no external tools required) and integration tests that exercise the `xtb` and `crest` binaries.

**Unit tests only** (fast, no external tools needed):

```bash
pytest tests/
```

**Including integration tests** (requires `xtb` and `crest` on `PATH`):

```bash
pytest tests/ --integration
```

## Citations

For `xTB`:
- C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme
  *WIREs Comput. Mol. Sci.*, **2020**, 11, e01493.
  DOI: [10.1002/wcms.1493](https://doi.org/10.1002/wcms.1493)
- S. Grimme, C. Bannwarth, P. Shushkov,
  *J. Chem. Theory Comput.*, **2017**, 13, 1989-2009.
  DOI: [10.1021/acs.jctc.7b00118](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- C. Bannwarth, S. Ehlert and S. Grimme.,
  *J. Chem. Theory Comput.*, **2019**, 15, 1652-1671.
  DOI: [10.1021/acs.jctc.8b01176](https://dx.doi.org/10.1021/acs.jctc.8b01176)
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme,
  *ChemRxiv*, **2019**, preprint.
  DOI: [10.26434/chemrxiv.8326202.v1](https://dx.doi.org/10.26434/chemrxiv.8326202.v1)
- S. Ehlert, M. Stahn, S. Spicher, S. Grimme,
  *J. Chem. Theory Comput.*, **2021**, 17, 4250-4261
  DOI: [10.1021/acs.jctc.1c00471](https://doi.org/10.1021/acs.jctc.1c00471)

For `CREST`:
 - P. Pracht, S. Grimme, C. Bannwarth, F. Bohle, S. Ehlert, G. Feldmann, J. Gorges, M. Müller, T. Neudecker, C. Plett, S. Spicher, P. Steinbach, P. Wesołowski, F. Zeller,
   *J. Chem. Phys.*, **2024**, *160*, 114110.
   DOI: [10.1063/5.0197592](https://doi.org/10.1063/5.0197592)
