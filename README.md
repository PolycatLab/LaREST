# LaREST: Lactone Ring-opening Energetics Sorting Tool

LaREST is a computational chemistry pipeline for predicting the thermodynamics of lactone ring-opening polymerization reactions. Given a set of monomer SMILES strings and polymer chain lengths, it computes Boltzmann-averaged enthalpies (H), entropies (S), and free energies (G) by chaining together four levels of theory:

1. **RDKit** — MMFF conformer generation, ranked by xTB free energy
2. **CREST** — conformer/rotamer ensemble exploration (iMTD-GC)
3. **CENSO** — DFT refinement of the ensemble via ORCA (four sub-stages of increasing accuracy)
4. **CREST entropy** — conformational entropy correction applied to the CENSO results

Two reaction types are supported:

- **RER** (Ring Equilibrium Reaction) — no initiator; polymer chain only
- **ROR** (Ring-Opening polymerization Reaction) — includes an initiating alcohol

Final reaction thermodynamics are computed as:

```
ΔX = (X_polymer - n·X_monomer - X_initiator) / n
```

where *n* is the polymer chain length.

## Installation

```bash
git clone https://github.com/Ryan-Reese/LaREST.git
cd LaREST
conda env create -f environment.yaml
conda activate larest
```

This installs `xTB`, `CREST`, and `LaREST` itself (along with all Python dependencies) into the `larest` conda environment.

**ORCA** must be installed separately and available on your `PATH`. See the [ORCA website](https://orcaforum.kofo.mpg.de) for download instructions.

> [!IMPORTANT]
> Update the `orcaversion` field in [config/config.toml](./config/config.toml) to match your installed ORCA version.

## Usage

### Running locally

```bash
larest -o <output_dir> -c <config_dir>
```

- `-o` / `--output` — directory where all results are written (default: `./output`)
- `-c` / `--config` — directory containing `config.toml` (default: `./config`)
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

### `[reaction]`

The most important section — defines what to compute.

```toml
[reaction]
type = "RER"                  # "RER" (no initiator) or "ROR" (with initiator)
monomers = ["O=C1OCC1"]      # list of monomer SMILES
lengths = [2, 3, 4]           # polymer chain lengths to evaluate
initiator = "C1=CC=C(C=C1)CO" # initiator SMILES (required for ROR only)
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

Options passed directly as CLI flags to the `xtb` binary. The `etemp` key is used internally to derive entropy via S = (H − G) / T and is not passed as a flag.

| Key | Description |
|---|---|
| `etemp` | Electronic temperature in K (required) |
| `gfn` | GFN-xTB method level (0, 1, or 2) |
| `alpb` | Implicit solvent (e.g. `"toluene"`) |
| `ohess` | Hessian level after optimisation (e.g. `"vtight"`) |
| `parallel` | Number of CPU threads |

### `[crest.confgen]` and `[crest.entropy]`

Options passed directly as CLI flags to `crest`. Separate sections control the conformer generation and entropy calculation runs. Key options:

| Key | Description |
|---|---|
| `T` | Number of CPU threads |
| `gfn2` / `gfnff` | Level of theory for CREST |
| `alpb` | Implicit solvent |
| `ewin` | Energy window for conformer ensemble (kcal/mol) |
| `rthr` | RMSD threshold for conformer deduplication (Å) |

### `[censo.*]`

Controls the four CENSO DFT sub-stages. Each sub-stage (`prescreening`, `screening`, `optimization`, `refinement`) has its own section with `func`, `basis`, `sm`, and `threshold` keys. The `[censo.general]` section sets global CENSO settings including `temperature` and `solvent`.

> [!IMPORTANT]
> `[censo.paths] orcaversion` must match your installed ORCA version exactly.

## Output

All results are written under the output directory with the following structure:

```
output/
  Monomer/<smiles_slug>/
    rdkit/
      conformers.sdf          # all MMFF conformers
      xtb/rdkit/results.csv   # xTB H, S, G for each conformer
    crest_confgen/
      crest_best.xyz          # lowest-energy CREST conformer
      crest_conformers.xyz    # full CREST ensemble
    censo/
      censo.txt               # raw CENSO output
      results.json            # parsed H, S, G for each sub-stage
      3_REFINEMENT.xyz        # refined conformer ensemble
    crest_entropy/
      crest.txt               # raw CREST entropy output
      results.json            # S_conf, S_rrho, S_total
    results.json              # merged results across all sections
    summary/
      rdkit.csv               # ΔH, ΔS, ΔG vs polymer length (RDKit level)
      crest.csv               # ΔH, ΔS, ΔG vs polymer length (CREST level)
      3_REFINEMENT.csv        # ΔH, ΔS, ΔG vs polymer length (DFT level)
      censo_corrected.csv     # ΔH, ΔS, ΔG with CREST entropy correction

  Polymer/<smiles_slug>_<length>/
    ...                       # same structure as Monomer

  Initiator/<smiles_slug>/
    ...                       # same structure as Monomer

  larest.log                  # full pipeline log
```

The `summary/` CSVs are the primary output, containing per-polymer-length reaction thermodynamics (ΔH, ΔS, ΔG in J/mol) at each level of theory. The `censo_corrected` section combines CENSO `3_REFINEMENT` enthalpies with the CREST conformational entropy.

## Tested versions

Dependency | Version
--- | ---
`xTB` | 6.7.1
`CREST` | 3.0.2
`CENSO` | 2.1.4
`ORCA` | 6.1.0

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
