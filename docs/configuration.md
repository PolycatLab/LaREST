# Configuration

Pipeline behaviour is controlled by a `config.toml` file. User settings are **deep-merged** on top of the built-in defaults in `src/larest/defaults.toml`, so you only need to set values that differ from defaults.

A minimal starting point is `config/example.config`. The full reference is `config/reference.toml`.

## `[reaction]`

The most important section — defines what to compute.

```toml
[reaction]
type = "RER"                         # "RER" (no initiator) or "ROR" (with initiator)
monomers = ["O=C1OCC1"]             # list of monomer SMILES
lengths = [2, 3, 4]                  # polymer chain lengths to evaluate
initiator = "C1=CC=C(C=C1)CO"       # initiator SMILES (required for ROR only)
```

- **RER** (Ring Equilibrium Reaction) requires `lengths >= 2`.
- **ROR** (Ring-Opening polymerization Reaction) requires `lengths >= 1` and an `initiator` SMILES.

## `[parallelisation]`

Set the core count once; it propagates to every stage automatically.

```toml
[parallelisation]
n_cores = 16
```

| Stage | Config key set | CLI flag |
|---|---|---|
| RDKit | `[rdkit].n_cores` | (internal) |
| xTB | `[xtb].parallel` | `--parallel N` |
| CREST confgen | `[crest.confgen].T` | `--T N` |
| CREST entropy | `[crest.entropy].T` | `--T N` |
| CENSO | `[censo.cli].maxcores` | `--maxcores N` |

Override a specific stage by setting its key directly:

```toml
[parallelisation]
n_cores = 16

[censo.cli]
maxcores = 8    # CENSO uses fewer cores (e.g. limited by ORCA memory)
```

## `[steps]`

Toggle individual pipeline stages on or off.

```toml
[steps]
rdkit = true          # RDKit conformer generation + xTB ranking
crest_confgen = true  # CREST conformer ensemble
censo = true          # CENSO DFT refinement
crest_entropy = true  # CREST conformational entropy correction
xtb = true            # xTB re-ranking after rdkit and crest_confgen
```

## `[rdkit]`

Controls MMFF conformer generation.

| Key | Default | Description |
|---|---|---|
| `n_conformers` | `50` | Number of conformers to generate |
| `n_cores` | `1` | CPU cores for parallel embedding |
| `mmff` | `"MMFF94"` | MMFF variant: `"MMFF94"` or `"MMFF94s"` |
| `random_seed` | `42` | Seed for reproducibility |

## `[xtb]`

All keys (except `etemp`) are passed directly as CLI flags to the `xtb` binary.

| Key | Default | Description |
|---|---|---|
| `etemp` | `298.15` | Electronic temperature (K); used to derive S = (H − G) / T |
| `gfn` | `2` | GFN-xTB method level (0, 1, 2) |
| `alpb` | `"toluene"` | ALPB implicit solvent |
| `ohess` | `"vtight"` | Hessian level after optimisation |
| `parallel` | `1` | Number of CPU threads |

## `[crest.confgen]` and `[crest.entropy]`

All keys are passed directly as CLI flags to the `crest` binary. Separate sections control conformer generation and entropy calculation.

| Key | Default | Description |
|---|---|---|
| `T` | `1` | Number of CPU threads |
| `gfn2` / `gfnff` | | Level of theory |
| `alpb` | `"toluene"` | ALPB implicit solvent |
| `ewin` | `6.0` | Energy window for ensemble (kcal/mol) |
| `rthr` | `0.125` | RMSD threshold for deduplication (Å) |

## `[censo.*]`

Controls the four CENSO DFT sub-stages. Each sub-stage (`prescreening`, `screening`, `optimization`, `refinement`) has its own section.

| Key | Description |
|---|---|
| `func` | DFT functional |
| `basis` | Basis set |
| `sm` | Solvation model (`"smd"`, `"cpcm"`, `"cosmo"`) |
| `threshold` | Energy window for ensemble pruning (kcal/mol) |

The `[censo.general]` section sets global CENSO settings:

| Key | Default | Description |
|---|---|---|
| `temperature` | `298.15` | Temperature (K) |
| `solvent` | `"toluene"` | Solvent for implicit solvation |

## Boolean flags

For `[xtb]` and `[crest.*]` sections, boolean values are translated to CLI flags:

- `true` → `--flag`
- `false` → flag is omitted entirely
- scalar → `--key value`
