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

### Supported initiator types

| Type | Example SMILES | Bond broken | Chain ends |
|---|---|---|---|
| Alcohol (R–OH) | `CCO` (ethanol) | O–H | `HO–[chain]–C(=O)–OR` |
| Formate ester (HC(=O)–OR) | `COC=O` (methyl formate) | C(formyl)–O | `HC(=O)–O–[chain]–C(=O)–OCH3` |

Methyl formate is a non-symmetric initiator: it breaks its acyl C–O bond rather than an O–H bond, so **both** chain termini become ester groups (a formate ester at the ring-O end, a methyl ester at the acyl end). Any formate ester (`HC(=O)–OR`) is supported — the alkyl group `R` caps the acyl end.

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

## `[thermo]`

Sets the temperature once; it propagates to every stage that evaluates
thermostatistical quantities (H, S, G) automatically.

```toml
[thermo]
temperature = 298.15    # K; the temperature at which H, S and G are evaluated
```

| Stage | Config key set | How it is applied |
|---|---|---|
| xTB | `[xtb].temperature` | injected via a `$thermo` xcontrol block (see `[xtb]`) |
| CENSO | `[censo.general].temperature` | written to `.censo2rc` `[general]` |
| CREST confgen | `[crest.confgen].temp` | `--temp N` |
| CREST entropy | `[crest.entropy].temp` | `--temp N` |

`temperature` is the single source of truth for every H/S/G calculation: it sets
the temperature at which each tool evaluates the free energy **and** the `T` used
to derive `S = (H − G) / T`. Keeping it in one place ensures these stay
consistent across stages.

Override a specific stage by setting its key directly (note the per-stage key
names differ — `temperature` for xTB/CENSO, `temp` for CREST):

```toml
[thermo]
temperature = 298.15

[crest.entropy]
temp = 310.0    # run the entropy search at a different temperature
```

> **Note.** `[thermo].temperature` is distinct from `[xtb].etemp`, the electronic
> (Fermi-smearing) temperature, which is a numerical SCF parameter and does *not*
> set the thermostatistical temperature.

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

All keys except `temperature` are passed directly as CLI flags to the `xtb`
binary. `temperature` is **not** a CLI flag: it is injected via a `$thermo`
xcontrol block (`xtb --input …`) so that xTB evaluates its free energy at that
temperature, and it is used to derive `S = (H − G) / T`.

| Key | Default | Description |
|---|---|---|
| `temperature` | from `[thermo]` | Thermostatistical temperature (K) for H/S/G; injected via `$thermo`. Inherited from `[thermo].temperature` unless set here |
| `etemp` | `298.15` | Electronic (Fermi-smearing) temperature (K) — a numerical SCF parameter, **not** the thermostatistical temperature |
| `gfn` | `2` | GFN-xTB method level (0, 1, 2) |
| `alpb` | `"toluene"` | ALPB implicit solvent |
| `ohess` | `"vtight"` | Hessian level after optimisation |
| `parallel` | `1` | Number of CPU threads |

## `[crest.confgen]` and `[crest.entropy]`

All keys are passed directly as CLI flags to the `crest` binary. Separate sections control conformer generation and entropy calculation.

| Key | Default | Description |
|---|---|---|
| `T` | `1` | Number of CPU threads |
| `temp` | from `[thermo]` | Temperature (K) for Boltzmann weighting (`--temp`). Inherited from `[thermo].temperature` unless set here |
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
| `temperature` | from `[thermo]` | Temperature (K) for mRRHO/thermostatistics; used to derive `S = (H − G) / T`. Inherited from `[thermo].temperature` unless set here |
| `solvent` | `"toluene"` | Solvent for implicit solvation |
| `sm_rrho` | `"alpb"` | Solvation model for the mRRHO contribution |
| `evaluate_rrho` | `true` | Include the GFN2-xTB mRRHO (thermal + entropy) contribution |
| `gas_phase` | `false` | Run all calculations in the gas phase, overriding **all** solvation settings |

> **Note — key names must match CENSO's fields exactly.** Keys in `[censo.*]`
> are written verbatim into CENSO's `.censo2rc`, and CENSO **silently ignores**
> keys it does not recognise. CENSO field names use underscores, never hyphens.
> LaREST normalises hyphens to underscores when writing the file (so `gas-phase`
> is corrected to `gas_phase`), but prefer the underscore form. In particular,
> use `gas_phase` (not `gas-phase`) — with the hyphen, gas-phase would be
> silently dropped and solvation would remain active.

> **Note — where entropy comes from.** In CENSO the entropy (mRRHO) is computed
> by GFN2-xTB (`evaluate_rrho`), not by ORCA; ORCA only supplies the electronic
> single-point energy. The prescreening sub-stage has no mRRHO, so the `S` it
> reports is not a true thermodynamic entropy — treat `censo_screening` and later
> as the meaningful entropy stages.

## Boolean flags

For `[xtb]` and `[crest.*]` sections, boolean values are translated to CLI flags:

- `true` → `--flag`
- `false` → flag is omitted entirely
- scalar → `--key value`
