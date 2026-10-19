# Configuration

Pipeline behaviour is controlled by a `config.toml` file. User settings are **deep-merged** on top of the built-in defaults in `src/larest/defaults.toml`, so you only need to set values that differ from defaults.

A minimal starting point is `config/example.toml`. The full reference is `config/reference.toml`.

## `[reaction]`

The most important section — defines what to compute.

```toml
[reaction]
type = "ROR"                         # "RER" (no initiator) or "ROR" (with initiator)
monomers = [                         # monomer SMILES, or tables with a group override
    "O=C1OCC1",
    { smiles = "O=C1OCCC(=O)N1", group = "oA" },
]
lengths = [1, 2, 3]                  # polymer chain lengths to evaluate
initiator = "auto"                   # initiator SMILES or "auto" (required for ROR only)
```

- **RER** (Ring Expansion Reaction) builds cyclic oligomers and requires `lengths >= 2`.
- **ROR** (Ring-Opening Reaction) requires `lengths >= 1` and an `initiator`.

### Monomer check

Before any calculation, each monomer is checked: it must be cyclic and contain
at least one supported ring-opening functional group. LaREST reports the number
and sizes of the rings, the functional groups found, the bond that will be
broken and the size of the ring it belongs to, and any stereocentres without an
assigned configuration (conformer generation would otherwise mix
stereoisomers, so specify stereochemistry in the SMILES). Unsuitable monomers
are skipped. The results appear in `larest_summary.txt`.

The bond broken is between the acyl carbon (C=O or C=S) and the ring
heteroatom. When several breakable bonds are present, the highest priority
group is opened:

| Priority | Group | Key |
|---|---|---|
| 1 | lactone | `L` |
| 2 | thionolactone | `tnL` |
| 3 | thiolactone | `tL` |
| 4 | dithiolactone | `dtL` |
| 5 | lactam | `Lm` |
| 6 | thiolactam | `tLm` |
| 7 | carbonates and derivatives (not yet benchmarked) | `CC`, `CtC`, `CdtC`, `CtnC`, `CX`, `CtX`, `oA` |

Give a monomer as a table with a `group` key to open a different group.

### Supported initiator types

| Type | Example SMILES | Bond broken | Product (γ-butyrolactone, n = 1) |
|---|---|---|---|
| Water | `O` | O–H | `HO–C(=O)–(CH2)3–OH` (hydrolysis) |
| Alcohol (R–OH) | `CO` (methanol) | O–H | `MeO–C(=O)–(CH2)3–OH` (methanolysis) |
| Ester / thioester / amide | `CC(=O)OC` (methyl acetate) | acyl C–X | `MeO–C(=O)–(CH2)3–O–C(=O)Me` (isodesmic) |
| `"auto"` | — | — | isodesmic initiator chosen per monomer |

With `initiator = "auto"`, each monomer gets the acyclic analogue of its
ring-opening group, so the reaction is isodesmic (the types of bonds broken and
formed are conserved), as recommended by Houk, Hall and co-workers: methyl
acetate for lactones, O-methyl thioacetate for thionolactones, S-methyl
thioacetate for thiolactones, methyl dithioacetate for dithiolactones,
N-methylacetamide for lactams, dimethyl carbonate for carbonates, and so on.
Monomers needing different initiators can therefore be run together; each
initiator is computed once and shared between monomers.

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

`solvent` sets the implicit solvation of every stage in the same way:

```toml
[thermo]
solvent = "toluene"     # or "gas" for gas-phase calculations
```

| Stage | Config key set |
|---|---|
| xTB | `[xtb].alpb` |
| CREST confgen / entropy | `[crest.confgen].alpb`, `[crest.entropy].alpb` |
| CENSO | `[censo.general].solvent` and `gas_phase` |

With `solvent = "gas"` (or `"none"`, `"vacuum"`, `false`), ALPB is switched
off for xTB and CREST and `gas_phase = true` is set for CENSO. A stage for which
you set a solvation key yourself (e.g. `[xtb].alpb`, `[xtb].gbsa`,
`[censo.general].gas_phase`) keeps your setting; LaREST logs a warning if the
stages end up with different solvation. Note that xTB/CREST (ALPB) and CENSO
(SMD by default) may name some solvents differently; override per stage if
needed.

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
| `alpb` | from `[thermo].solvent` | ALPB implicit solvent |
| `ohess` | `"vtight"` | Hessian level after optimisation |
| `parallel` | `1` | Number of CPU threads |

## `[crest.confgen]` and `[crest.entropy]`

All keys are passed directly as CLI flags to the `crest` binary. Separate sections control conformer generation and entropy calculation.

| Key | Default | Description |
|---|---|---|
| `T` | `1` | Number of CPU threads |
| `temp` | from `[thermo]` | Temperature (K) for Boltzmann weighting (`--temp`). Inherited from `[thermo].temperature` unless set here |
| `gfn2` / `gfnff` | | Level of theory |
| `alpb` | from `[thermo].solvent` | ALPB implicit solvent |
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
| `solvent` | from `[thermo].solvent` | Solvent for implicit solvation |
| `sm_rrho` | `"alpb"` | Solvation model for the mRRHO contribution |
| `evaluate_rrho` | `true` | Include the GFN2-xTB mRRHO (thermal + entropy) contribution |
| `gas_phase` | from `[thermo].solvent` | Run all calculations in the gas phase, overriding **all** solvation settings |

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
