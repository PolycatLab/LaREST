# Output

All results are written under the directory specified by `-o` (default: `./output`).

## Directory structure

```
output/
  larest_summary.txt              # human-readable overview of the run (see below)
  larest.log                      # full pipeline log
  results/<monomer_dir>/          # key results for each monomer, in one place
    monomer.xyz                   # lowest free energy structure (final stage)
    initiator.xyz                 # (ROR only)
    polymer_n<n>.xyz              # one per chain length
    summary/                      # copy of the monomer's summary/ directory
    results.json                  # machine-readable results (see below)

  Monomer/<monomer_dir>/
    rdkit/
      conformers.sdf              # all MMFF conformers
    xtb/rdkit/results.csv         # xTB E, H, S, G for each conformer
    crest_confgen/
      crest_best.xyz              # lowest-energy CREST conformer
      crest_conformers.xyz        # full CREST ensemble
    censo/
      censo.txt                   # raw CENSO output
      <n>_<PART>.json             # CENSO per-conformer energies for each sub-stage
      results.json                # parsed E, H, S, G for each sub-stage
      3_REFINEMENT.xyz            # refined conformer ensemble
      censo_best.xyz              # lowest free energy refined conformer
    crest_entropy/
      crest.txt                   # raw CREST entropy output
      results.json                # S_conf, S_rrho, S_total
    <stage>/larest_settings.json  # settings each stage was run with
    results.json                  # merged results across all sections
    summary/
      <section>.csv               # ΔE, ΔH, ΔS, ΔG vs polymer length
      energy_matrix/<section>_<E|H|G>.csv

  Polymer/<monomer_slug>_<length>_<hash>/
    ...                           # same structure as Monomer/ (no summary/)

  Initiator/<initiator_dir>/
    ...                           # same structure as Monomer/ (no summary/)
```

Molecule directories are named `<slug>_<hash>`, where the slug is a
filesystem-safe version of the SMILES and the hash is the first 8 characters of
the SHA-1 of the full SMILES. The hash keeps molecules whose slugs coincide
(e.g. an alkene and the corresponding alkane, since `=` is dropped from the
slug) in separate directories. Polymer directories use the monomer slug, the
chain length and a hash of the polymer SMILES (which includes the initiator).

### Migrating output from older versions

Older versions named directories `<slug>` (and `Polymer/<monomer_slug>_<length>`),
so their output is not found on resume. `larest-migrate` renames them to the
current naming, given the `config.toml` of the old run (old polymer directories
are keyed only by monomer and length, so the reaction type and initiator must
match the old run):

```bash
larest-migrate config.toml output        # dry run: print planned moves
larest-migrate config.toml output --apply
```

Then rerun `larest config.toml -o output`: completed stages are resumed and CENSO
results are re-parsed from the raw CENSO output. Existing directories are never
overwritten. A legacy directory claimed by several configured molecules (e.g. an
alkene and the corresponding alkane, which had the same slug) is not moved, as it
is unclear which molecule its contents belong to; recompute it or move it by hand.
After upgrading an existing install, rerun `pip install .` so the `larest-migrate`
command is available.

`summary/` and `results/<monomer_dir>/` are keyed by monomer only: running the
same monomer with a different reaction type or initiator into the same output
directory overwrites them (the per-molecule results are kept). Use a separate
output directory for each reaction setup.

## Run summary

`larest_summary.txt` is rewritten after each monomer and contains:

- the key settings (reaction, temperature, solvation, levels of theory);
- for each monomer: the monomer check (rings, functional groups, bond opened,
  unassigned stereochemistry), the initiator used, which stages completed for
  each molecule, the reaction thermodynamics for every section, and the energy
  matrices for the final stage;
- warnings, e.g. failed calculations or a poor 1/n extrapolation.

## Energies

| Quantity | Definition |
|---|---|
| `E` | Electronic energy, including implicit solvation where used |
| `H` | Enthalpy: E + thermal (mRRHO) enthalpy correction |
| `G` | Free energy: E + thermal (mRRHO) free energy correction |
| `S` | Entropy: (H − G) / T |

For the RDKit and CREST stages these come from xTB (`TOTAL ENERGY`,
`TOTAL ENTHALPY`, `TOTAL FREE ENERGY`). For the CENSO stages, `E` is the DFT
energy plus the solvation free energy ΔGsolv, `G` is CENSO's `Gtot`, and `H`
uses the enthalpy correction from the same xTB mRRHO calculation CENSO uses
for `G` (`TOTAL ENTHALPY − TOTAL ENERGY` in `xtb_rrho.out`), assuming the
solvation enthalpy equals ΔGsolv. The entropy therefore does not include the
solvation entropy. The CENSO prescreening stage has no mRRHO contribution, so
only `E` is reported for it.

## Sections

| Section | Description |
|---|---|
| `rdkit` | Lowest free energy RDKit conformer (xTB) |
| `crest` | Lowest free energy CREST conformer (xTB) |
| `censo_<stage>` | CENSO sub-stage, Boltzmann average over the conformer ensemble |
| `censo_<stage>_min` | CENSO sub-stage, lowest free energy conformer |
| `censo_corrected` | `censo_refinement_min` with the CREST conformational entropy added (G = H − T·S) |

where `<stage>` is `prescreening`, `screening`, `optimization` or
`refinement`. Boltzmann weights use CENSO's free energies (`Gtot`) at the
configured temperature.

## Summary CSVs

Each `summary/<section>.csv` contains one row per polymer length:

| Column | Description |
|---|---|
| `polymer_length` | Chain length *n* |
| `monomer_E/H/S/G` | Monomer thermodynamic parameters (J/mol, J/mol/K) |
| `initiator_E/H/S/G` | Initiator parameters; zero for RER |
| `polymer_E/H/S/G` | Polymer parameters |
| `delta_E/H/S/G` | Reaction energy per monomer (J/mol, J/mol/K) |

```
delta_param = (polymer_param - n * monomer_param - initiator_param) / n
```

A final row (`polymer_length = inf`) gives each delta extrapolated to infinite
chain length by linear regression against 1/n.

## Energy matrices

For each section and each of E, H and G, `summary/energy_matrix/<section>_<X>.csv`
(kJ/mol) has rows *n* (chain lengths) and columns *m*:

| Entry | Descriptor | Definition |
|---|---|---|
| (n, 0) | X<sup>n</sup><sub>0</sub> | X<sup>n</sup><sub>n</sub> / n, reaction energy per monomer |
| (n, n) | X<sup>n</sup><sub>n</sub> | X(P<sub>n</sub>) − n·X(M) − X(I), total reaction energy |
| (n, m), 0 < m < n | IX<sup>n</sup><sub>m</sub> | (X(P<sub>n</sub>) − X(P<sub>m</sub>)) / (n − m) − X(M), incremental energy per monomer |
| (inf, 0) | X<sup>∞</sup><sub>0</sub> | X<sup>n</sup><sub>0</sub> extrapolated linearly in 1/n; R² in column `R2` |

The initiator term is omitted for RER, and for RER the cyclic chain of degree 1
is the monomer itself. Incremental energies are only available where P<sub>m</sub>
was computed (i.e. *m* is one of the configured lengths, or 1 for RER). They
need no initiator energy, so end-group contributions cancel. A warning is
given if R² < 0.8.

## `results.json`

`results/<monomer_dir>/results.json` contains the monomer check, the initiator,
the raw E/H/S/G of every molecule and section (J/mol, J/mol/K), and all energy
matrices (kJ/mol), with missing values as `null`.

## Log file

`larest.log` is written to the output directory and contains the full pipeline log at DEBUG level. Console output respects the `--verbose` flag.
