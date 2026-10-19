# Pipeline architecture

LaREST chains four computational stages to produce thermodynamic parameters (H, S, G) for ring-opening polymerization reactions. Each stage feeds into the next, progressively increasing the level of theory.

## Stages

### 1. RDKit + xTB

**Purpose:** Generate a set of initial molecular conformers and rank them by free energy.

1. RDKit generates MMFF conformers (default: 50 per molecule).
2. Each conformer is optimised at the xTB GFN2 level and a Hessian is computed to obtain H, S, and G.
3. The lowest-free-energy conformer is passed to the next stage.

Controlled by `[rdkit]` and `[xtb]` in `config.toml`. Set `steps.rdkit = false` to skip.

### 2. CREST conformer ensemble

**Purpose:** Explore the conformational space more thoroughly using metadynamics.

1. CREST's iMTD-GC algorithm generates a conformer/rotamer ensemble from the best RDKit conformer.
2. The ensemble is deduplicated and sorted by energy (CREGEN).
3. xTB re-ranks the ensemble by free energy; the lowest-energy conformer is passed forward.

Controlled by `[crest.confgen]` in `config.toml`. Set `steps.crest_confgen = false` to skip.

### 3. CENSO + ORCA (DFT refinement)

**Purpose:** Refine the CREST ensemble with density functional theory using four sub-stages of increasing accuracy.

| Sub-stage | Label | Default functional | Default basis |
|---|---|---|---|
| Prescreening | `censo_prescreening` | PBE-D4 | def2-SV(P) |
| Screening | `censo_screening` | r2SCAN-3c | def2-TZVP |
| Optimisation | `censo_optimization` | r2SCAN-3c | def2-TZVP |
| Refinement | `censo_refinement` | wB97X-V | def2-TZVP |

Each sub-stage applies an energy window threshold to prune the ensemble before passing it to the next stage. ORCA is used as the QM backend throughout.

For each sub-stage, LaREST reads the per-conformer energies from CENSO's
`<n>_<PART>.json` files and reports both the Boltzmann-averaged values
(`censo_<stage>`) and those of the lowest free energy conformer
(`censo_<stage>_min`). The enthalpy is reconstructed from the xTB mRRHO
calculation CENSO uses for the free energy (see [Output](output.md#energies)).

Controlled by `[censo.*]` in `config.toml`. Set `steps.censo = false` to skip.

### 4. CREST entropy

**Purpose:** Compute the conformational entropy correction using CREST's entropy mode.

1. Starting from the lowest free energy refined CENSO conformer, CREST re-explores the conformational space to obtain a well-converged entropy estimate.
2. Following Grimme and co-workers, the resulting conformational entropy (S_conf + δS_rrho) is added to the entropy of the lowest free energy CENSO refinement conformer to produce the `censo_corrected` section.

If CENSO is disabled, the most refined structure available (CREST, then RDKit) is used as the reference and the raw entropies are stored as `crest_entropy`.

Controlled by `[crest.entropy]` in `config.toml`. Set `steps.crest_entropy = false` to skip.

## Checkpointing

At the start of each molecule's run, LaREST walks through the output directory and identifies the first missing result file. All stages up to that point are skipped; execution resumes from there. This means interrupted runs can be restarted without any manual intervention.

Each completed stage records the settings it was run with in
`<stage>/larest_settings.json`. If those settings have changed (other than
core counts), the stage and all later stages are re-run rather than reusing
stale results. CENSO results are always re-derived from CENSO's raw output when
resuming, so runs from older LaREST versions pick up the current parsing.

## Molecules processed

For each monomer SMILES in the config, LaREST first checks the monomer (see
[Configuration](configuration.md#monomer-check)), skipping it if it cannot be
ring-opened, then runs the full pipeline for:

1. The **monomer** itself.
2. The **initiator** (ROR only), either as configured or the isodesmic initiator for the monomer (`initiator = "auto"`).
3. Each **polymer** at every requested chain length.

Polymer SMILES are constructed automatically from the monomer, initiator and chain length. Final reaction thermodynamics (ΔE, ΔH, ΔS, ΔG) are computed in `compile_results` and written to `summary/` CSVs, the energy matrices to `summary/energy_matrix/`, and an overview of the run to `larest_summary.txt`.
