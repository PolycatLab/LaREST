# Changelog

## Unreleased

### Features

- **Monomer check** (`chem.analyse_monomer`): rings, ring sizes, functional groups,
  the bond opened and its ring size, and unassigned stereochemistry are reported;
  acyclic monomers or monomers without a supported group are skipped. Bond priority:
  lactone > thionolactone > thiolactone > dithiolactone > lactam > thiolactam
  (new) > carbonates and derivatives. Monomers can be given as
  `{ smiles = "...", group = "..." }` to open a different group.
- **More initiators**: water and esters/thioesters/amides (e.g. methyl acetate) in
  addition to alcohols and formates.
- **`initiator = "auto"`**: the isodesmic initiator is chosen per monomer (methyl
  acetate for lactones, N-methylacetamide for lactams, dimethyl carbonate for
  carbonates, ...), so monomers needing different initiators can be run together.
- **Electronic energies**: `E` is reported alongside H, S and G at every stage.
- **Lowest-conformer CENSO results**: each CENSO stage has a Boltzmann-averaged
  section (`censo_<stage>`) and a lowest free energy conformer section
  (`censo_<stage>_min`).
- **Energy matrices** (`energy_matrix.py`): Eⁿ₀, Eⁿₙ and incremental energies IEⁿₘ for
  E, H and G, with the 1/n extrapolation and its R² (warning below 0.8), written to
  `summary/energy_matrix/`.
- **Run summary and collected results** (`report.py`): `larest_summary.txt` gives an
  overview of the settings, monomer checks, initiators, completed calculations,
  thermodynamics and energy matrices; `results/<monomer>/` collects the final
  structures (`.xyz`), summary CSVs and a machine-readable `results.json`.
- **`[thermo].solvent`**: one setting for implicit solvation in every stage (`"gas"`
  for gas phase), like `[thermo].temperature`.
- **`larest-migrate`** (`migrate.py`): renames molecule directories from older LaREST
  versions (`<slug>`, `<monomer slug>_<n>`) to the current `<slug>_<hash>` naming,
  using the config of the old run, so old runs resume and CENSO results are re-derived
  without re-running DFT. Dry run by default (`--apply` to move); never overwrites, and
  directories shared by several molecules (legacy slug collisions) are left in place
  with a warning.

### Bug Fixes

- **Solution-phase free energies mixed standard states** (`standard_state.py`): xTB and
  CENSO combine 1 atm ideal-gas thermostatistics with 1 M → 1 M solvation free
  energies, so reaction ΔG and ΔS were off by Δn·7.9 kJ/mol and Δn·26.6 J/mol/K (Δn =
  change in the number of molecules, e.g. −1 for ROR with n = 1). With implicit
  solvation, G and S are now reported for a 1 M standard state
  (`[thermo].standard_state`, default `"auto"`).
- **Resume could reuse results built on outdated upstream stages** (`checkpoint.py`):
  each stage's settings stamp now includes the settings of all earlier enabled stages,
  so e.g. a changed RDKit setting invalidates CREST, CENSO and entropy results even if
  a re-run was interrupted. Stamps from the previous version are therefore seen as
  changed once.
- **Stale CENSO sub-stage files were read after a re-run** (`censo.py`): the CENSO
  directory is now cleared before CENSO runs, so e.g. an old `3_REFINEMENT.json` is not
  taken as current when the refinement is disabled.
- **CREST entropy failed, and re-ran on every resume, without CENSO enthalpies**
  (`checkpoint.combine_entropy`): if no CENSO stage has H and S, the CREST entropies
  are stored as `crest_entropy` instead.
- **CREST conformer generation re-ran on every resume with `[steps].xtb = false`**: the
  CREST ensemble is now the checkpoint in that case.
- **Old CENSO structures and sections were reused after disabling a stage**: the final
  structure ignores CENSO output when CENSO is disabled, and sections are no longer
  carried over from an earlier `results.json`.
- **`[thermo].solvent = "gas"` overrode an explicit `[censo.general].solvent`.**
- **`larest-migrate` could attach results to a different polymer**: it now checks that
  the molecule stored in each legacy directory (`rdkit/conformers.sdf`) is the one the
  current version builds, and skips the directory otherwise.
- **R² of the 1/n extrapolation is no longer reported for two chain lengths**, where
  the fit is always exact.
- **CENSO enthalpy and entropy were wrong** (`censo.py`): CENSO's `avE` (the DFT
  electronic energy, without thermal corrections or solvation) was stored as `H`, so
  `S = (H - G) / T = -(ΔGsolv + GmRRHO) / T` was meaningless (typically negative) and
  ΔH was really ΔE. Results are now read per conformer from CENSO's `<n>_<PART>.json`
  files: `E = E(DFT) + ΔGsolv`, `G = Gtot`, and `H = E + H_mRRHO`, with the enthalpy
  correction taken from the same xTB mRRHO calculation CENSO uses for `G`. ΔG is
  unchanged. For the β-propiolactone example (RER, n = 2) ΔH changes from −54.5 to
  −44.5 kJ/mol and ΔS from −128 to −94 J/mol/K. Existing runs are corrected on resume
  without re-running CENSO.
- **Different molecules could share an output directory** (`output.py`, `main.py`):
  `slugify` drops `=`, `#`, `/` and `\`, so e.g. an alkene and the matching alkane (M1
  and M1-H2, or BiL= and BiL) mapped to the same directory and the second silently
  reused the first's checkpointed results. Directories are now `<slug>_<hash>` with a
  hash of the full SMILES; polymer directories hash the polymer SMILES, so changing
  the initiator or reaction type no longer reuses old polymers. **Old output
  directories are not resumed** unless migrated with `larest-migrate` (see Features).
- **Wrong bond broken for α-alkylidene lactones in ROR** (`chem.py`): for e.g.
  α-methylene-γ-butyrolactone the C–C(=O) bond was cleaved instead of the acyl–O
  bond. RER results were unaffected.
- **Unrecognised CENSO settings were silently ignored** (`censo.py`): CENSO drops
  unknown `.censo2rc` keys, and `config/reference.toml` listed several that CENSO 3 no
  longer supports (`bhess`, `scale`, `consider_sym`, `rmsdbias`, `retry_failed`, each
  part's `run`, `implicit`). LaREST now warns about unrecognised `[censo.*]` keys at the
  start of a run (checked against CENSO's own configuration models), and the reference
  config lists the CENSO 3.0.6 keys (sub-stages are enabled through `[censo.cli]`).
- **Carbonate/carbamate ROR products had the initiator on the wrong end**
  (`chem.py`): e.g. trimethylene carbonate + methanol gave a peroxide. Breakable bonds
  are now defined with atom-mapped SMARTS (acyl carbon first, ring bonds only).
- **Changed settings silently reused old results** (`checkpoint.py`): each stage now
  records its settings in `<stage>/larest_settings.json` and is re-run if they change
  (core counts excluded).
- **`censo_corrected` paired the CREST conformational entropy with Boltzmann-averaged
  CENSO values**; it now uses the lowest free energy refinement conformer (the CREST
  entropy reference structure), following Grimme and co-workers.
- **`crest_entropy` failed when CENSO was disabled**; it now falls back to the CREST or
  RDKit best structure.

## v0.0.1

### Testing

- **Added pytest test suite** (`tests/`): 132 unit tests covering all modules except CENSO/ORCA
  execution. Tests are split into fast unit tests (mocked subprocess, no external tools) and
  integration tests that invoke `xtb`/`crest` directly. Run with `pytest tests/`; pass
  `--integration` to include the integration tests.

  Modules covered: `data`, `output`, `setup`, `chem`, `xtb`, `crest`, `censo` (parsing only),
  `checkpoint`, `rdkit` (stage), `main`.

  Sample fixture files added to `tests/data/`: `xtb_output.txt`, `crest_entropy_output.txt`,
  `rdkit_results.csv`, `censo_output.txt`, `censo_conformers.xyz`.

- **`pytest` configuration added to `pyproject.toml`**: `testpaths = ["tests"]` and `-v` set
  as default options so `pytest` can be invoked without arguments from the project root.

### Features

- **Unified `[thermo].temperature` setting** (`setup.py`, `defaults.toml`): a single
  temperature (default 298.15 K) now sets the temperature at which every stage evaluates
  H/S/G and the `T` used to derive `S = (H − G) / T`. It propagates to `[xtb].temperature`,
  `[censo.general].temperature`, `[crest.confgen].temp`, and `[crest.entropy].temp` unless
  overridden per-stage. For xTB the temperature is injected via a `$thermo` xcontrol block
  (`--input`) so xTB evaluates its free energy at that T; `[xtb].etemp` is now the electronic
  (Fermi-smearing) temperature only and no longer doubles as the thermostatistical T.

### Bug Fixes

- **`gas-phase` CENSO setting was silently ignored** (`censo.py`): CENSO's field is
  `gas_phase` (underscore) and its config loader silently drops unrecognised keys, so the
  hyphenated `gas-phase` written to `.censo2rc` never disabled solvation — SMD(solvent) stayed
  active even when gas-phase was requested. `create_censorc` now normalises hyphens to
  underscores when writing `.censo2rc`, and the reference config uses `gas_phase`.

- **`censo_corrected` free energy was left stale after the entropy correction**
  (`checkpoint.py`): `apply_entropy_correction` overwrote the entropy with the CREST
  conformational entropy but left `G` at the CENSO refinement value, so `G ≠ H − T·S`. `G`
  is now recomputed as `H − T·S` from the corrected entropy.

- **`MolToXYZFile` wrote the same conformer to every `.xyz` file** (`base.py`): Missing
  `confId=conformer_id` argument meant all 50 conformer files were identical, making
  downstream xTB ranking meaningless. Fixed by passing the correct conformer ID.

- **`Monomer` crashed on construction for RER reactions** (`base.py`): `Initiator` was
  always instantiated in `Monomer.__init__`, even for RER where `initiator = ""`. Calling
  `get_mol("")` raised `PolymerBuildError` immediately. `Initiator` is now only created
  when `reaction.type == "ROR"`.

- **`Exception.__init__` packed args into a tuple** (`exceptions.py`): `super().__init__(args)`
  passed the whole `*args` tuple as a single argument, so `str(PolymerBuildError("message"))`
  rendered as `"('message',)"` instead of `"message"`. Fixed to `super().__init__(*args)`.

- **Failed external tool runs were silently ignored** (`base.py`): `subprocess.Popen(...).wait()`
  never checked the return code, so a crashed xTB, CREST, or CENSO process would be treated
  as success. Replaced with `subprocess.run(..., check=True)` at all five call sites, which
  raises `CalledProcessError` on non-zero exit.

- **Falsy energy check treated `0.0` as missing** (`parsers.py`): `if xtb_output["H"] and
  xtb_output["G"]` evaluates `0.0` as falsy, so entropy would not be computed if either
  energy happened to be exactly zero. Changed to explicit `is not None` checks in
  `parse_xtb_output` and `parse_censo_output`.

- **`sdfstream` could leak on exception** (`base.py`): The SDF file was opened with a bare
  `open()` and closed manually at the end of the conformer loop. An exception mid-loop would
  skip `sdfstream.close()`. Replaced with a `with` context manager.

### Design / Architecture

- **Removed logger-passing anti-pattern** (`parsers.py`, `chem.py`, `checkpoint.py`,
  `output.py`, `setup.py`): Every helper function previously required a `logger: Logger`
  parameter, polluting all signatures. Each module now owns a
  `logger = logging.getLogger(__name__)` at module level, consistent with standard Python
  practice. All call sites updated.

- **Eliminated double exception logging** (`base.py`, `parsers.py`, `chem.py`, `setup.py`):
  The pervasive pattern of `logger.exception(err); logger.exception("Failed to X"); raise`
  logged two entries (both with tracebacks) per exception, and caused cascading duplicates as
  exceptions bubbled up. Replaced with a single `logger.exception("Failed to X")` per catch
  site, only where additional context is added. Redundant try/except-log-raise wrappers with
  no added context were removed entirely.

- **Removed runtime config mutation** (`base.py`, `setup.py`): `_setup_pipeline()` stored
  `temp_config_dir` into `self._config["temp_config_dir"]`, turning the config dict into a
  hybrid config/state object. The temp directory path is now passed explicitly as a
  `temp_dir: Path` parameter to `create_censorc`, and derived directly from `self._args.config`
  at each call site.

- **Refactored repetitive checkpoint loaders** (`checkpoint.py`): Five near-identical
  `_load_*_results` functions (~275 lines total) reduced to a single
  `_load_stage(path, load_fn, label)` helper plus two small stage-specific parse helpers
  (~80 lines). Each stage now specifies only its file path and parse logic.

- **Removed no-op `Initiator.__init__` override** (`base.py`): The override did nothing
  except call `super().__init__()` with the same arguments. Deleted; `LarestMol.__init__`
  is now used directly.

- **Separated pipeline execution from data model** (`base.py` → `data.py`, `pipeline.py`):
  `LarestMol` and its subclasses carried results state, path management, config, and pipeline
  execution methods all in one place. Split into two concerns: `data.py` contains plain
  dataclasses (`Monomer`, `Initiator`, `Polymer`, `MolResults`) with no IO or side effects;
  `pipeline.py` contains `MolPipeline` which owns execution, path resolution, and checkpointing.
  `MolPipeline(mol, output_dir, config).run()` returns a `MolResults` dataclass.

- **Removed path and config state from molecule classes** (`data.py`): `output_dir`,
  `config_dir`, `verbose`, and `config` fields removed from molecule dataclasses. Paths are
  now computed in `MolPipeline._dir_path()` via a `match` on molecule type. The `config_dir`
  parameter was eliminated entirely — the CENSO temp directory now uses `output_dir / "temp"`.

- **Moved `compile_results` out of `Monomer`** (`main.py`): Result aggregation and CSV writing
  was a method on `Monomer`, mixing post-processing concerns into the data class. Now a
  standalone function in `main.py` that takes explicit `MolResults` arguments.

- **Removed custom exceptions** (`exceptions.py`): `PolymerBuildError` and `NoResultsError`
  were trivial `Exception` subclasses adding no behaviour. Replaced with `ValueError` at all
  raise and catch sites. `exceptions.py` is now empty.

- **Checkpoint restoration moved into `MolPipeline`** (`checkpoint.py`): `restore_results` no
  longer accepts a pre-initialised results dict; it initialises one internally and returns
  `(results, stage)`. Construction of molecule objects no longer triggers filesystem access.

- **Polymer construction moved to `main.py`** (`main.py`, `data.py`): `Monomer.__post_init__`
  previously called `build_polymer` and `restore_results` for every polymer length during
  construction. Polymer SMILES are now built explicitly in `main.py` before constructing
  `Polymer` objects, making construction side-effect free.

### Minor

- **`os.makedirs` replaced with `Path.mkdir`** (`output.py`): `os.makedirs(path, exist_ok=False)`
  with a `FileExistsError` catch replaced by `path.mkdir(parents=True, exist_ok=True)`.
  Removed `import os`.

- **`config` variable shadowing in `parse_command_args`** (`parsers.py`): The loop that
  traversed the config rebound the parameter name `config` at each step. Renamed the
  traversal variable to `cfg`.

- **`summary_headings` moved outside the section loop** (`base.py`): The list was identical
  on every iteration of `for section in self._results`. Moved above the loop.
