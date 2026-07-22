# Changelog

## Unreleased

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
