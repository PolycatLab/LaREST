# Usage

## Running locally

```bash
larest [config.toml] -o <output_dir>
```

| Argument | Description |
|---|---|
| `config` / `-c` / `--config` | Path to the config file, as the first argument or with `-c` (default: `config.toml` in the current directory) |
| `-o` / `--output` | Directory where all results are written (default: `./output`) |
| `-v` / `--verbose` | Increase console log verbosity |

**Example** — run the example config and write results to `./output`:

```bash
larest config/example.toml -o output/
```

## Running on HPC (Imperial College PBS)

A PBS job script (`pipeline.sh`) is included for Imperial College's HPC service. Edit the resource directives and the `MODIFY ME` block at the top of the file, then submit from the LaREST root directory:

```bash
qsub pipeline.sh
```

Key variables to configure in `pipeline.sh`:

| Variable | Description |
|---|---|
| `CONDA_DIR` | Path to the directory containing your conda binary |
| `CONDA_ENV` | Conda environment name (default: `larest`) |
| `N_CORES` | Number of cores, matching the PBS `ncpus` directive |

## Checkpointing

LaREST automatically checkpoints after each pipeline stage. If a run is interrupted, re-running the same command resumes from the first incomplete stage — completed stages are not re-executed.

Checkpoint state is inferred from the presence of result files on disk; no separate checkpoint file is written. To force a full re-run, delete the relevant molecule directory under the output folder.
