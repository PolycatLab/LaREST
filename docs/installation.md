# Installation

## Prerequisites

LaREST requires [conda](https://docs.conda.io) to manage its environment. The `env.yaml` file installs xTB, CREST, and all Python dependencies automatically.

**ORCA** must be installed separately — it is not available via conda. Download it from the [ORCA forum](https://orcaforum.kofo.mpg.de) and ensure the `orca` binary is on your `PATH`, ahead of `/usr/bin` (some Linux systems ship an unrelated `orca` screen reader there). The shared-library ORCA builds need the matching OpenMPI version for parallel runs, e.g. `conda install -c conda-forge "openmpi>=4.1,<5"` for the `openmpi418` builds, with ORCA's directory and OpenMPI's `lib` directory on `LD_LIBRARY_PATH`.

## Steps

```bash
git clone https://github.com/PolycatLab/LaREST.git
cd LaREST
conda env create -f env.yaml
conda activate larest
pip install .
```

This installs Python ≥ 3.12, xTB and CREST (via conda) and the `larest` Python package with CENSO (via pip) into the `larest` conda environment. Python 3.12 or later is required (CENSO 3 needs it).

## Development setup

To contribute to LaREST, install the dev dependencies and register the pre-commit hooks:

```bash
pip install -e ".[dev]"
pre-commit install
```

This installs ruff, ty, pre-commit, and pytest. The hooks run ruff (lint + format) and ty (type checking) automatically before each commit. To run them manually against all files:

```bash
pre-commit run --all-files
```
