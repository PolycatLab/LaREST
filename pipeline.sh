#!/usr/bin/bash
#PBS -N LaREST
#PBS -l walltime=23:59:00
#PBS -l select=1:ncpus=16:mem=512gb:mpiprocs=16

# MODIFY ME
# ---------
CONDA_DIR="${HOME}/anaconda3/bin"
CONDA_ENV="larest"
N_CORES=16
# Config file and output directory. Relative paths are relative to the
# directory the job was submitted from. Both can also be set per job, e.g.
#   qsub -v CONFIG_FILE=configs/gbl.toml,OUTPUT_DIR=runs/gbl pipeline.sh
CONFIG_FILE="${CONFIG_FILE:-config.toml}"
OUTPUT_DIR="${OUTPUT_DIR:-output}"
# ---------

# run from the submission directory (PBS starts jobs in $HOME)
cd "${PBS_O_WORKDIR:-$(pwd)}" || exit 1

# resolve relative paths against the submission directory
[[ "${CONFIG_FILE}" = /* ]] || CONFIG_FILE="$(pwd)/${CONFIG_FILE}"
[[ "${OUTPUT_DIR}" = /* ]] || OUTPUT_DIR="$(pwd)/${OUTPUT_DIR}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "Config file not found: ${CONFIG_FILE}" >&2
    exit 1
fi
echo "Config: ${CONFIG_FILE}"
echo "Output: ${OUTPUT_DIR}"

# xtb-required options
ulimit -s unlimited
ulimit -l unlimited
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS="${N_CORES},1"
export OMP_MAX_ACTIVE_LEVELS=1
export OPENBLAS_NUM_THREADS=1
export XTBPATH="$(pwd)"

# activate conda env
eval "$("${CONDA_DIR}"/conda shell.bash hook)"
conda activate "${CONDA_ENV}"

# load orca and crest (censo is installed in the conda env with larest)
module load tools/prod
module load ORCA/6.1.1-gompi-2024a-avx2 CREST/3.0.2-gfbf-2023b

# create run directory
mkdir -p "${OUTPUT_DIR}"

# run LaREST
larest "${CONFIG_FILE}" -o "${OUTPUT_DIR}"
