#!/usr/bin/bash
#PBS -N LaREST
#PBS -l walltime=23:59:00
#PBS -l select=1:ncpus=16:mem=512gb:mpiprocs=16

# MODIFY ME
# ---------
CONDA_DIR="${HOME}/anaconda3/bin"
CONDA_ENV="larest"
N_CORES=16
# ---------

OUTPUT_DIR="output"
CONFIG_FILE="config.toml"

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
conda activate ${CONDA_ENV}

# load orca and censo and crest
module load tools/prod
module load ORCA/6.1.1-gompi-2024a-avx2 CREST/3.0.2-gfbf-2023b

# create run directory
mkdir -p "${PBS_O_WORKDIR}/${OUTPUT_DIR}"

# run LaREST
larest "${PBS_O_WORKDIR}/${CONFIG_FILE}" -o "${PBS_O_WORKDIR}/${OUTPUT_DIR}"
