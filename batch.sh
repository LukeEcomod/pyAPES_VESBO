#!/bin/bash -l
# created: Apr 26, 2019 1:28 PM
# author: haahtike
#SBATCH -J pyAPES_lettosuo
#SBATCH -o lettosuo_out_%j
#SBATCH -e lettosuo_err_%j
#SBATCH --mail-type=END
#SBATCH --mail-user=kersti.haahti@luke.fi
#SBATCH --mem-per-cpu=3500
#SBATCH --time=17:00:00
#SBATCH --ntasks=1
#SBATCH --partition=serial

module load python-env/3.5.3

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

# run commands
srun python pyAPES.py

seff $SLURM_JOBID
