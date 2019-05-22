#!/bin/bash -l
#SBATCH -J pyAPES_lettosuo
#SBATCH -o lettosuo_out_%j
#SBATCH -e lettosuo_err_%j
#SBATCH --mail-type=END
#SBATCH --mail-user=kersti.haahti@luke.fi
#SBATCH --mem-per-cpu=500
#SBATCH --time=00:01:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --partition=serial

module load python-env/3.5.3

export OMP_NUM_THREADS=1

# run commands
srun python parallelAPES.py