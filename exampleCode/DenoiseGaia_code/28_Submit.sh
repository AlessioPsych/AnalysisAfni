#!/bin/bash
# Set number of tasks to run
#SBATCH --ntasks=1
# Set the number of CPU cores for each task
#SBATCH --cpus-per-task=28
# Walltime format hh:mm:ss
#SBATCH --time=12:00:00
# Output and error files
#SBATCH -o job_denoise.%J.out
#SBATCH -e job_denoise.%J.err

# Load modules here (safety measure)
module purge
module load afni
module load R
source /share/apps/NYUAD/miniconda/3-4.11.0/bin/activate
conda activate myenv

# run script
#Rscript /home/af4887/Programs/Projects/DenoiseGaia/full_28_01_testCode.R
Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/full_28_01_testCode.R

# to run from shell: screen.sh 'sbatch 28_Submit.sh'
