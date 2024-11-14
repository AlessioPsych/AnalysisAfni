#!/bin/bash
# Set number of tasks to run
#SBATCH --ntasks=1
# Set the number of CPU cores for each task
#SBATCH --cpus-per-task=28
# Walltime format hh:mm:ss
#SBATCH --array=1-30 
#SBATCH --time=10:00:00
# Output and error files
##SBATCH -o job_denoise.%J.out
##SBATCH -e job_denoise.%J.err

# Load modules here (safety measure)
module purge
module load afni
module load R
source /share/apps/NYUAD/miniconda/3-4.11.0/bin/activate
conda activate myenv

# run script
Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/28_parallel_testCode.R ${SLURM_ARRAY_TASK_ID}

# to run from shell: screen.sh 'sbatch parallel_28_Submit.sh'
