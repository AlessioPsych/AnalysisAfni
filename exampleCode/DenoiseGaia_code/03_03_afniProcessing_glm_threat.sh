#!/bin/bash
# Set number of tasks to run
#SBATCH --ntasks=1
# Set the number of CPU cores for each task
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
# Walltime format hh:mm:ss
#SBATCH --time=22:00:00
# Output and error files
##SBATCH -o job_denoise.%J.out
##SBATCH -e job_denoise.%J.err

# Load modules here (safety measure)
module purge
module load afni
module load R
module load freesurfer
source /share/apps/NYUAD/miniconda/3-4.11.0/bin/activate
conda activate myenv
export OMP_NUM_THREADS=2

# run script
#Rscript /home/af4887/Programs/Projects/DenoiseGaia/full_28_01_testCode.R
Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/03_03_afniProcessing_glm_threat.R

# to run from shell: screen.sh 'sbatch 03_02_afniProcessing_glm_safe.sh'
