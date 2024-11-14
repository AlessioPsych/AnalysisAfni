#!/bin/bash
# Set number of tasks to run
#SBATCH --ntasks=1
# Set the number of CPU cores for each task
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
# Walltime format hh:mm:ss
#SBATCH --array=1
#SBATCH --time=04:00:00
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
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# run script
#Rscript /home/af4887/Programs/Projects/DenoiseGaia/full_28_01_testCode.R
Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/06_afni_glm.R ${SLURM_ARRAY_TASK_ID}

# to run from shell: screen.sh 'sbatch 06_afni_glm.sh'
