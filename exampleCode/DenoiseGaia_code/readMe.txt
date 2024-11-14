step 01:
  denoise data using dwidenoise. There are two ways: 
  1. slower, ~3 days: run code 28_submit.sh via shell, using the following call:
    screen.sh 'sbatch 28_Submit.sh'
    this code runs each participant iteratively, using 28 threads in the nthreads field for dwidenoise
  2. USE THIS APPROACH: faster, ~2 hours: run code parallel_28_Submit.sh via the shell, using the following call:
    screen.sh 'sbatch parallel_28_Submit.sh'
    this code runs each participant in parallel, using slurm. Each participant and file uses 28 threads in the nthreads         field for dwidenoise
  
step 02: - update, we need this step, but we are not running fmriprep, just afni, see below
  copies the remaining folders (anat and fmap) from rawdata folder to rawdata_denoised data. This way we can run the whole    preprocessing pipeline using fMRI prep afterwards. This step uses the simple R code: 
  02_copy_anat_func.R
  
step 03:
  start running time slice correction, This step uses the simple R code: 
  03_time_slice_correction.R. This code can be called from the shell, using the following call:
  screen.sh 'sbatch 03_time_slice_correction.sh'

step 03_01:
  for each participant, split the data into 'threat' and 'safe' folders, then the prepocessing will be performed on these separate folders. Data will be coupled with stim timing information (wetransfer_stimTiming)
  simply run the R code, no need for screen.sh and sbatch commands 
  03_01_splitDataintoSafeThreat.R

step 03_02:
  perform despike, blip up - down, and motcorr and glm based on 'Safe' condition using afni_proc_py. 
  remember to CLEAN the folder 'afni_processed_Safe' before running it
  This code can be called from the shell, using the following call:
  screen.sh 'sbatch 03_02_afniProcessing_glm_safe.sh' 
  
step 03_03:
  perform despike, blip up - down, and motcorr and glm based on 'Threat' condition using afni_proc_py. 
  remember to CLEAN the folder 'afni_processed_Threat' before running it
  This code can be called from the shell, using the following call:
  screen.sh 'sbatch 03_02_afniProcessing_glm_threat.sh' 

step 04_01:
  coregistration, safe condition.
  This code can be run via jubail terminal, module load afni, freesurfer and R, no need for screen.sh from the      shell, just run 'Rscript 04_01_afni_coregistration_safe.R'

step 04_02:
  coregistration, threat condition.
  This code can be run via jubail terminal, module load afni, freesurfer and R, no need for screen.sh from the      shell, just run 'Rscript 04_02_afni_coregistration_threat.R'

step05: 
  derive suma folders from freesurfer output - This step uses the simple R code: 
  05_afni_processing_suma.R. 
  We can run ~10 participants at a time, to limit the amount of time requested for the interactive Jubail Desktop
  This code can be called from the Jubail Desktop shell (not the login shell), using the following call: 
  
  singularity shell /share/apps/NYUAD/afni/Ubuntu.sif/
  module load freesurfer
  module load python
  module load R
  Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/05_afni_processing_suma.R

step 06
  bring stats into suma folderd for each individual participant;
  convert stats volume into .niml and .1D surfaces
  run file: 06_bringStatsOnSuma_suma.R from Jubail terminal directly:
  module load freesurfer
  module load python
  module load R
  Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/06_bringStatsOnSuma_suma.R
  
step 07
  dump average ROI stats into 1D files for anova analyses  
  run file: 07_3dROIstats_for_safe_and_threat.sh from shell:
  screen.sh 'sbatch  07_3dROIstats_for_safe_and_threat.sh' 

step 08
  brings the stats into a single folder for anova analysis, it creates the folder: dataFrom_3dROIstats
  run file: 07_3dROIstats_for_safe_and_threat.R from Rstudio itself;




















step04: 
  perform despike, blip up - down, and motcorr using afni_proc_py. 
  This code can be called from the shell, using the following call:
  screen.sh 'sbatch 04_afni_processing.sh' 
  ##### check participant 13, run 02, only 10 dynamics, removed  from the EPI_timeSlicedCorrected folder and moved to the   folder participant13_run02_tooShort_epiSliceTimingCorrect ####
  
step05: 
  derive suma folders from freesurfer output - This step uses the simple R code: 
  05_afni_processing_suma.R. 
  We can run ~10 participants at a time, to limit the amount of time requested for the interactive Jubail Desktop
  This code can be called from the Jubail Desktop shell (not the login shell), using the following call: 
  
  singularity shell /share/apps/NYUAD/afni/Ubuntu.sif/
  module load freesurfer
  module load python
  module load R
  Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/05_afni_processing_suma.R
  
step06: 
  run glm analysis based on Gaia's timings;
  use code 06_afni_glm.R, need to test it, it requires a wrap with sh to run on batch
  
  it uses
  module load afni
  module load R
  Rscript /scratch/af4887/Proj_Gaia_David/DenoiseGaia_code/06_afni_glm.R


  
