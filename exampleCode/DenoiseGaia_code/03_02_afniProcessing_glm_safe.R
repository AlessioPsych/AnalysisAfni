rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
inputDirEpi <- '/scratch/af4887/Proj_Gaia_David/EPI_timeSlicedCorrected_Safe'
inputDirBlipData <- '/scratch/af4887/Proj_Gaia_David/rawdata_denoised'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
outputDir <- '/scratch/af4887/Proj_Gaia_David/afni_processed_Safe'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- dir( inputDirEpi )
runCodeFlag <- 1

# remember, to run this code with slurm (in parallel), the output folder need to be present,
# I cannot create it within the code, otherwise it is going to try to create as many times as 
# the number of parallel processes performed

#args <- commandArgs(trailingOnly = TRUE)
#print(args[1])
#var<- args[1]

#for ( nSubj in var : var  ) { # nSubj <- 1 length( singleSubjectFolders )

for ( nSubj in 26 : 30  ) { # nSubj <- 1 length( singleSubjectFolders )
    
  setwd( mainDir )
  print( getwd() )
  
  # get blip files
  setwd( sprintf( '%s/%s/fmap', inputDirBlipData, singleSubjectFolders[ nSubj ] ) ) 
  print( getwd() )
  filesBlip <- dir( pattern = '*.nii.gz' )
  forwardBlipDataset <- sprintf( '%s/%s/fmap/%s', inputDirBlipData, singleSubjectFolders[ nSubj ], filesBlip[2] )
  reverseBlipDataset <- sprintf( '%s/%s/fmap/%s', inputDirBlipData, singleSubjectFolders[ nSubj ], filesBlip[1] )
  print( sprintf('forward blip dataset: %s', forwardBlipDataset ) )
  print( sprintf('reverse blip dataset: %s', reverseBlipDataset ) )
  
  # define input dataset folder
  dsetsFolder <- sprintf('%s/%s', inputDirEpi, singleSubjectFolders[ nSubj ] )
  
  # convert brain.mgz from freesurfer folder
  setwd( sprintf( '%s/%s/mri', inputFreesurfer, singleSubjectFolders[ nSubj ] ) ) 
  print( getwd() )
  if ( file.exists( 'brain.nii' ) ) {
    print( sprintf('file brain.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm brain.nii' ) }
  }
  instr <- 'mri_convert brain.mgz brain.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # move to the output folder
  setwd( outputDir )
  print( getwd() )
  
  instr <- paste( 'afni afni_proc.py',
                  sprintf('-subj_id %s', singleSubjectFolders[ nSubj ] ),
                  sprintf('-copy_anat %s/%s/mri/brain.nii', inputFreesurfer, singleSubjectFolders[ nSubj ] ),
                  sprintf('-dsets %s/*.nii.gz', dsetsFolder ),
                  sprintf('-blocks despike volreg blur scale regress'),
                  sprintf('-blip_forward_dset %s', forwardBlipDataset ),
                  sprintf('-blip_reverse_dset %s', reverseBlipDataset ),
                  sprintf('-volreg_align_to MIN_OUTLIER'),
                  '-blur_size 4',
                  sprintf('-regress_stim_times %s/stimTiming/*.1D', dsetsFolder),
                  '-regress_local_times',
                  '-regress_polort 5',
                  '-regress_stim_types AM2',
                  '-regress_stim_labels anx_rat color fix mem_resp mem_test T1 T2',
                  sprintf('-regress_basis \u0027dmUBLOCK\u0027'),
                  '-regress_opts_3dD -jobs 2',
                  #'-mask_apply epi',
                  '-regress_compute_fitts',
                  '-regress_run_clustsim no',
                  sprintf('-scr_overwrite'),
                  sprintf('-no_epi_review'),
                  sprintf('-execute') 
  )
  
  # instr <- paste( 'afni afni_proc.py',
  #                 sprintf('-subj_id %s', singleSubjectFolders[ nSubj ] ),
  #                 sprintf('-copy_anat %s/%s/mri/brain.nii', inputFreesurfer, singleSubjectFolders[ nSubj ] ),
  #                 sprintf('-dsets %s/*.nii.gz', dsetsFolder ),
  #                 sprintf('-blocks despike align tlrc volreg'),
  #                 sprintf('-blip_forward_dset %s', forwardBlipDataset ),
  #                 sprintf('-blip_reverse_dset %s', reverseBlipDataset ),
  #                 sprintf('-tlrc_base MNI152_2009_template.nii.gz'),
  #                 sprintf('-volreg_align_to MIN_OUTLIER'),
  #                 sprintf('-volreg_align_e2a'),
  #                 sprintf('-execute')
  #          )
  
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
} 


# possibility, convert brain from free surfer folder, then remove to keep freesurfer folder clear,
# then move atlases from freesurfer folder

#@SUMA_Make_Spec_FS -NIFTI -fspath %s/FreeSeg_result -sid FreeSeg_result

