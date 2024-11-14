rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
inputDirEpi <- '/scratch/af4887/Proj_Gaia_David/EPI_timeSlicedCorrected'
inputDirBlipData <- '/scratch/af4887/Proj_Gaia_David/rawdata_denoised'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
outputDir <- '/scratch/af4887/Proj_Gaia_David/afni_processed'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- dir( inputDirEpi )
runCodeFlag <- 1

# # clean up output folder, if it exists
# if ( dir.exists( outputDir ) ) {
#   instr <- sprintf('rm -R %s', outputDir )
#   print( instr )
#   if (runCodeFlag==1) { system( instr ) }
#   instr <- sprintf('mkdir %s', outputDir )
#   print( instr )
#   if (runCodeFlag==1) { system( instr ) }
# }else{
#   instr <- sprintf('mkdir %s', outputDir )
#   print( instr )
#   if (runCodeFlag==1) { system( instr ) }
# }

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
var<- args[1]

for ( nSubj in var : var  ) { # nSubj <- 1 length( singleSubjectFolders )
  
  setwd( mainDir )
  
  # make subject specific output folder
  #instr <- sprintf('mkdir %s/%s', outputDir, singleSubjectFolders[ nSubj ] )
  #print( instr )
  #if (runCodeFlag==1) { system( instr ) }

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
                  sprintf('-blocks despike volreg'),
                  sprintf('-blip_forward_dset %s', forwardBlipDataset ),
                  sprintf('-blip_reverse_dset %s', reverseBlipDataset ),
                  sprintf('-volreg_align_to MIN_OUTLIER'),
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

