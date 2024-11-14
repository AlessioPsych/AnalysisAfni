rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
dataDir <- '/scratch/af4887/Proj_Gaia_David/afni_processed_Safe'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- list.dirs( dataDir, recursive = FALSE )
runCodeFlag <- 1

for ( nSubj in 1 : 30  ) { # nSubj <- 1 length( singleSubjectFolders )
  
  setwd( mainDir )
  print( getwd() )
  
  # define input dataset folder
  dsetsFolder <- sprintf('%s', singleSubjectFolders[ nSubj ] )
  setwd( dsetsFolder )
  print( '###################' )
  print( '###################' )
  print( '###################' )
  print( getwd() )
  print( '###################' )
  print( '###################' )
  print( '###################' )
        
  targetAnatomy_head <- 'brain+orig.HEAD'
  startingEpi_head <- 'vr_base_min_outlier+orig.HEAD'
  targetAnatomy_brik <- 'brain+orig.BRIK'
  startingEpi_brik <- 'vr_base_min_outlier+orig.BRIK'
  
  if ( dir.exists('Coregistration') ) {
    instr <- sprintf('rm -R Coregistration')
    print( sprintf('cleaning coregistration folder...') )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- sprintf('mkdir Coregistration')
  print( sprintf('make coregistration folder...') )
  if (runCodeFlag==1) { system( instr ) }

  instr <- sprintf('cp %s Coregistration/%s', targetAnatomy_head, targetAnatomy_head)
  print( sprintf('copy anatomy file head into coregistration folder...') )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('cp %s Coregistration/%s', startingEpi_head, startingEpi_head)
  print( sprintf('copy starting epi file head into coregistration folder...') )
  if (runCodeFlag==1) { system( instr ) }
  
  instr <- sprintf('cp %s Coregistration/%s', targetAnatomy_brik, targetAnatomy_brik)
  print( sprintf('copy anatomy file brik into coregistration folder...') )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('cp %s Coregistration/%s', startingEpi_brik, startingEpi_brik)
  print( sprintf('copy starting epi file brik into coregistration folder...') )
  if (runCodeFlag==1) { system( instr ) }
  
  setwd('Coregistration')
  print( getwd() )
  currentCoregFolder <- getwd()
  
  instr <- paste( 'afni align_epi_anat.py', 
                    '-epi_base 0', 
                    sprintf('-epi %s', startingEpi_brik ), 
                    sprintf('-anat %s', targetAnatomy_brik ), 
                    '-epi2anat',
                    sprintf('-master_epi %s', targetAnatomy_brik ),
                    '-anat_has_skull no',
                    '-epi_strip None',
                    '-cmass cmass'
    )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  instr <- 'afni cat_matvec -ONELINE vr_base_min_outlier_al_mat.aff12.1D -I > coregMat.1D'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  instr <- paste('afni 3dAllineate',
                 sprintf('-base %s', targetAnatomy_brik ), 
                 sprintf('-input %s', startingEpi_brik ), 
                 '-prefix singleShot01.nii.gz',
                 '-1Dmatrix_apply coregMat.1D',
                 '-final linear' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  instr <- paste('afni 3dAllineate',
                 sprintf('-base %s', targetAnatomy_brik ), 
                 sprintf('-input %s', startingEpi_brik ), 
                 '-prefix singleShot02.nii.gz',
                 '-1Dmatrix_apply vr_base_min_outlier_al_mat.aff12.1D',
                 '-final linear' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  

  # convert from freesurfer folder
  currentSubjFolder <- strsplit( singleSubjectFolders[ nSubj ], '[/]' )[[1]][6]
  currentSubj <- strsplit( currentSubjFolder, '[.]' )[[1]][1]
  setwd( sprintf( '%s/%s/mri', inputFreesurfer, currentSubj ) ) 
  print( getwd() )
  if ( file.exists( 'aparc.a2009s+aseg.nii' ) ) {
    print( sprintf('file aparc.a2009s+aseg.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm aparc.a2009s+aseg.nii' ) }
  }
  instr <- 'mri_convert aparc.a2009s+aseg.mgz aparc.a2009s+aseg.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv %s %s/%s', 'aparc.a2009s+aseg.nii', currentCoregFolder, 'aparc.a2009s+aseg.nii' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  if ( file.exists( 'aparc.DKTatlas+aseg.nii' ) ) {
    print( sprintf('file aparc.DKTatlas+aseg.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm aparc.DKTatlas+aseg.nii' ) }
  }
  instr <- 'mri_convert aparc.DKTatlas+aseg.mgz aparc.DKTatlas+aseg.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv %s %s/%s', 'aparc.DKTatlas+aseg.nii', currentCoregFolder, 'aparc.DKTatlas+aseg.nii' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  if ( file.exists( 'aparc+aseg.nii' ) ) {
    print( sprintf('file aparc+aseg.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm aparc+aseg.nii' ) }
  }
  instr <- 'mri_convert aparc+aseg.mgz aparc+aseg.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv %s %s/%s', 'aparc+aseg.nii', currentCoregFolder, 'aparc+aseg.nii' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  if ( file.exists( 'aparc+aseg.nii' ) ) {
    print( sprintf('file aparc+aseg.nii exists, remove...') )
    if (runCodeFlag==1) { system( 'rm aparc+aseg.nii' ) }
  }
  instr <- 'mri_convert aseg.mgz aseg.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mv %s %s/%s', 'aseg.nii', currentCoregFolder, 'aseg.nii' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
} 

