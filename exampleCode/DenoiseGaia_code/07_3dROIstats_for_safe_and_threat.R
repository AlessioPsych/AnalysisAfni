rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
dataDirSafe <- '/scratch/af4887/Proj_Gaia_David/afni_processed_Safe'
dataDirThreat <- '/scratch/af4887/Proj_Gaia_David/afni_processed_Threat'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders_Safe <- list.dirs( dataDirSafe, recursive = FALSE )
singleSubjectFolders_Threat <- list.dirs( dataDirThreat, recursive = FALSE )
singleSubjectFolders_Freesurfer <- list.dirs( inputFreesurfer, recursive=FALSE )
runCodeFlag <- 1

for ( nSubj in 1 : 30  ) { # nSubj <- 1 length( singleSubjectFolders )
  
  setwd( mainDir )
  print( getwd() )
  
  # define input dataset folder
  dsetsFolder <- sprintf('%s/SUMA', singleSubjectFolders_Freesurfer[ nSubj ] )
  setwd( dsetsFolder )
  print( '###################' )
  print( '###################' )
  print( '###################' )
  print( getwd() )
  print( '###################' )
  print( '###################' )
  print( '###################' )
  
  currentSubjFolder <- strsplit( singleSubjectFolders_Safe[ nSubj ], '[/]' )[[1]][6]
  currentSubj <- strsplit( currentSubjFolder, '[.]' )[[1]][1]
  
  # upsample stats to anatomy (For Desikan-Killiany Atlas, Destrieux Atlas)
  if ( file.exists( sprintf('%s/SUMA/stats_safe_upsampled.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats_safe_upsampled.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  if ( file.exists( sprintf('%s/SUMA/stats_threat_upsampled.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats_threat_upsampled.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  
  # upsample safe
  instr <- paste('afni 3dresample', 
                 '-inset stats_safe.nii.gz',
                 '-prefix stats_safe_upsampled.nii.gz',
                 '-master aparc+aseg_REN_all.nii.gz',
                 '-rmode NN'
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # upsample threat
  instr <- paste('afni 3dresample', 
                 '-inset stats_threat.nii.gz',
                 '-prefix stats_threat_upsampled.nii.gz',
                 '-master aparc+aseg_REN_all.nii.gz',
                 '-rmode NN'
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # stats safe, Desikan-Killiany
  instr <- paste('afni 3dROIstats', 
                 '-mask aparc+aseg_REN_all.nii.gz',
                 '-nzmean',
                 '-nzmedian',
                 '-nzsigma',
                 sprintf( 'stats_safe_upsampled.nii.gz > %s_stats_safe_Desikan-Killiany.1D', currentSubj  )
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # stats safe, Destrieux
  instr <- paste('afni 3dROIstats', 
                 '-mask aparc.a2009s+aseg_REN_all.nii.gz',
                 '-nzmean',
                 '-nzmedian',
                 '-nzsigma',
                 sprintf( 'stats_safe_upsampled.nii.gz > %s_stats_safe_Destrieux.1D', currentSubj )
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # stats threat, Desikan-Killiany
  instr <- paste('afni 3dROIstats', 
                 '-mask aparc+aseg_REN_all.nii.gz',
                 '-nzmean',
                 '-nzmedian',
                 '-nzsigma',
                 sprintf( 'stats_threat_upsampled.nii.gz > %s_stats_threat_Desikan-Killiany.1D', currentSubj )
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # stats threat, Destrieux
  instr <- paste('afni 3dROIstats', 
                 '-mask aparc.a2009s+aseg_REN_all.nii.gz',
                 '-nzmean',
                 '-nzmedian',
                 '-nzsigma',
                 sprintf( 'stats_threat_upsampled.nii.gz > %s_stats_threat_Destrieux.1D', currentSubj )
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # clean up
  instr <- sprintf('rm stats_safe_upsampled.nii.gz' ) 
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('rm stats_threat_upsampled.nii.gz' )
  if (runCodeFlag==1) { system( instr ) }
  
} 

