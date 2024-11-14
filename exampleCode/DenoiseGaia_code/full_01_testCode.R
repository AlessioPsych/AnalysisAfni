rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David/rawdata'
copyDir <- '/scratch/af4887/Proj_Gaia_David/rawdata_denoised'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- dir()
runCodeFlag <- 1
flagResample <- 0

for ( nSubj in 1 : 1  ) {# length( singleSubjectFolders ) # nSubj <- 1
  
  sprintf('Move to the individual participant folder:')
  setwd( mainDir )
  setwd( singleSubjectFolders[ nSubj ] )  
  sprintf( 'Current folder: %s', print( getwd() ) )
  setwd( 'func' )
  json_files_toCopy <- dir( pattern='*sbref.json' )
  nifti_files_toCopy <- dir( patter='*sbref.nii*' )
  json_files_toProcess <- dir( pattern='*bold.json' )
  nifti_files_toProcess <- dir( patter='*bold.nii*' )
  #json_files_toProcess <- json_files_toProcess[1:3]
  #nifti_files_toProcess <- nifti_files_toProcess[1:3]
  
  # copy also directories anat and fmap in rawdata_denoise #
  
  # clean up target folder and generate output folder
  if ( dir.exists( sprintf('%s/%s', copyDir, singleSubjectFolders[ nSubj ] ) ) ) {
    instr <- sprintf('rm -R %s/%s', copyDir, singleSubjectFolders[ nSubj ] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  
  # generate output folders (rawdata_denoised/subjID/func)
  instr <- sprintf('mkdir %s/%s', copyDir, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- sprintf('mkdir %s/%s/func', copyDir, singleSubjectFolders[ nSubj ] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # copy nifty files that just need to be copied to (rawdata_denoised/subjID/func)
  for ( nFilesToCopy in 1 : length( nifti_files_toCopy )  ) { # nFilesToCopy <- 1
    print( sprintf('nifti file to copy:') )
    print( sprintf( '%s', nifti_files_toCopy[ nFilesToCopy ] ) )
    targetFile <- sprintf('%s/%s/func/%s', copyDir,  singleSubjectFolders[ nSubj ], nifti_files_toCopy[ nFilesToCopy ] )
    print( targetFile )
    instr <- sprintf('cp %s %s', nifti_files_toCopy[ nFilesToCopy ], targetFile )
    print( instr ) 
    if (runCodeFlag==1) { system( instr ) }
  }

  # copy json files that just need to be copied to (rawdata_denoised/subjID/func)
  for ( nFilesToCopy in 1 : length( json_files_toCopy )  ) { # nFilesToCopy <- 1
    print( sprintf('json file to copy:') )
    print( sprintf( '%s', json_files_toCopy[ nFilesToCopy ] ) )
    targetFile <- sprintf('%s/%s/func/%s', copyDir,  singleSubjectFolders[ nSubj ], json_files_toCopy[ nFilesToCopy ] )
    print( targetFile )
    instr <- sprintf('cp %s %s', json_files_toCopy[ nFilesToCopy ], targetFile )
    print( instr ) 
    if (runCodeFlag==1) { system( instr ) }
  }

  # copy json files, from files that need to be processed in (rawdata_denoised/subjID/func) 
  for ( nFilesToCopy in 1 : length( json_files_toProcess )  ) { # nFilesToCopy <- 1
    print( sprintf('json file to copy (from files to process):') )
    print( sprintf( '%s', json_files_toProcess[ nFilesToCopy ] ) )
    targetFile <- sprintf('%s/%s/func/%s', copyDir,  singleSubjectFolders[ nSubj ], json_files_toProcess[ nFilesToCopy ] )
    print( targetFile )
    instr <- sprintf('cp %s %s', json_files_toProcess[ nFilesToCopy ], targetFile )
    print( instr ) 
    if (runCodeFlag==1) { system( instr ) }
  }

  # denoise and copy nifti files
  for ( nFilesToProcess in 1 : length( nifti_files_toProcess ) ) { # nFilesToProcess <- 1
    if ( file.exists( 'ttt_denoised.nii.gz' ) ) {
      instr <- 'rm ttt_denoised.nii.gz'
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
    print( sprintf('file to process:') )
    print( sprintf( '%s', nifti_files_toProcess[ nFilesToProcess ] ) )
    
    splitNiftiOutput <- strsplit( nifti_files_toProcess[ nFilesToProcess ], '-' )
    niftiOutputName <- paste( splitNiftiOutput[[1]][1], '-',
                              splitNiftiOutput[[1]][2],
                              '-denoised-',
                              splitNiftiOutput[[1]][3], '-',
                              splitNiftiOutput[[1]][4], sep='' )

    # 3dresample instruction to test and speed up things
    inputToDenoise <- 'ttt_inputToDenoise.nii.gz'
    if ( file.exists( inputToDenoise ) ) {
      instr <- sprintf( 'rm %s', inputToDenoise )
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
    if (flagResample==1) {
    instr <- sprintf( 'afni 3dresample -dxyz 20 20 20 -prefix %s -input %s -rmode NN',
                      inputToDenoise,
                      nifti_files_toProcess[ nFilesToProcess ] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    }
    if (flagResample==0) {
    instr <- sprintf( 'afni 3dresample -prefix %s -input %s -rmode NN',
                      inputToDenoise,
                      nifti_files_toProcess[ nFilesToProcess ] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    }
    # 3dautomask instruction to test and speed up things
    inputToDenoiseAutomask <- 'ttt_inputToDenoiseAutomask.nii.gz'
    if ( file.exists( inputToDenoiseAutomask ) ) {
      instr <- sprintf( 'rm %s', inputToDenoiseAutomask )
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
    instr <- sprintf( 'afni 3dAutomask -apply_prefix %s %s',
                      inputToDenoiseAutomask,
                      inputToDenoise )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # denoise instruction
    instr <- sprintf('dwidenoise -nthreads 4 %s %s/%s/func/%s', 
                     inputToDenoiseAutomask, 
                     copyDir, 
                     singleSubjectFolders[ nSubj ], 
                     nifti_files_toProcess[ nFilesToProcess ] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # clean up
    if ( file.exists( inputToDenoise ) ) {
      instr <- sprintf( 'rm %s', inputToDenoise )
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
    if ( file.exists( inputToDenoiseAutomask ) ) {
      instr <- sprintf( 'rm %s', inputToDenoiseAutomask )
      print( instr )
      if (runCodeFlag==1) { system( instr ) }
    }
  }  
}

setwd( mainDir )
sprintf( 'done, current folder: %s', print( getwd() ) )

