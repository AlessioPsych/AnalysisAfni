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

for ( nSubj in 1 : length( singleSubjectFolders )  ) {# nSubj <- 1
  
  sprintf('Move to the individual participant folder:')
  setwd( mainDir )
  setwd( singleSubjectFolders[ nSubj ] )  
  sprintf( 'Current folder: %s', print( getwd() ) )
  
  # copy 'anat' folder
  targetFolder <- 'anat'
  pathToCopy <- sprintf( '%s/%s/%s', mainDir, singleSubjectFolders[ nSubj ], targetFolder )
  pathToDestination <- sprintf( '%s/%s/%s', copyDir, singleSubjectFolders[ nSubj ], targetFolder )
  # clean up destination folder, if it exists
  if ( dir.exists( pathToDestination ) ) {
    instr <- sprintf('rm -R %s', pathToDestination )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- sprintf('cp -R %s %s', pathToCopy, pathToDestination )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # copy 'fmap' folder
  targetFolder <- 'fmap'
  pathToCopy <- sprintf( '%s/%s/%s', mainDir, singleSubjectFolders[ nSubj ], targetFolder )
  pathToDestination <- sprintf( '%s/%s/%s', copyDir, singleSubjectFolders[ nSubj ], targetFolder )
  # clean up destination folder, if it exists
  if ( dir.exists( pathToDestination ) ) {
    instr <- sprintf('rm -R %s', pathToDestination )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- sprintf('cp -R %s %s', pathToCopy, pathToDestination )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
}

setwd( mainDir )
sprintf( 'done, current folder: %s', print( getwd() ) )

