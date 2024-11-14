rm( list=ls() ); gc();

# need to enter afni singularity, load freesurfer and load python beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders <- dir( inputFreesurfer )
runCodeFlag <- 1

#args <- commandArgs(trailingOnly = TRUE)
#print(args[1])
#var<- args[1]

for ( nSubj in 15 : 20  ) { # nSubj <- 1 length( singleSubjectFolders )
  
  setwd( mainDir )
  
  # check subject specific output folder
  folderToCheck <- sprintf('%s/freesurfer/%s/SUMA', mainDir, singleSubjectFolders[ nSubj ] )
  print( sprintf( 'folder to check: %s', folderToCheck ) )
  if ( dir.exists( folderToCheck ) ) {
    instr <- sprintf('rm -R %s', folderToCheck )
    print( sprintf( 'run instruction: %s', instr ) )
    if (runCodeFlag==1) { system( instr ) }
  }

  fsPath <- sprintf('%s/freesurfer/%s', mainDir, singleSubjectFolders[ nSubj ] )
  instr <- sprintf('@SUMA_Make_Spec_FS -NIFTI -fspath %s -sid FreeSeg_results', fsPath )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
} 


