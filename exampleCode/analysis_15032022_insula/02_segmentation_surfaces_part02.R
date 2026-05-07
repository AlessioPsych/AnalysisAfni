rm( list=ls() )
mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
partDir <- 'part2/'

setwd( mainDir )
setwd( partDir )
getwd()
participantsFolders <- dir()
for ( participantsFoldersIndex in seq( 1,length( participantsFolders ) ) ) { #
  
  setwd( mainDir )
  setwd( partDir )
  print( '...' )
  print( '...' )
  print( sprintf( '%s...', participantsFolders[ participantsFoldersIndex ] ) )
  print( '...' )
  print( '...' )
  setwd( sprintf( '%s', participantsFolders[ participantsFoldersIndex ] ) )
  setwd( dir()[1] )
  setwd( dir()[1] )
  print( dir() )
  
  # clean up resamped volume ad freesurfer segmentation
  if ( file.exists( 'anatCopy_1mm.nii.gz' ) ) { system('rm anatCopy_1mm.nii.gz') }
  if ( file.exists( 'fsaverage' ) ) { system('rm fsaverage') }
  if ( dir.exists( 'AHEAD_test' ) ) { system('rm -R AHEAD_test') }
  
  # get t1w file
  participantFiles <- dir()
  fileIdx <- grep( 't1w', participantFiles )
  print( '...' )
  participantT1w <- participantFiles[ fileIdx ]; print( participantT1w )
  print( '...' )
  
  #resample t1w file to 1mm iso
  print( '...' )
  instr <- sprintf( '3dresample -dxyz 1 1 1 -prefix anatCopy_1mm.nii.gz -rmode Linear -input %s', participantT1w ); print( instr ); system( instr )
  print( '...' )
  
  #redefine SUBJECTS_DIR to current directory
  currentDir <- getwd()
  print( '...' )
  Sys.setenv( SUBJECTS_DIR = paste( currentDir ) )
  print( Sys.getenv( 'SUBJECTS_DIR') )
  print( '...' )

  #recon-all instruction
  print( '...' )
  instr <- sprintf( 'recon-all -subjid AHEAD_test -i anatCopy_1mm.nii.gz -all -parallel -openmp 8' ); print( instr ); system( instr )
  print( '...' )
  
  #suma instruction
  print( '...' )
  instr <- sprintf( '@SUMA_Make_Spec_FS -NIFTI -fspath %s/AHEAD_test -sid AHEAD_test', currentDir ); print( instr ); system( instr )
  print( '...' )
  
}


