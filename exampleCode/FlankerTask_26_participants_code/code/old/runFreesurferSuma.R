rm( list=ls() )
mainDir <- '/media/alessiofracasso/DATADRIVE1/Flanker'
setwd( mainDir )
outputDir <- 'Freesurfer_output'

partDirs <- list.dirs(full.names = FALSE, recursive = FALSE)
partDirs <- partDirs[grepl("^*sub-", partDirs)]
print( 'participant dirs...' )
print(partDirs)

print( 'clean up and create main output folder in derivatives...' )
dirToCheck <- sprintf('%s/derivatives/%s', mainDir, outputDir)
flagDir <- dir.exists( dirToCheck  )
if ( flagDir==TRUE ) { 
  #system( sprintf('rm -R %s', dirToCheck ) )
  #dir.create( dirToCheck )
}
if ( flagDir==FALSE ) { dir.create( dirToCheck ) }

nCores <- 6
runFlag <- 1

##### run freesurfer and suma ####
for ( i in  7:length( partDirs ) ) { # i <- 1 length( partDirs )
  setwd( mainDir )
  setwd( partDirs[ i ] )
  setwd( 'anat' )
  print( getwd() )
  anatNii <- list.files(pattern='*.nii.gz')
  pathToFileAnatomy <- sprintf('%s/%s', getwd(),anatNii)  
  
  print( sprintf('make anat folder in derivatives, participant: %s', partDirs[ i ] ) )
  partOutputDir <- sprintf('%s/derivatives/%s/%s', mainDir, outputDir, partDirs[ i ] )
  print( sprintf('mkdir %s', partOutputDir ) )
  system( sprintf('mkdir %s', partOutputDir ) )
  
  setwd( partOutputDir )
  #redefine SUBJECTS_DIR to output directory
  currentDir <- getwd()
  print( currentDir )
  print( '...' )
  Sys.setenv( SUBJECTS_DIR = paste( currentDir ) )
  print( Sys.getenv( 'SUBJECTS_DIR') )
  print( '...' )
  
  #recon-all instruction
  print( '...' )
  instr <- sprintf( 'recon-all -subjid Freesurfer_result -i %s -all -parallel -openmp %1.0f', pathToFileAnatomy, nCores )
  print( instr )
  if ( runFlag==1 ) {
    system( instr )
  }
  print( '...' )
  
  #suma instruction
  print( '...' )
  instr <- sprintf( '@SUMA_Make_Spec_FS -NIFTI -fspath %s/Freesurfer_result -sid Freesurfer_result', currentDir )
  print( instr )
  if ( runFlag==1 ) {
    system( instr )
  }
  print( '...' )
  
}
