args <- commandArgs(T)
print( args )

dataDir <- args[1]
fileAnatomy <- args[2]
nCores <- as.numeric( args[3] )
runFlag <- as.numeric( args[4] )

setwd( dataDir )

##### get all directories ####
allDirsToAnalyze_full <- list.dirs( path = dataDir )
allDirsToAnalyze <- allDirsToAnalyze_full[ 2:length( allDirsToAnalyze_full  ) ] 

##### run freesurfer and suma ####
for ( i in  1:length( allDirsToAnalyze ) ) { #i <- 1
  setwd( allDirsToAnalyze[ i ] )

  #redefine SUBJECTS_DIR to current directory
  currentDir <- getwd()
  print( currentDir )
  print( '...' )
  Sys.setenv( SUBJECTS_DIR = paste( currentDir ) )
  print( Sys.getenv( 'SUBJECTS_DIR') )
  print( '...' )
  
  #recon-all instruction
  print( '...' )
  instr <- sprintf( 'recon-all -subjid FreeSeg_result -i %s -all -parallel -openmp %1.0f', fileAnatomy, nCores )
  print( instr )
  if ( runFlag==1 ) {
    system( instr )
  }
  print( '...' )
  
  #suma instruction
  print( '...' )
  instr <- sprintf( '@SUMA_Make_Spec_FS -NIFTI -fspath %s/FreeSeg_result -sid FreeSeg_result', currentDir )
  print( instr )
  if ( runFlag==1 ) {
    system( instr )
  }
  print( '...' )
  
}
