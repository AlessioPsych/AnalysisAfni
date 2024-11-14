rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders_Freesurfer <- list.dirs( inputFreesurfer, recursive=FALSE )
runCodeFlag <- 1

if ( dir.exists( sprintf('%s', 'dataFrom_3dROIstats' ) ) ) {
  instr <- sprintf('rm -R %s', 'dataFrom_3dROIstats' )
  if (runCodeFlag==1) { system( instr ) }
}
instr <- sprintf('mkdir %s', 'dataFrom_3dROIstats' )
if (runCodeFlag==1) { system( instr ) }
setwd('dataFrom_3dROIstats')
targetDir <- getwd()

for ( nSubj in 1 : 30  ) { # nSubj <- 1 length( singleSubjectFolders )
  
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
  
  currentSubj <- strsplit( singleSubjectFolders_Freesurfer[ nSubj ], '[/]' )[[1]][6]

  instr <- paste( sprintf( 'cp %s*Desikan-Killiany.1D %s', currentSubj, targetDir ) )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  instr <- paste( sprintf( 'cp %s*Destrieux.1D %s', currentSubj, targetDir ) )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
} 

