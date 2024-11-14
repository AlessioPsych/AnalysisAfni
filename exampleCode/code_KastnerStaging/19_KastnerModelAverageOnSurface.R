rm(list=ls())

# toolbox and functions
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

# folders 
mainFolder <- '/analyse/Project0226/KastnerModel/staging_area_Kastner'
surfaceFolder <- '/analyse/Project0226/KastnerModel/staging_area_Kastner/anatomies_KastnerClassic_Freesurfer/surfaceAtlases/suma_MNI152_2009'
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
testSubjectDir <- '/analyse/Project0226/KastnerModel/staging_area_Kastner/anatomies_KastnerClassic_Freesurfer/AFH28_ANATOMY/FreeSeg_result/SUMA/KastnerNewModelSurfaces1D'

setwd( testSubjectDir )
rawFilenames <- list.files( path=getwd(), recursive = FALSE, full.names = FALSE )
testSubject <- 'AFH28'
cleanFilenames <- rep('a',length(rawFilenames))
for ( nModeName in 1:length(rawFilenames) ) {
  testFilename <- rawFilenames[ nModeName ]
  cleanFilenames[ nModeName ] <- strsplit( testFilename, testSubject )[[1]][2] 
}

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

#### create outputFolder folder, kastnerClassic ####
dirToCheck <- sprintf('%s/kastnerNewModelAverageModelsSurfaces1D', surfaceFolder )
if ( dir.exists( dirToCheck ) ) {
  instr <- sprintf( 'rm -R %s', dirToCheck ) 
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
}  
if ( !dir.exists( dirToCheck ) ) {
  instr <- sprintf('mkdir %s', dirToCheck)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
}         

for ( nAverages in 1:length( cleanFilenames ) ) { #length( cleanFilenames ) #nAverages <- 1
  
  dirTest <- list.files( path=anatomyDir, pattern=cleanFilenames[ nAverages ], all.files = TRUE, recursive = TRUE )

  setwd( anatomyDir )
  print( getwd() )
  print( dirTest[ 1 ] )
  fileTemp <- read.table( dirTest[ 1 ], as.is = TRUE )
  
  fileStore <- array( 0, dim( fileTemp ) )
  for ( nFiles in 1:length( dirTest ) ) { #nFiles <- 1
    setwd( anatomyDir )
    print( getwd() )
    print( dirTest[ nFiles ] )
    fileTemp <- read.table( dirTest[ nFiles ], as.is = TRUE )
    fileStore <- fileStore + fileTemp
  }
  fileStore <- fileStore/length( dirTest )
  fileNameOutput <- paste( 'average', cleanFilenames[ nAverages ], sep = '' )

  write.table( x=round(fileStore,4), file=sprintf('%s/kastnerNewModelAverageModelsSurfaces1D/%s', surfaceFolder, fileNameOutput ), col.names = FALSE, row.names = FALSE )
  
}




