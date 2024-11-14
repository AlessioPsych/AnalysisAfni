rm(list=ls())

generalFolder <- '/analyse/Project0226/tests/Daniela_pilotData'
setwd( generalFolder )
if ( dir.exists( sprintf( '%s/data4Daniela', generalFolder ) ) ) {
  instr <- sprintf('rm -R %s/data4Daniela', generalFolder )
  print( instr )
  system( instr )
}
instr <- sprintf('mkdir %s/data4Daniela', generalFolder )
print( instr )
system( instr )

for (nSubj in 1:8) { #nSubj <- 1
  
  if ( nSubj==1 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject03082023_Daniela'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject03082023_Daniela/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==2 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Holly'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Holly/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==3 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Belinda'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Belinda/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==4 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_MTR13'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_MTR13/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==5 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_NWS30'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject28092023_NWS30/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==6 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject10102023_Rosie'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject10102023_Rosie/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==7 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_LucasEdward'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_LucasEdward/ANATOMY/anat_test/SUMA'
  }
  
  if ( nSubj==8 ) {
    mainDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_Theo'
    anatomyDir <- '/analyse/Project0226/tests/Daniela_pilotData/DanielaProject31102023_Theo/ANATOMY/anat_test/SUMA'
  }
  
  setwd( mainDir )
  getwd()
  
  subjFilename <- strsplit( mainDir, '[_]')[[1]][3]
  
  # create folders
  if ( !dir.exists( sprintf( '%s/data4Daniela/%s', generalFolder, subjFilename ) ) ) {
    instr <- sprintf( 'mkdir %s/data4Daniela/%s', generalFolder, subjFilename )
    print( instr )
    system( instr )
  }
  if ( !dir.exists( sprintf( '%s/data4Daniela/%s/motionCorrect_topUp_Epi', generalFolder, subjFilename ) ) ) {
    instr <- sprintf( 'mkdir %s/data4Daniela/%s/motionCorrect_topUp_Epi', generalFolder, subjFilename )
    print( instr )
    system( instr )
  }
  if ( !dir.exists( sprintf( '%s/data4Daniela/%s/motionEstimates', generalFolder, subjFilename ) ) ) {
    instr <- sprintf( 'mkdir %s/data4Daniela/%s/motionEstimates', generalFolder, subjFilename )
    print( instr )
    system( instr )
  }
  
  setwd( mainDir )
  setwd( 'motionCorrect_topUp_Epi' )    
  getwd()
  filesToconvert <- dir( pattern='*.BRIK')
  # instruction to convert files:
  for ( nFiles in 1:length( filesToconvert ) ) {
    filenameVolreg <- strsplit( filesToconvert[ nFiles ], '[.]' )
    filenameVolreg <- paste( filenameVolreg[[1]][1], '_', filenameVolreg[[1]][2], '_volreg', sep='' )
    instr <- sprintf( '3dcopy %s %s/data4Daniela/%s/motionCorrect_topUp_Epi/%s.nii.gz', filesToconvert[ nFiles ], generalFolder, subjFilename, filenameVolreg )
    print( instr )
    system( instr )
  }

  setwd( mainDir )
  setwd( 'motionCorrectEpi' )    
  getwd()
  filesToCopy <- dir( pattern='*.lin.volreg')
  # instruction to copy motion estimates:
  for ( nFiles in 1:length( filesToCopy ) ) {
    instr <- sprintf( 'cp %s %s/data4Daniela/%s/motionEstimates/%s.txt', filesToCopy[ nFiles ], generalFolder, subjFilename, filesToCopy[ nFiles ] )
    print( instr )
    system( instr )
  }
  
  #copy stim times
  setwd( mainDir )
  getwd()
  instr <- sprintf('cp -R stims/ %s/data4Daniela/%s', generalFolder, subjFilename )
  print( instr )
  system( instr )
  
  # copy anatomy
  setwd( mainDir )
  getwd()
  instr <- sprintf('cp ANATOMY/anatCopy.nii.gz %s/data4Daniela/%s/anatCopy.nii.gz', generalFolder, subjFilename )
  print( instr )
  system( instr )
  
}


  