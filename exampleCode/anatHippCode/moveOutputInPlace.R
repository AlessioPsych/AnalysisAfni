rm(list=ls())
mainDir <- '/analyse/Project0370'
codeDir <- '/analyse/Project0370/hippunfold_output/anatHippSimon/anatHippCode'
inputDir <- '/analyse/Project0370/hippunfold_output/anatHippSimon/anatHippData'
setwd( mainDir )
print( getwd() )

subjDirsFull <- dir( pattern='2023*' )
subjDirs <- subjDirsFull[ 3:length(subjDirsFull) ]
inputDirs <- dir( inputDir, pattern = '*_output')

data.frame( subjDirs, inputDirs )
runCodeFlag <- 1

for (i in 1:length(inputDirs)) { #i <- 1
  
  if ( dir.exists( sprintf( '%s/%s/ANATOMY/FreeSeg_result/SUMA/%s_hippUnfold', mainDir, subjDirs[i], inputDirs[i] ) )  ) {
    instr <- sprintf( 'rm -R %s/%s/ANATOMY/FreeSeg_result/SUMA/%s_hippUnfold', mainDir, subjDirs[i], inputDirs[i] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- sprintf( 'mkdir %s/%s/ANATOMY/FreeSeg_result/SUMA/%s_hippUnfold', mainDir, subjDirs[i], inputDirs[i] )  
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  setwd( inputDir )
  print( getwd() )
  instr <- sprintf('cp -R %s/hippunfold/sub-000/anat %s/%s/ANATOMY/FreeSeg_result/SUMA/%s_hippUnfold', inputDirs[i], mainDir, subjDirs[i], inputDirs[i] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  setwd( sprintf( '%s/%s/ANATOMY/FreeSeg_result/SUMA', mainDir, subjDirs[i] ) )
  print( getwd() )
  instr <- sprintf('cp anatCopy.nii.gz %s/%s/ANATOMY/FreeSeg_result/SUMA/%s_hippUnfold/anat/anatCopy.nii.gz', mainDir, subjDirs[i], inputDirs[i] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  setwd( sprintf( '%s/%s/ANATOMY/FreeSeg_result/SUMA', mainDir, subjDirs[i] ) )
  print( getwd() )
  instr <- sprintf('chmod 777 -R %s_hippUnfold/', inputDirs[i] )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
}
