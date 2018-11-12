args <- commandArgs(T)
print( args )

#args <- c( 'EPI/', 'TOPUP/', '2', '2' )
#setwd('/media/alessiofracasso/storage2/SpinozaTest/HighRes/AF_HighRes_04112016/topUpDataset2')

# get actual dir
mainDir <- getwd()

# get EPIs
setwd( args[1] )
epiFiles <- dir( pattern=sprintf('*.nii') )

# get TOPUPs
setwd( mainDir )
setwd( args[2] )
topFiles <- dir( pattern=sprintf('*.nii') )


# get selected EPI 
selectedEpi <- as.numeric(args[3])
if (is.numeric(selectedEpi)==FALSE) { # selected epi is not numeric!
  msg <- sprintf( 'the third argument MUST be a number (which EPI do you want to realign all the data to??)' )
  warning( msg )
  stopifnot(flagDir)  
}

# get selected TOP UP 
selectedTop <- as.numeric(args[4])
if (is.numeric(selectedTop)==FALSE) { # selected TOP is not numeric!
  msg <- sprintf( 'the third argument MUST be a number (which TOP UP do you want to realign all the data to??)' )
  warning( msg )
  stopifnot(flagDir)  
}

# get EPI n of TRs
setwd( mainDir )
instr <- sprintf('3dinfo %s%s > infoEPI.1D',args[1],epiFiles[selectedEpi])
system( instr )
trString <- scan( file='infoEPI.1D', what=c('character'), skip=17, nlines=1 )
nTRs <- as.numeric( trString[6] )
system('rm infoEPI.1D')

# get TOPUP n of TRs
setwd( mainDir )
instr <- sprintf('3dinfo %s%s > infoTOPUP.1D',args[2],topFiles[selectedTop])
system( instr )
trString <- scan( file='infoTOPUP.1D', what=c('character'), skip=17, nlines=1 )
nTRsTOPUP <- as.numeric( trString[6] )
system('rm infoTOPUP.1D')


# motionCorrect EPIs
print('##################')
print('##################')
print('motionCorrect EPIs')
print('##################')
print('##################')
setwd( args[1] )
for ( nEpi in 1:length(epiFiles) ) {
  filename <- strsplit(epiFiles[nEpi],'.nii')[[1]][1]
  prefixName <- sprintf('pb.%s.volreg+orig', filename)
  motion1DfileAff <- sprintf('pb.%s.volreg', filename)  
  motion1DfileLin <- sprintf('pb.%s.lin.volreg', filename)  
  instr <- sprintf('3dvolreg -verbose -zpad 1 -base %s[%1.0f] -1Dfile %s -1Dmatrix_save %s -prefix %s -Fourier %s ', epiFiles[selectedEpi], nTRs-1, motion1DfileLin, motion1DfileAff, prefixName, epiFiles[nEpi] )
  print( instr )
  system( instr )
}

setwd( mainDir )
targetDir <- 'motionCorrectEpi'
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}

setwd( args[1] )
filesToMove <- dir(pattern='*volreg*')
setwd( mainDir )
for ( nFiles in 1:length(filesToMove) ) {
  instr <- sprintf('mv %s%s %s/%s', args[1], filesToMove[nFiles], targetDir, filesToMove[nFiles] )
  system( instr )
}


# EPI for Top Up
print('####################')
print('####################')
print('copy EPI for top up')
print('####################')
print('####################')
setwd( args[1] ) 
instr <- sprintf( '3dTcat -prefix ../epiForTopUp.nii %s[%1.0f]', epiFiles[selectedEpi], nTRs-1 )
system( instr )
setwd( mainDir )

# TOP UP for Top Up
print('##################################')
print('##################################')
print('copy top up volume EPIs for top up')
print('##################################')
print('##################################')
setwd( args[2] ) 
instr <- sprintf( '3dTcat -prefix ../topUp.nii %s[%1.0f]', topFiles[selectedTop], 0 )
system( instr )
setwd( mainDir ) 

# Compute non-linear transformation
print('#########################')
print('#########################')
print('compute top up transform')
print('#########################')
print('#########################')
setwd( mainDir )
print('estimate non-linear transformation, it might take a while...')
instr <- sprintf('3dQwarp -source epiForTopUp.nii -base topUp.nii -prefix warpTop -verb -iwarp -pblur 0.05 0.05 -blur -1 -1 -noweight -minpatch 9 -plusminus' )
print( instr )
system( instr )


# Apply non-linear transformation, create target dir
print('#######################')
print('#######################')
print('apply top up transform')
print('#######################')
print('#######################')
targetDir <- 'motionCorrect_topUp_Epi'
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}

# apply non-linear transformation to the files and save them in target dir
setwd('motionCorrectEpi')
motionCorrAffine <- dir(pattern='*aff12.1D')
setwd( mainDir )
for ( nEpi in 1:length(epiFiles) ) {
  namePrefix <- strsplit( epiFiles[nEpi], '.nii' )[[1]][1]
  filenameBlip <- sprintf( 'pb.%s.volreg+orig', namePrefix )
  instr <- sprintf('3dNwarpApply -master %s%s -source %s%s -nwarp motionCorrectEpi/%s warpTop_PLUS_WARP+orig -interp wsinc5 -prefix %s/%s', args[1], epiFiles[nEpi], args[1], epiFiles[nEpi], motionCorrAffine[nEpi], targetDir, filenameBlip)
  print( instr )
  system( instr )
}

# storage directory
print('###################')
print('###################')
print('clean up transform')
print('###################')
print('###################')
targetDir <- 'topUpDir'
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}

system('mv warpTop* topUpDir/')
system('mv epiForTopUp.nii topUpDir/')
system('mv topUp.nii topUpDir/')
