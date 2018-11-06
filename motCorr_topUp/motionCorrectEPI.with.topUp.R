args <- commandArgs(T)
print( args )

#args <- c( 'Copy_of_EPI/', 'TopUp_attempt08/', '2')
#setwd('/home/jelle/Documents/Projects/contrast_transfer_function/JeDi_high_res_13122016_contrasts03')

# get actual dir
mainDir <- getwd()

# get EPIs
setwd( args[1] )
epiFiles <- dir( pattern=sprintf('*.nii') )

# get EPI n of TRs, first EPI
selectedEpi <- as.numeric(args[3])
if (is.numeric(selectedEpi)==FALSE) { # selected epi is not numeric!
  msg <- sprintf( 'the third argument MUST be a number (which EPI do you want to realign all the data to??)' )
  warning( msg )
  stopifnot(flagDir)  
}
setwd( mainDir )
instr <- sprintf('3dinfo %s%s > infoEPI.1D',args[1],epiFiles[selectedEpi])
system( instr )
trString <- scan( file='infoEPI.1D', what=c('character'), skip=17, nlines=1 )
nTRs <- as.numeric( trString[6] )
system('rm infoEPI.1D')

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
  instr <- sprintf('3dNwarpApply -master %s%s -source %s%s -nwarp motionCorrectEpi/%s %swarpTop_PLUS_WARP+orig -interp wsinc5 -prefix %s/%s', args[1], epiFiles[nEpi], args[1], epiFiles[nEpi], motionCorrAffine[nEpi], args[2], targetDir, filenameBlip)
  print( instr )
  system( instr )
}

