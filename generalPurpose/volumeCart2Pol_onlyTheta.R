args <- commandArgs(T)
print( args )

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )

inputVol <- args[1]
xIdx <- args[2]
yIdx <- args[3]
labelOut <- args[4]

print( sprintf( 'input volume: %s', inputVol ) )
print( sprintf( 'x idx (in afni): %s', xIdx ) )
print( sprintf( 'y idx (in afni): %s', yIdx ) )
print( sprintf( 'output volume: %s', labelOut ) )

volumeCart2Pol <- function( inputVol, xIdx, yIdx, labelOut  ) {
  instr <- sprintf( '3dcalc -a %s[%s] -b %s[%s] -expr \u0027 mod( atan2(a,b)*180/PI + 360, 360 ) \u0027 -prefix del_ThetaAfterBlur.nii.gz', inputVol, xIdx, inputVol, yIdx ); system( instr )
  instr <- sprintf( '3dTcat -prefix %s %s del_ThetaAfterBlur.nii.gz', labelOut, inputVol ); system( instr )
  instr <- sprintf( '3dinfo -nt %s > _ttt_nvols.1D', inputVol ); system( instr )
  nVols <- as.numeric( read.table( '_ttt_nvols.1D' ) )
  instr <- sprintf( '3drefit -sublabel %1.0f theta_conv %s', nVols, labelOut ) ; system( instr ) 
  system('rm _ttt_nvols.1D del_ThetaAfterBlur.nii.gz')
}

volumeCart2Pol( inputVol, xIdx, yIdx, labelOut  )
