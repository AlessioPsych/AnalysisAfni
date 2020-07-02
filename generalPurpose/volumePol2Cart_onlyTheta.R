args <- commandArgs(T)
print( args )

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )

inputVol <- args[1]
thetaIdx <- args[2]
labelOut <- args[3]

print( sprintf( 'input volume: %s', inputVol ) )
print( sprintf( 'theta idx (in afni): %s', thetaIdx ) )
print( sprintf( 'output volume: %s', labelOut ) )

volumePol2Cart <- function( inputVol, thetaIdx, labelOut  ) {
  instr <- sprintf( '3dcalc -a %s[%s] -expr \u0027 sin( a*PI/180 ) \u0027 -prefix del_x.nii.gz', inputVol, thetaIdx, inputVol ); system( instr )
  instr <- sprintf( '3dcalc -a %s[%s] -expr \u0027 cos( a*PI/180 ) \u0027 -prefix del_y.nii.gz', inputVol, thetaIdx, inputVol ); system( instr )
  instr <- sprintf( '3dTcat -prefix %s %s del_x.nii.gz del_y.nii.gz', labelOut, inputVol ); system( instr )
  instr <- sprintf( '3dinfo -nt %s > _ttt_nvols.1D', inputVol ); system( instr )
  nVols <- as.numeric( read.table( '_ttt_nvols.1D' ) )
  instr <- sprintf( '3drefit -sublabel %1.0f x_conv %s', nVols, labelOut ) ; system( instr ) 
  instr <- sprintf( '3drefit -sublabel %1.0f y_conv %s', (nVols+1), labelOut ) ; system( instr ) 
  system('rm _ttt_nvols.1D del_x.nii.gz del_y.nii.gz')
}


volumePol2Cart( inputVol, thetaIdx, labelOut  )
