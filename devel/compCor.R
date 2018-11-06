args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/home/fracasso/data/hrfProject/AK_03021017_hrf/AK_03022017_hrf_WIP_pRFHRF_SENSE_09_1')
#args <- c('AK_03022017_hrf_WIP_pRFHRF_SENSE_09_1-t000-d0340.nii', '4', 'stimOn340.1D', '/packages/afni/16.1.13','/home/fracasso/analysisAfni/surfaces')

source( sprintf('%s/AFNIio.R', args[4] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[5] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[5] ) )
tsFilename <- args[1]

# preprocess
instr <- sprintf( '3dcopy %s _ttt_ts.nii.gz', tsFilename )
system( instr )
instr <- '3dAutomask -apply_prefix _ttt_ts_msk.nii.gz _ttt_ts.nii.gz'
system( instr )
instr <- sprintf( '3dDetrend -prefix _ttt_ts_msk_det.nii.gz -polort %s _ttt_ts_msk.nii.gz', args[2] )
system( instr )
instr <- '3dTstat -prefix _ttt_ts_msk_det_stdev.nii.gz -stdev _ttt_ts_msk_det.nii.gz'
system( instr )

# standard deviation & ts selection
stdevFile <- read.AFNI('_ttt_ts_msk_det_stdev.nii.gz')
stdevVolume <- stdevFile$brk[,,,1]
thrStdev <- quantile( array( stdevVolume[ stdevVolume>0 ] ), 0.98 )
voxelIndex <- which( stdevVolume>=thrStdev )

tsFile <- read.AFNI( '_ttt_ts_msk_det.nii.gz' )
extractTs <- function( niiFile, vIdx ) {
  coordsIdx <- coordinateFromLinearIndex( vIdx, dim( tsFile$brk )[ seq(1,3) ]  )
  tsOut <- tsFile$brk[ coordsIdx[1], coordsIdx[2], coordsIdx[3], ]
  return( ( tsOut-mean(tsOut) ) / sd(tsOut) )
}
tsSel <- sapply( voxelIndex, extractTs, niiFile=tsFile )

# ts exclusion due to significant voxels
stimSequence <- read.table( args[3] )
stimArray <- stimSequence[,1]
glmTs <- function( tsMat, stimulus, matIdx  ) {
  x <- stimulus
  y <- tsMat[,matIdx]
  mod <- summary( lm( y~x ) )
  return( mod$coefficients[2,4] )
}
arrayPval <- sapply( seq( 1, dim(tsSel)[2] ), glmTs, tsMat=tsSel, stimulus=stimArray  )

# save mask of selected voxels
emptyVol <- array( 0, dim( stdevVolume ) )
emptyVol[ voxelIndex[ arrayPval>=0.25 ] ] <- 1
volFileName <- sprintf( '_ttt_selMask.nii.gz' )
write.AFNI(volFileName,
           brk=emptyVol,
           label=NULL,
           view='+orig',
           orient=stdevFile$orient,
           origin=stdevFile$origin,
           defhead=stdevFile$NI_head )

#compute principal components
instr <- '3dpc -pcsave 6 -nscale -mask _ttt_selMask.nii.gz -prefix _ttt_ts_msk_det_pc _ttt_ts_msk_det.nii.gz'
system( instr )

pcMatrix <- as.matrix( read.table('_ttt_ts_msk_det_pc_vec.1D', as.is=TRUE) )
# meanPcMatrix <- apply( pcMatrix, 1, mean  )
# 
# ## plot
# instr <- '3dinfo -tr _ttt_ts_msk_det.nii.gz > _ttt_tr.1D'
# system( instr )
# trValue <- read.table( '_ttt_tr.1D' )
# 
# fs <- 1/trValue[1,1]
# N <-length( meanPcMatrix )
# index_t <- seq(1,N,1)
# f <- index_t*fs/N
# xf <- 2*abs(fft(meanPcMatrix ))/N
# x11(width=5, height=4)
# plotIndex <- index_t[ 1:floor(length(index_t)/2) ]
# xfPlot <- xf[ plotIndex  ]^2
# fPlot <- f[ plotIndex  ]
# plot( xfPlot~fPlot , type='l', bty='n', lwd=2, col='blue', xlab='frequency (hz)', ylab='power (dB)')


## clean up 
instr <- sprintf( '3dinfo -prefix_noext %s > _ttt_prefixNOEXT.1D', args[1] )
system( instr )
filenameBase <- scan( '_ttt_prefixNOEXT.1D', n=1, what=list('character') )[[1]]
#dev.copy2pdf( file=sprintf('%s_ps.pdf',filenameBase), width=5, height=4 )
write.table( round( pcMatrix, 6 ), file = sprintf( '%s_compCor.1D', filenameBase ), row.names = FALSE, col.names = FALSE )
instr <- sprintf( 'mv _ttt_selMask.nii.gz %s_mask_compCor.nii.gz', filenameBase )
system( instr )
system('rm _ttt_*')



