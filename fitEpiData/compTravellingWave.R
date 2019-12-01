args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/media/alessiof/Data/tests/phaseEncodedTest')
#args <- c('meanWedges_cat.nii.gz', '7', '1' )

AFNI_INSTALLDIR <- Sys.getenv( x = 'AFNI_INSTALLDIR' )
generalPurpose_DIR <- Sys.getenv( x = 'AFNI_TOOLBOXDIRGENERALPURPOSE' )
surfaces_Dir <- Sys.getenv( x = 'AFNI_TOOLBOXDIRSURFACES' )
source( sprintf('%s/AFNIio.R', AFNI_INSTALLDIR ) )
source( sprintf('%s/coordinateFromLinearIndex.r', surfaces_Dir ) )
source( sprintf('%s/linearIndexFromCoordinate.r', surfaces_Dir ) )
source( sprintf('%s/scaleData.R', generalPurpose_DIR ) )
library( pracma )
#ph_shift_input <- as.numeric( args[4] )
tsFilename <- args[1]
instr <- sprintf('3dcopy %s ttt_data.nii.gz',args[1])
system( instr )
instr <- '3dTstat -mean -prefix ttt_mean.nii.gz ttt_data.nii.gz'
system( instr )
instr <- '3dcalc -a ttt_data.nii.gz -b ttt_mean.nii.gz -expr \u0027min((a/b)*100,200)\u0027 -prefix ttt_percent.nii.gz'
system( instr )
instr <- sprintf( '3dDetrend -polort %s -prefix ttt_data_detrend.nii.gz ttt_percent.nii.gz', args[3] )
print(instr)
system( instr )

dataMean <- read.AFNI( filename='ttt_data.nii.gz')
dataBrkMean <- dataMean$brk 
dataFile <- read.AFNI( filename='ttt_data_detrend.nii.gz' )
dataBrk <- dataFile$brk
#dataFile <- read.AFNI( filename='ttt_data.nii.gz' )
#dataBrk <- dataFile$brk

progress <- 0.05
emptyVol <- array( 0, c( dim(dataBrk)[c(1:3)], 3 ) )
nCycles <- as.numeric( args[2] )
#selIdx <- seq( 1, ceiling( dim(dataBrk)[4]/2 ) )
selIdx <- seq( 1, 1+round( dim(dataBrk)[4]/2 ) )
emptyFFT <- array( 0, c( dim(dataBrk)[c(1:3)], length(selIdx) ) )
for ( k in 1:dim(dataBrk)[3] ) {
  
  dataLoop <- dataBrk[,,k,]
  dataMeanLoop <- dataBrkMean[,,k,]
  dimData <- dim( dataLoop )
  arraySlice <- array( dataLoop, c( prod( dimData[1:2] ), dimData[3] ) )
  arrayMeanSlice <- array( dataMeanLoop, c( prod( dimData[1:2] ), dimData[3] ) )
  
  #ft <- apply( arraySlice, 1, fft )
  #selIdx <- seq( 1, 1+round( dim(ft)[1]/2 ) )
  #ft <- ft[ selIdx , ];
  #scaledAmp <- abs( ft );
  #amp <- 2*( scaledAmp[nCycles+1, ] ) / dim(arraySlice)[1];
  #sqrtsummagsq <- sqrt( apply( scaledAmp^2, 2, sum ) )
  #coSlice <- scaledAmp[ nCycles+1, ] / sqrtsummagsq 
  #ph <- -( pi/2 ) - angle( ft[ nCycles+1, : ] );
  #ph(ph<0) = ph(ph<0)+pi*2;
  
  
  sliceFFT <- apply( arraySlice, 1, fft )
  sliceFFT <- t(sliceFFT)
  sliceMean <- apply( arrayMeanSlice, 1, mean )
  sliceMean <- t(sliceMean)
  sliceFFT_h <- sliceFFT[ , selIdx ] 
  
  scaledSliceFFT_h <- abs( sliceFFT_h )
  amp <- 2*( scaledSliceFFT_h[ ,nCycles+1 ] ) / dim(dataBrk)[4]
  sqrtsummagsq = sqrt( apply( scaledSliceFFT_h^2, 1, sum ) )
  coSlice <- scaledSliceFFT_h[ ,nCycles+1 ] / sqrtsummagsq 
  
# Calculate phase:
# 1) add pi/2 so that it is in sine phase.
# 2) minus sign because sin(x-phi) is shifted to the right by phi.
# 3) Add 2pi to any negative values so phases increase from 0 to 2pi.

  #ph_slice <- -(pi/2) - atan2( Im( sliceFFT[,nCycles+1] ), Re( sliceFFT[,nCycles+1] ) ) 
  ph_slice <- -(pi/2) - angle( sliceFFT[,nCycles+1] )
  ph_slice[ph_slice<0] = ph_slice[ph_slice<0]+pi*2;
  #ph_slice <- ( ph_slice + ph_shift_input ) %% pi*2
  
  #emptyVol[,,k,1] <- round( amp/sliceMean*100, 4)
  emptyVol[,,k,1] <- round( amp, 4)
  emptyVol[,,k,2] <- round( coSlice, 4)
  emptyVol[,,k,3] <- round( ph_slice, 4)
  
  periodogram <- 2*( abs( sliceFFT ) / dim(dataBrk)[4] )
  emptyFFT[,,k,] <- round( periodogram[,selIdx], 4 )                          
  
  if (k>dim(dataBrk)[3]*progress) {
    cat(paste('*',""))
    progress <- progress + 0.05
  }
  
}

print('prepare and save output...')

dataFileOut <- read.AFNI( filename='ttt_mean.nii.gz')
emptyOut <- array(0,dim(dataFileOut$brk))
emptyOut[,,,1] <- emptyVol[,,,1]
volFileName <- sprintf( 'ttt_out_amp.nii.gz' )
write.AFNI(volFileName,
           brk=emptyOut,
           label=c('amp'),
           view='+orig',
           orient=dataFileOut$orient,
           origin=dataFileOut$origin,
           defhead=dataFileOut$NI_head )

emptyOut[,,,1] <- emptyVol[,,,2]
volFileName <- sprintf( 'ttt_out_co.nii.gz' )
write.AFNI(volFileName,
           brk=emptyOut,
           label=c('co'),
           view='+orig',
           orient=dataFileOut$orient,
           origin=dataFileOut$origin,
           defhead=dataFileOut$NI_head )

emptyOut[,,,1] <- emptyVol[,,,3]
volFileName <- sprintf( 'ttt_out_ph.nii.gz' )
write.AFNI(volFileName,
           brk=emptyOut,
           label=c('ph'),
           view='+orig',
           orient=dataFileOut$orient,
           origin=dataFileOut$origin,
           defhead=dataFileOut$NI_head )

volFileName <- sprintf( 'periodogram.nii.gz' )
write.AFNI(volFileName,
           brk=emptyFFT,
           label=NULL,
           view='+orig',
           orient=dataFileOut$orient,
           origin=dataFileOut$origin,
           defhead=dataFileOut$NI_head )


instr <- '3dTcat -prefix trWave.nii.gz ttt_out_amp.nii.gz ttt_out_co.nii.gz ttt_out_ph.nii.gz'
system( instr )
system('applyObliqueMatrix.sh trWave.nii.gz ttt_mean.nii.gz')
system('rm trWave.nii.gz')
system('mv trWave_refit.nii.gz trWave.nii.gz')
system('applyObliqueMatrix.sh periodogram.nii.gz ttt_mean.nii.gz')
system('rm periodogram.nii.gz')
system('mv periodogram_refit.nii.gz periodogram.nii.gz')
system('rm ttt_*')

# instr <- '3dPeriodogram -prefix ttt_data_periodogram.nii.gz ttt_data.nii.gz'
# system( instr )
# 
# selVol <- strsplit( args[2], '[-]' )[[1]]
# selVolCorr <- as.character( as.numeric(selVol)-1 )
# selVolChr <- paste( selVolCorr, collapse=',' )
# inputChr <- sprintf('ttt_data_periodogram.nii.gz[%s]',selVolChr)
# 
# instr <- '3dinfo -nv ttt_data.nii.gz > ttt_nsubbrik.1D'
# system( instr )
# nsubbrik <- read.table( 'ttt_nsubbrik.1D' )
# ampChr <- sprintf( '2*abs(a)/%s', as.character(nsubbrik) ) 
# instr <- sprintf('3dcalc -a %s -expr \u0027%s\u0027 -prefix ttt_data_periodogram_scale.nii.gz', inputChr, ampChr)
# system( instr )
# instr <- '3dcalc -a ttt_data_periodogram.nii.gz -expr \u0027abs(a)*abs(a)\u0027 -prefix ttt_data_periodogram_squared.nii.gz'
# system( instr )
# instr <- '3dTstat -prefix ttt_data_periodogram_squared_sum.nii.gz -sum ttt_data_periodogram_squared.nii.gz'
# system( instr )
# instr <- '3dcalc -a ttt_data_periodogram_squared_sum.nii.gz -expr \u0027sqrt(a)\u0027 -prefix ttt_data_periodogram_squared_sum_sqrt.nii.gz'
# system( instr )
# instr <- sprintf('3dcalc -a %s -b ttt_data_periodogram_squared_sum_sqrt.nii.gz -expr \u0027abs(a)/b\u0027 -prefix ttt_data_co.nii.gz', inputChr)
# system( instr )
# instr <- sprintf('3d')
# 
# system('rm ttt_*')
# 
# 
# ph = -(pi/2) - angle(ft(nCycles+1,:));
# ph(ph<0) = ph(ph<0)+pi*2;

#amp = 2*(scaledAmp(nCycles+1,:))/size(ptSeries,1);
#sqrtsummagsq = sqrt(sum(scaledAmp(noiseIndices, :).^2));

