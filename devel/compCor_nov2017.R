
args <- commandArgs(T)
args

#rm(list=ls())
#setwd('/home/fracasso/data/compCorTest')
#args <- c('EPI/','0.5','EPI_noDetrend_nocompCor/','*.nii','0','stimTimingCompCor/',0)

#load libraries
library(ANTsR)
library(neuRosim)
library(pracma)
source( sprintf('%s/generalPurpose/powerSpectra.R', args[8] ) )

#define variables
mainDir <- getwd()
epiDir <- args[1]
clfracPar <- args[2]
epiWriteDir <- args[3]
filesToGet <- args[4]
polortParameter <- args[5]
stimulusDir <- args[6]
justDetrendFlag <- args[7]

#get input epi files
setwd( mainDir )
setwd( epiDir )
filesInput <- dir( pattern=sprintf( '%s', filesToGet ) )

setwd(mainDir)
system( sprintf('mkdir %s', epiWriteDir) )

if (stimulusDir!=0) {
  #define the function to identify the the TRs where a stimuli was presented
  findNumbers <- function( idxIn, targetVector, sampleVector ) {
    tol <- 0.001
    return( which( abs(targetVector-sampleVector[idxIn]) < tol ) )
  }
  #define the function to compute the correlation between the expected TS and the compCor PCA components
  testCor <- function( idx, xInt, yInt  ) {
    cOut <- cor.test( xInt[,idx], yInt )
    return( cOut$p.value )
  }
}

for ( nFiles in 1:length(filesInput) ) { #loop through EPI files
  
  # get the file, mask it and compute the components (with ANTsR)
  setwd( mainDir )
  setwd( epiDir )
  instr <- sprintf('3dTstat -prefix _ttt_meanEpi.nii.gz %s', filesInput[ nFiles ] )
  system( instr )
  instr <- sprintf('3dAutomask -clfrac %s -prefix _ttt_maskEpi.nii.gz %s',  clfracPar, '_ttt_meanEpi.nii.gz' )
  system( instr )
  maskEpi <- antsImageRead( '_ttt_maskEpi.nii.gz' )
  inputEPI <- antsImageRead( filesInput[ nFiles ] )
  compMat <- compcor( inputEPI, mask=maskEpi, variance_extreme=0.99  )

  if (stimulusDir!=0) {
    #get the 1D file with the timings for the stimuli
    setwd( mainDir )
    setwd( stimulusDir )
    stimFile <- dir(pattern = '.1D')
    dataTime <- as.numeric( read.table( stimFile[nFiles], as.is=TRUE ) )
  }
  
  # get the TR of the EPI sequence and compute the expected TS plus its power spectrum
  setwd( mainDir )
  setwd( epiDir )
  instr <- sprintf('3dinfo -TR %s > _ttt_tr.1D', filesInput[ nFiles ] )
  system( instr )
  TR <- as.numeric( read.table('_ttt_tr.1D', as.is=TRUE) )
  system('rm _ttt_tr.1D')

  x11(width=12, height=5)
  par(mfrow=c(2,3))
  if (stimulusDir!=0) {
    tsStepX <- seq( TR, TR*dim(compMat)[1], length.out = dim(compMat)[1] )
    tsStepY <- rep(0,length(tsStepX))
    idxTime <- sapply( 1:length(dataTime), findNumbers, targetVector=tsStepX, sampleVector=dataTime  )
    tsStepY[ idxTime ] <- 1
    hrf <- canonicalHRF( seq(1,25,by=TR) )
    outTsTemp <- conv( tsStepY, hrf )
    outTs <- outTsTemp[1:length(tsStepX)]
    expOut <- powerSpectra( outTs, TR )
  
    #compute the p value of the correlation between the expected TS and the estimated components
    corPvals <- sapply( 1:4, testCor, xInt=compMat, yInt=outTs  )
    #plot the expected TS plus its power spectrum
    plot( outTs~tsStepX, type='l', xlab='time (s)', ylab='activity (A.U.)', bty='n', main='expected TS' )
    plot( expOut$amp~expOut$freq, type='l', xlab='frequency (hz)', ylab='amplitude', bty='n', main='power spectra, expected TS' )
  }
  if (stimulusDir==0) {
    corPvals <- rep(999,4)
  }
  
  for (nPlot in 1:4) {
    psOut <- powerSpectra( compMat[,nPlot], TR )
    plot( psOut$amp ~ psOut$freq, type='l',
          xlab='frequency (hz)', ylab='amplitude', bty='n',
          main=sprintf('power spectra, compcor component \n %1.2f',round(corPvals[nPlot],3))  )
  }
  dev.copy2pdf( file=sprintf( '../%s%s_fileReport.pdf', epiWriteDir, strsplit( filesInput[ nFiles ], '[.]' )[[1]][1] ), width=12, height=5 )
  
  #remove the components / detrend, depending on the correlation / input parameters
  if (justDetrendFlag==0) {
    if (length(which( corPvals >= 0.05 ))==0) { #==0 #if all components are correlated with expected TS, do not remove anything, just detrend if requested
      outName <- sprintf( '%s_corr_compCor.nii', strsplit( filesInput[ nFiles ], '[.]' )[[1]][1] )
      instr <- sprintf('3dDetrend -prefix _ttt_compCor.nii -polort %s %s', polortParameter, filesInput[ nFiles ] )
      system( instr )
      instr <- sprintf('3dcalc -a _ttt_compCor.nii -b _ttt_meanEpi.nii.gz -expr \u0027a+b\u0027 -prefix %s', outName )
      system( instr )
    }
    if (length(which( corPvals >= 0.05 ))>0) { #>0 #remove only those components that are not correlated with expected TS
      selectedComponentsIdx <- which( corPvals >= 0.05 )
      selectedComponents <- round( compMat[,selectedComponentsIdx], 4 )
      write.table( selectedComponents, file='_ttt_pca.1D', row.names = FALSE, col.names = FALSE )
      outName <- sprintf( '%s_corr_compCor.nii', strsplit( filesInput[ nFiles ], '[.]' )[[1]][1] )
      instr <- sprintf('3dDetrend -prefix _ttt_compCor.nii -polort %s -vector _ttt_pca.1D %s', polortParameter, filesInput[ nFiles ] )
      system( instr )
      instr <- sprintf('3dcalc -a _ttt_compCor.nii -b _ttt_meanEpi.nii.gz -expr \u0027a+b\u0027 -prefix %s', outName )
      system( instr )
    }
  }
  if (justDetrendFlag==1) { #if only detrend is requested, then just detrend data according to the parameter at input
    outName <- sprintf( '%s_corr_compCor.nii', strsplit( filesInput[ nFiles ], '[.]' )[[1]][1] )
    instr <- sprintf('3dDetrend -prefix _ttt_compCor.nii -polort %s %s', polortParameter, filesInput[ nFiles ] )
    system( instr )
    instr <- sprintf('3dcalc -a _ttt_compCor.nii -b _ttt_meanEpi.nii.gz -expr \u0027a+b\u0027 -prefix %s', outName )
    system( instr )
  }
    
  system('rm _ttt*') #clean up

}

# move files to the designated folder (defined at the input)
setwd(mainDir)
system( sprintf('mv %s*_corr_compCor* %s', epiDir, epiWriteDir ) )
#system( sprintf('mv Rplots.pdf %scompCor_report.pdf', epiDir, epiWriteDir ) )

