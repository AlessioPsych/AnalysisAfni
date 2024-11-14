
debug <- 0

if ( debug==0 ) {
  args <- commandArgs(T)
  print( args )
  mainDir <- getwd()
  setwd( mainDir )
}
if ( debug==1 ) {
  rm( list=ls() ); gc();
  mainDir <- '/analyse/Project0226/tests/nilsTest/data/20230622_SHI27' 
  setwd( mainDir )
  args <- c('stim_data/feat_out.txt', 'ppp_EPI_predictions', '300','2')
}

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( neuRosim )
library( parallel )

# arrange inputs
inputFile <- args[1]
outDirectory <- args[2]
nDynamics <- as.numeric( args[3] )
trValue <- as.numeric( args[4] )

print( sprintf('input file = %s', inputFile ) )
print( sprintf('output directory = %s', outDirectory ) )
print( sprintf('n. trs = %1.0f', nDynamics ) )
print( sprintf('tr = %1.0f', trValue ) )

setwd( mainDir )
print( getwd() )

# create output directory
system( sprintf('mkdir %s', outDirectory) )

# load stimulus data
stimData <- read.table( inputFile, as.is=TRUE, header=FALSE )

#this part of the code builds a matrix with all the possible prediction tested
hrfDelayOnsetElements <- 3
hrfDelayUnderShootElements <- 3
hrfDispersion1Elements <- 3
hrfDispersion2Elements <- 3
hrfCElements <- 3

print('build prediction...')
par_hrf_a1 <- seq( 6, 9, length.out=hrfDelayOnsetElements )
par_hrf_a2 <- seq( 12, 15, length.out=hrfDelayUnderShootElements )
par_hrf_b1 <- seq( 0.8, 1.2, length.out=hrfDispersion1Elements )
par_hrf_b2 <- seq( 0.8, 1.2, length.out=hrfDispersion2Elements )
par_hrf_c <- seq( 0.3, 0.5, length.out=hrfCElements )

# generate the global prediction grid
predictionGridGlobal <- expand.grid( par_hrf_a1, par_hrf_a2, par_hrf_b1, par_hrf_b2, par_hrf_c )
print( sprintf('total number of predictions to generate: %1.0f', dim( predictionGridGlobal )[1] ) )

# this part builds a matrix with the starting and ending prediction index to fit for each loop in the next for loop
if (dim( predictionGridGlobal )[1] <= 900 ) {
  totalIterations <- 3
}
if (dim( predictionGridGlobal )[1] > 900 ) {
  totalIterations <- ceil( dim( predictionGridGlobal )[1]/500 )
}

limitsPrediction <- round( seq(1,dim(predictionGridGlobal)[1], length.out=totalIterations ) ) 
limitsPrediction[1] <- 1
limitsPrediction[ length(limitsPrediction) ] <- dim(predictionGridGlobal)[1]
limitsPredictionMatrix <- array( 0, c( (length(limitsPrediction)-1) , 2 ) )
limitsPredictionMatrix[,1] <- limitsPrediction[1:(length(limitsPrediction)-1)]
limitsPredictionMatrix[,2] <- limitsPrediction[2:(length(limitsPrediction))]
limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] <- limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
runIndexPredictions <- seq( 1:dim(limitsPredictionMatrix)[1] )
print( sprintf('total number of iterations: %1.0f', length( runIndexPredictions ) ) )

# MRI time and stimuli time
mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) # + as.numeric( trValue[1] )/2
samplingTime <- trValue

#### function to generate the predictions ####
generatePrediction <- function( indexPrediction, inputPredictionGrid, samplingTime, mriTime, scaleData, stimData ) {
  
  library(neuRosim)
  library(pracma)

  prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
  
  testPlot <- 0
  if (testPlot==1) {
    prfPar <- c( 6, 12, 0.9, 0.9, 0.35)
  }
  
  hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[1], a2=prfPar[2], b1=prfPar[3], b2=prfPar[4], c=prfPar[5] ), verbose=FALSE ); #plot( hrf ~ seq(0,30, samplingTime ) )
  
  stimSequence <- rep(0,length(mriTime))
  for (nStimuli in 1:dim(stimData)[1] ) { #nStimuli <- 1
    idxStimuliOn <- which( stimData[nStimuli,1]==mriTime )[1]
    stimuliLength <- stimData[nStimuli,2]
    stimSequence[ idxStimuliOn:(idxStimuliOn+stimuliLength-1) ] <- 1
  }
  
  if ( sum( is.na( stimSequence ) ) == 0 ) {
    
    pConv01 <- conv( stimSequence, hrf )
    pConvTrim01 <- pConv01[ 1 : length( mriTime ) ]; 
    returnPrediction <- round( scaleData( pConvTrim01, 1, 0 ), 6 ); 

  }
  
  if (testPlot==1) {
    par(mfrow=c(2,3))
    plot( hrf ~ seq(0,30, samplingTime ), type='l' )
    plot( stimSequence ~ mriTime, type='l' )
    plot( stimSequence ~ mriTime, type='l' )
    plot( pConvTrim01 ~ mriTime,type='l')
    plot( returnPrediction ~ mriTime,type='l' )
  }
  
  if ( sum( is.na( stimSequence ) ) > 0 ) {
    returnPrediction <- rep(0,length(mriTime))
  }
  
  return( returnPrediction ) 
  
} 

predictionWrapper <- function( passIdx, limitsMat, inputGrid ) {
  
  print( sprintf('iteration: %1.0f of %1.0f', passIdx, dim( limitsMat )[1] ) )
  selectedRows <- limitsMat[passIdx, ]
  
  # generate predictions:
  nCores <- 6
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( tsPredictions <- parSapply( cl, selectedRows[1]:selectedRows[2], generatePrediction,
                                                          inputPredictionGrid=inputGrid,
                                                          samplingTime=samplingTime,
                                                          mriTime=mriTime,
                                                          scaleData=scaleData,
                                                          stimData=stimData ) )
  stopCluster(cl)
  showConnections()
  rm(cl); gc();
  print( storeTimePar )
  #storeTimePar <- system.time( tsPredictions <- sapply( selectedRows[1]:selectedRows[2], generatePrediction, inputPredictionGrid=inputGrid ) )
  #print( storeTimePar )
  
  # concatenate parameters and predictions:
  outMatrix <- rbind( t( predictionGridGlobal[ selectedRows[1]:selectedRows[2], ] ), tsPredictions )
  
  # save parameters and predictions:
  save( outMatrix, file = sprintf('%s/ttt_predictionFile_%1.0f.RData', outDirectory, passIdx ) )
  
  return(1)
}

sapply( 1:dim(limitsPredictionMatrix)[1], predictionWrapper, limitsMat = limitsPredictionMatrix, inputGrid = predictionGridGlobal)

rm( list=ls() ); gc()

