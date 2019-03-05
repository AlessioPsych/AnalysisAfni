args <- commandArgs(T)
print( args )

## prepare the stimuli with the portrait as well

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')


#args <- c('maxVarEye.nii.gz', 'meanTs_eye_topUp_res.nii', 'output_del', '0', '3','0.166')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir ) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
source( sprintf('%s/vonMisesDist.R', generalPurposeDir ) )

## test vonMises distribution ##
#x <- seq(-pi,pi,0.1)
#mu <- pi*0.9
#k <- 1
#vonValues <- vonMisesDist( x, mu, k, 2, 0 )
#plot( vonValues~x )

library( pracma )
library( abind )
library( neuRosim )
library( parallel )
library( matrixcalc )

epiFile <- args[1]
inputFile <- args[2]
outSuffix <- args[3]
flagSurround <- as.numeric( args[4] )
polortArg <- as.numeric( args[5] )
samplingTime <- as.numeric( args[6] )

# load the data
# get orientation
print('detrend and load data...')
instr <- sprintf('3dinfo -orient %s > __tttt.1D', inputFile)
system( instr )
orientValue <- read.table('__tttt.1D')
system('rm __tttt.1D')
# get TR
instr <- sprintf('3dinfo -tr %s > __tttt.1D', inputFile)
system( instr )
trValue <- read.table('__tttt.1D')
system('rm __tttt.1D')

instr <- sprintf( '3dresample -prefix _ttt_orient.nii.gz -inset %s', inputFile )
system( instr)
instr <- sprintf('3dDetrend -polort %1.0f -prefix _ttt_data.nii.gz _ttt_orient.nii.gz', polortArg)
system( instr )
instr <- sprintf('3dresample -orient %s -prefix %s_detrendedData.nii.gz -inset _ttt_data.nii.gz', orientValue[1,1], outSuffix )
system( instr )
tsOriginal <- read.AFNI( '_ttt_orient.nii.gz' )
tsOriginalArray <- array( tsOriginal$brk, c( prod( dim( tsOriginal$brk )[1:3] ), dim(tsOriginal$brk)[4] ) )
ts <- read.AFNI( '_ttt_data.nii.gz' )
tsArray <- array( ts$brk, c( prod( dim( ts$brk )[1:3] ), dim(ts$brk)[4] ) )
system( 'rm _ttt_orient.nii.gz' )
system( 'rm _ttt_data.nii.gz' )

instr <- sprintf( '3dresample -prefix _ttt_mask.nii.gz -inset %s', epiFile )
system( instr)
meanEpi <- read.AFNI( '_ttt_mask.nii.gz' )
system( 'rm _ttt_mask.nii.gz' )


# load stimuli definition
print('get stimuli...')
setwd('/analyse/Project0226/dataSummary')
arrayStim <- scan( 'eyeMovingStim.txt' )
stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
stimMatFlip_eyeMoving <- aperm( stimMat[ dim(stimMat)[1]:1,, ], c(2,1,3) )

arrayStim <- scan( 'eyeFixStim_border.txt' )
stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
stimMatFlip_eyeFix <- aperm( stimMat[ dim(stimMat)[1]:1,, ], c(2,1,3) )

# load saccade direction definition
directionFile <- read.table( file='directionFile.1D', header = FALSE, as.is=TRUE )
directionFile[,3] <- directionFile[,3]*10000
directionCart <- pol2cart( cbind( deg2rad( directionFile[,1] ), directionFile[,2] ) )
directionPol <- cart2pol( directionCart )
directionPol[,1] <- round( directionPol[,1], 3 )
data.frame( directionFile, directionPol[,1] ) 
stimuliSequenceTr <- 200 #in milliseconds
sampleRepetitions <- directionFile[,3] / stimuliSequenceTr
baselines <- ifelse( sampleRepetitions>10, 1, 0 )
directionArray <- directionPol[,1]
directionArray_expanded_radiants <- rep( directionArray, sampleRepetitions )
saccadeOnsetBin <- which( directionArray_expanded_radiants[1:length(directionArray_expanded_radiants)-1] !=
                            directionArray_expanded_radiants[2:length(directionArray_expanded_radiants)] ) + 1
for (k in 1:length(saccadeOnsetBin)) {
  if (k==1) {
    saccStart <- seq( saccadeOnsetBin[k], saccadeOnsetBin[k]+(1600/stimuliSequenceTr)*7, by=(1600/stimuliSequenceTr) ) 
  }
  if (k>1) {
    tempArray <- seq( saccadeOnsetBin[k], saccadeOnsetBin[k]+(1600/stimuliSequenceTr)*7, by=(1600/stimuliSequenceTr) ) 
    saccStart <- c( saccStart, tempArray )
  }
}

setwd(mainDir)

## visualize stimuli
x11( width=3, height=3 )
for ( snap in 1:dim(stimMatFlip_eyeFix)[3] ) {
 image( stimMatFlip_eyeFix[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01)
}
x11( width=3, height=3 )
for ( snap in 1:dim(stimMatFlip_eyeMoving)[3] ) {
  image( stimMatFlip_eyeMoving[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01)
}
stimSeq_eyeFix <- stimMatFlip_eyeFix
stimSeq_eyeMoving <- stimMatFlip_eyeMoving
x <- seq(-10,10,length.out = dim(stimSeq_eyeFix)[1] )
y <- seq(-5.5,5.5,length.out = dim(stimSeq_eyeFix)[2] )
outMesh <- meshgrid(y, x)
outMesh$X <- scaleData( outMesh$X, 1, 0 )
outMesh$Y <- scaleData( outMesh$Y, 1, 0 )

xElements <- 6
yElements <- 6
sigmaArrayPositiveElements <- 4
hrfDelayOnsetElements <- 1
hrfDelayUnderShootElements <- 1
multParElements <- 2
centerVonElements <- 4
kVonElements <- 4
xElementsGain <- 4
yElementsGain <- 4
sigmaArrayPositiveElementsGain <- 4

print('build prediction...')
xPosFit <- seq( -9, 9, length.out=xElements )
yPosFit <- seq( -4.5, 4.5, length.out=yElements )
sigmaArrayPositive <- seq( 0.25, 6, length.out=sigmaArrayPositiveElements )
if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
par_hrf_a1 <- seq( 6, 10, length.out=hrfDelayOnsetElements )
par_hrf_a2 <- seq( 12, 16, length.out=hrfDelayUnderShootElements )
if (flagSurround==1) { multPar <- seq(0,0.8, length.out = multParElements) }
if (flagSurround==0) { multPar <- 0 }
centerVon <- seq( -0.95*pi, 0.95*pi, length.out=centerVonElements ) 
kVon <- seq( 0.8, 2.5, length.out=kVonElements ) 
xPosFit_gain <- seq( -9, 9, length.out=xElementsGain )
yPosFit_gain <- seq( -4.5, 4.5, length.out=yElementsGain )
sigmaArrayPositive_gain <- seq( 3, 6, length.out=sigmaArrayPositiveElementsGain )

predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, centerVon, kVon, xPosFit_gain, yPosFit_gain, sigmaArrayPositive_gain )
if (flagSurround==1) {
  keepPredictionIdx <- (predictionGridTemp[ ,3] < predictionGridTemp[ ,4]) & #center smaller than surround
    (predictionGridTemp[ ,4] < predictionGridTemp[ ,12]) #surround smaller than gain 
  predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]
}
if (flagSurround==0) {
  keepPredictionIdx <- (predictionGridTemp[ ,3] < predictionGridTemp[ ,12]) #center smaller than gain 
  predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]
}

dim( predictionGridGlobal )

# this part builds a matrix with the starting and ending prediction index to fit for each loop in the next for loop
if (dim( predictionGridGlobal )[1] <= 250 ) {
  totalIterations <- 3
}
if (dim( predictionGridGlobal )[1] > 250 ) {
  totalIterations <- ceil( dim( predictionGridGlobal )[1]/100 )
}
limitsPrediction <- round( seq(1,dim(predictionGridGlobal)[1], length.out=totalIterations) ) 
limitsPrediction[1] <- 1
limitsPrediction[ length(limitsPrediction) ] <- dim(predictionGridGlobal)[1]
limitsPredictionMatrix <- array( 0, c( (length(limitsPrediction)-1) , 2 ) )
limitsPredictionMatrix[,1] <- limitsPrediction[1:(length(limitsPrediction)-1)]
limitsPredictionMatrix[,2] <- limitsPrediction[2:(length(limitsPrediction))]
limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] <- limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
runIndexPredictions <- seq( 1:dim(limitsPredictionMatrix)[1] )

# select time series to fit based on the provided mask
indexVol <- meanEpi$brk[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
selIdxVoxel <- which( indexArray == 1 )
tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series

#modelFitCounter <- 1
for (modelFitCounter in 1:length(runIndexPredictions) ) { #length(runIndexPredictions)
  
  print( sprintf( 'iteration %1.0f of %1.0f, start...', modelFitCounter, length(runIndexPredictions)  ) )
  
  # here I select the portion of all the predictions that are going to be tested later
  predictionGrid <- predictionGridGlobal[ limitsPredictionMatrix[modelFitCounter,1]:limitsPredictionMatrix[modelFitCounter,2], ]
  
  print( sprintf( 'generate predictions in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  
  #### eye predictions ####
  eyeDirectionPrediction <- function(indexPrediction, inputPredictionGrid) {
    
    #indexPrediction <- 10
    #inputPredictionGrid <- predictionGrid

    prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
    muRad <- prfPar[8]
    kappaPar <- prfPar[9]
    stimuliSequenceTr <- 200 #in milliseconds, sample the sequence of direction changes every 200 ms
    sampleRepetitions <- directionFile[,3] / stimuliSequenceTr
    baselines <- ifelse( sampleRepetitions>10, 1, 0 ) #which are the baselines in the array?
    directionArray <- directionPol[,1]
    directionArray_expanded_radiants <- rep( directionArray, sampleRepetitions ) 
    baselines_expanded <- rep( baselines, sampleRepetitions )
    baselines_expanded <- rep(1,length(baselines_expanded))
    baselines_expanded[ saccStart ] <- 0 #code saccades with 0, baselines with 1, every 200 ms (stimuliSequenceTr)
    
    firstStepTs_orig <- suppressWarnings( vonMisesDist( directionArray_expanded_radiants, muRad, kappaPar, 0  ) ) # direction read throught the current VonMises distribution
    formattedBaselines <- abs(baselines_expanded-1) #code saccades as ones (1), all the rest as zero (0), every 200 ms (stimuliSequenceTr)
    formattedBaselines[ formattedBaselines==0 ] <- -1
    firstStepTs_orig_baseline <- firstStepTs_orig
    firstStepTs_orig_baseline[ formattedBaselines==-1 ] <- 0 #code saccades according to the vonMises distribution, every 200 ms (stimuliSequenceTr)
    
    convTs <- conv( firstStepTs_orig_baseline, hrf )
    convTsTrim <- convTs[1:length(firstStepTs_orig_baseline)]
    timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr/1000, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
    iTs <- interp1( timeArray, convTsTrim, seq(0,max(timeArray),length.out = dim(ts$brk)[4] ) )
    scaleConvTsTrim <- scaleData( round( convTsTrim, 5 ), 1, 0 )
    scaleIts <- scaleData( round( iTs, 5 ), 1, 0 )
    return( scaleIts )
    #iTs <- scaleData( round( iTs, 5 ), 1, 0 )

    # #to plot
    # par(mfrow=c(3,3))
    # piSpace <- seq(-pi*0.95, +pi*0.95, length.out = 200)
    # plot( formattedBaselines, type='l' )
    # plot( directionArray_expanded_radiants, type='l' )
    # plot( formattedBaselines * directionArray_expanded_radiants, type='l' )
    # plot( suppressWarnings( vonMisesDist( piSpace, muRad, kappaPar, 0  ) ) ~ piSpace )
    # firstStepTs_orig_plot_baseline <- firstStepTs_orig
    # firstStepTs_orig_plot_baseline[ formattedBaselines==-1 ] <- -4
    # plot( firstStepTs_orig_plot_baseline, type='l' )
    # plot( firstStepTs_orig, type='l' )
    # plot( firstStepTs_orig_baseline, type='l' ) # this is the line that matters, directions coded as the VonMises distribution, 0 everywhere else
    # plot( convTsTrim~timeArray )
    # lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), iTs, col='red' )
    # plot( scaleConvTsTrim~timeArray )
    # lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), scaleIts, col='red' )
    
  }
  eyePredictions_direction <- sapply( 1:dim(predictionGrid)[1], eyeDirectionPrediction, inputPredictionGrid=predictionGrid )
  
  #### prepare to generate the predictions, eye fix ####
  stimSeqMat <- array( stimSeq_eyeFix, c( length(x)*length(y), dim(stimSeq_eyeFix)[3] ) )
  incrementalCounter <- dim(predictionGrid)[1] * 0.05
  timeStimuli <- seq( 0, dim(stimSeq_eyeFix)[3]*samplingTime, length.out = dim(stimSeq_eyeFix)[3] )
  mriTime <- seq( 0, dim(stimSeq_eyeFix)[3]*samplingTime, length.out = dim(ts$brk)[4] )
  
  #### function to generate the predictions, eye fix ####
  generatePrediction_visual <- function( indexPrediction, inputPredictionGrid ) {
    
    #prfPar <- c(0, 0, 1.25, 1000, 6, 12, 0, 0.75, 1000, 135) #constraint: sigma pos 1 << sigma pos 2
    #indexPrediction <- 80
    #inputPredictionGrid <- predictionGrid
    
    prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
    
    a <- dnorm( x, prfPar[1], prfPar[3] )
    b <- dnorm( y, prfPar[2], prfPar[3] )
    imgCenter <- scaleData( tcrossprod(a,b), 1, 0 )
    
    a <- dnorm( x, prfPar[1], prfPar[4] )
    b <- dnorm( y, prfPar[2], prfPar[4] )
    imgSurround <- scaleData( tcrossprod(a,b), 1, 0)*prfPar[7]
    
    r <- imgCenter - imgSurround
    rMat <- array(r)
    
    predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
    pConv <- conv( predictionLoop, hrf )
    pConvTrim <- pConv[ 1 : dim(stimSeq_eyeFix)[3] ]
    tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('nearest') )
    returnPrediction <- round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 ) #### scale predictions betweeo 0 and 1 and round them to 5 digits
    
    # # to plot
    # par(mfrow=c(2,2))
    # plot( r[,70], type='l' )
    # image( r, col=gray.colors(500) )
    # plot( returnPrediction, type='l' )
    
    return( returnPrediction ) 
  }
  
  #### generate predictions in parallel, eye fix ####
  detectCores()
  nCores <- 6
  cl <- makeCluster(nCores, type='FORK')
  storeTimePar <- system.time( tsPredictionTransposed_eyeFix <- parSapply(cl, 1:dim(predictionGrid)[1], generatePrediction_visual, inputPredictionGrid=predictionGrid ) )
  stopCluster(cl)
  print( storeTimePar )

  #### prepare to generate the predictions, eye moving ####
  stimSeqMat <- array( stimSeq_eyeFix, c( length(x)*length(y), dim(stimSeq_eyeFix)[3] ) )
  incrementalCounter <- dim(predictionGrid)[1] * 0.05
  timeStimuli <- seq( 0, dim(stimSeq_eyeFix)[3]*samplingTime, length.out = dim(stimSeq_eyeFix)[3] )
  mriTime <- seq( 0, dim(stimSeq_eyeFix)[3]*samplingTime, length.out = dim(ts$brk)[4] )
  
  #### function to generate the predictions, eye moving ####
  generatePrediction_moving <- function( indexPrediction, inputPredictionGrid ) {
    
    #prfPar <- c(0, 0, 1.25, 1000, 6, 12, 0, 0.75, 1000, 135) #constraint: sigma pos 1 << sigma pos 2
    #indexPrediction <- 1
    #inputPredictionGrid <- predictionGrid
    
    prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
    
    a <- dnorm( x, prfPar[10], prfPar[12] )
    b <- dnorm( y, prfPar[11], prfPar[12] )
    imgCenter <- scaleData( tcrossprod(a,b), 1, 0 )
    
    #a <- dnorm( x, prfPar[1], prfPar[4] )
    #b <- dnorm( y, prfPar[2], prfPar[4] )
    #imgSurround <- scaleData( tcrossprod(a,b), 1, 0)*prfPar[7]
    
    r <- imgCenter
    rMat <- array(r)
    
    predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
    pConv <- conv( predictionLoop, hrf )
    pConvTrim <- pConv[ 1 : dim(stimSeq_eyeFix)[3] ]
    tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('nearest') )
    returnPrediction <- round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 ) #### scale predictions betweeo 0 and 1 and round them to 5 digits
    
    # # to plot
    # par(mfrow=c(2,2))
    # plot( r[,70], type='l' )
    # image( r, col=gray.colors(500) )
    # plot( returnPrediction, type='l' )
    
    return( returnPrediction ) 
  }
  
  #### generate predictions in parallel, eye moving ####
  detectCores()
  nCores <- 6
  cl <- makeCluster(nCores, type='FORK')
  storeTimePar <- system.time( tsPredictionTransposed_eyeMoving <- parSapply(cl, 1:dim(predictionGrid)[1], generatePrediction_moving, inputPredictionGrid=predictionGrid ) )
  stopCluster(cl)
  print( storeTimePar )
  
  #### clean up ts predictions and prediction grid ####
  tsPrediction_eyeFix <- t( tsPredictionTransposed_eyeFix )
  tsPrediction_eyeMoving <- t( tsPredictionTransposed_eyeMoving )
  tsPrediction_direction <- t( eyePredictions_direction )
  controlPredictions <- apply( tsPrediction_eyeFix, 1, sum ) != 0 & apply( tsPrediction_eyeMoving, 1, sum ) != 0 & apply( tsPrediction_direction, 1, sum ) != 0 
  tsPrediction_eyeFix <- tsPrediction_eyeFix[controlPredictions,]
  tsPrediction_eyeMoving <- tsPrediction_eyeMoving[controlPredictions,]
  tsPrediction_direction <- tsPrediction_direction[controlPredictions,]
  predictionGrid <- predictionGrid[controlPredictions,]
  print( dim( predictionGrid ) )
  
  
  #### linear fitting in parallel ####
  print( sprintf( 'linear fitting in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )

  totalVoxelsIterations <- ceil( length(selIdxVoxel) / 200 )
  limitsSelVoxels <- round( seq(1,length(selIdxVoxel), length.out=totalVoxelsIterations) ) #split between chuncks for parallel fitting (about 100 voxels each)
  limitsSelVoxels[1] <- 1
  limitsSelVoxels[ length(limitsSelVoxels) ] <- length(selIdxVoxel)
  limitsSelVoxelsMatrix <- array( 0, c( (length(limitsSelVoxels)-1) , 2 ) )
  limitsSelVoxelsMatrix[,1] <- limitsSelVoxels[1:(length(limitsSelVoxels)-1)]
  limitsSelVoxelsMatrix[,2] <- limitsSelVoxels[2:(length(limitsSelVoxels))]
  limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] <- limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk
  runIndex <- seq( 1:dim(limitsSelVoxelsMatrix)[1] )
  voxelModel <- function(passIdx) { #this fits the model on serveral voxels at a time, see limitsSelVoxelsMatrix
    selTsVoxel <- tsTransposedSel[ , limitsSelVoxelsMatrix[passIdx,1]:limitsSelVoxelsMatrix[passIdx,2] ]
    selTsVoxelMean <- apply( selTsVoxel, 2, mean ) #voxels average
    ssTot <- apply( (selTsVoxel-selTsVoxelMean)^2, 2, sum) #voxels total sum of squares
    runLinMod <- function(nIndex) { #get best fit for every prediction (nIndex = counter of predictions)
      dMat <- cbind( tsPrediction_eyeFix[nIndex,], tsPrediction_eyeMoving[nIndex,], tsPrediction_direction[nIndex,] )  #visual, movement and direction predictors    
      dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (4: intercept, beta eyeFix, beta eyeMoving, beta direction
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot) #r squares
      
      # dfTot <- dim(dMat01)[1]-dim(dMat01)[2]
      # seRoutine <- function(currIdxPred) {
      #    seComp <- sqrt( diag( ssRes[currIdxPred]/dfTot*matrix.power( ( t(dMat01) %*% dMat01 ) , -1 ) ) ) # standard error based without lm
      # }
      # seCoeff <- sapply( 1:length(ssRes), seRoutine )
      # tStats <- a/seCoeff
      # pVals <- 2 * pt(-abs(tStats), dfTot)
      #  
      # # to test
      # nVox <- 50
      # r2[nVox]
      # rbind( coef( lm( selTsVoxel[,nVox] ~ dMat ) ), a[,nVox] ) #these two lines should be VERY similar
      # fit = lm( selTsVoxel[,nVox] ~ dMat )
      # se <- sqrt(diag(vcov(fit)))
      # se # standard error based on the lm function in r, slow
      # seCoeff[,nVox] # standard error based on the analytical formula
      # # S^2{b} = MSE * (X^T * X)^-1 # the analytical formula for standard error in multiple linear regression
      # # t stats
      # summary(fit)$coefficients[,3]
      # tStats[,nVox]
      # # pvals
      # summary(fit)$coefficients[,4]
      # pVals[,nVox]
      
      #storePredictorVisual <- array( rep( dMat[,1], length(r2) ), c( length(dMat[,1]), length(r2) ) ) 
      #storePredictorMoving <- array( rep( dMat[,2], length(r2) ), c( length(dMat[,2]), length(r2) ) ) 
      #storePredictorDirection <- array( rep( dMat[,3], length(r2) ), c( length(dMat[,3]), length(r2) ) ) 
      return( rbind( r2, a, expectedTs ) )
    }
    outVoxelList <- lapply( 1:dim(tsPrediction_direction)[1], runLinMod  ) #apply the function for all predictions, get a list as output
    outVoxel3D <- abind( outVoxelList, along=3 ) #reshape the list in a 3dmatrix with fields: output X voxels tested X predictions
    betaPositiveMatrix <- outVoxel3D[3,,] > 0 & outVoxel3D[4,,] > 0 & outVoxel3D[5,,] > 0 # each beta must be positive
    r2Matrix <- outVoxel3D[1,,]	
    extractData <- function(nSelectedVoxels) {
      indexBetaZero <- betaPositiveMatrix[nSelectedVoxels,]
      if ( sum(indexBetaZero)>0 ) {
        indexBetaZero <- which( indexBetaZero )
        indexVarExp <- which.max( r2Matrix[nSelectedVoxels,indexBetaZero] )
        return( as.numeric( c( predictionGrid[indexBetaZero[indexVarExp],], outVoxel3D[,nSelectedVoxels,indexBetaZero[indexVarExp]] ) ) ) 
      }
      if ( sum(indexBetaZero)==0 ) {
        return( rep(0, dim(predictionGrid)[2]+5+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 5=r2,beta intercept, beta slope visual, beta slope movement, beta slope direction 
      }
    }	
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }
  #system.time( aaa <- lapply( 1:20, voxelModel ) )
  detectCores()
  nCores <- 6
  cl <- makeCluster(nCores, type='FORK')
  storeTimePar <- system.time( outModel <- parLapply(cl, runIndex, voxelModel ) )
  stopCluster(cl)
  print( storeTimePar )
  
  outMatrix <- array(0, c( dim(outModel[[1]])[1], length(selIdxVoxel) ) )
  for ( nElements in 1:length(outModel) ) {
    tempNElementsStart <- limitsSelVoxelsMatrix[nElements,1]
    tempNElementsEnd <- limitsSelVoxelsMatrix[nElements,2]
    outMatrix[,tempNElementsStart:tempNElementsEnd] <- outModel[[nElements]]
  }
  
  outModel <- outMatrix
  if (modelFitCounter==1) { outModelLoop <- outModel }
  if (modelFitCounter>1) { 
    selectedCols <- outModel[13,] > outModelLoop[13,] #where r2 is store (after 12 model parameters, there are: r2, beta intercept, beta slope visual, beta slope movement, beta slope direction 
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
}

# get single parameter contribution
print('get single parameter contribution...')
fittedPredictionGrid <- t( outModelLoop[seq(1:12),] )
singlePredictorStats <- function( idxLine, inputMat ) {
  
  #idxLine <- 57
  #inputMat <- fittedPredictionGrid
  
  if ( sum(inputMat[idxLine,])!=0 ) {
    
    directionPrediction <- eyeDirectionPrediction(idxLine,inputMat)
    gainPrediction <- generatePrediction_moving(idxLine,inputMat)
    eyeCenterPrediction <- generatePrediction_visual(idxLine,inputMat)
    tsToFit <- tsTransposedSel[,idxLine]
    
    ssTotSingleVoxel <- sum ( (tsToFit-mean(tsToFit))^2 )  #voxel total sum of squares
    dMat <- cbind( directionPrediction, gainPrediction, eyeCenterPrediction )  #visual, movement and direction predictors    
    dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
    a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,tsToFit) ) #beta coefficients (4: intercept, beta eyeFix, beta eyeMoving, beta direction
    expectedTs <- crossprod( t(dMat01), a ) #expected ts
    residualsSquared <- ( tsToFit - expectedTs )^2 
    ssRes <- sum( residualsSquared )
    r2 <- 1-(ssRes/ssTotSingleVoxel) #r squares
    
    dfTot <- dim(dMat01)[1]-dim(dMat01)[2]
    seCoeff <- sqrt( diag( ssRes/dfTot*matrix.power( ( t(dMat01) %*% dMat01 ) , -1 ) ) ) # standard error based without lm
    tStats <- a/seCoeff
    pVals <- 2 * pt(-abs(tStats), dfTot)
  }
  
  if ( sum(inputMat[idxLine,])==0 ) {
    r2 <- 0
    tStats <- rep(0,4)
    pVals <- rep(0,4)
  }
  
  return( c(r2, tStats, pVals) )
  
  ## to test
  # nVox <- 1
  # r2[nVox]
  # rbind( coef( lm( selTsVoxel[,nVox] ~ dMat ) ), a[,nVox] ) #these two lines should be VERY similar
  # fit = lm( selTsVoxel[,nVox] ~ dMat )
  # se <- sqrt(diag(vcov(fit)))
  # se # standard error based on the lm function in r, slow
  # seCoeff[nVox] # standard error based on the analytical formula
  # # S^2{b} = MSE * (X^T * X)^-1 # the analytical formula for standard error in multiple linear regression
  # # t stats
  # summary(fit)$coefficients[,3]
  # tStats
  # # pvals
  # summary(fit)$coefficients[,4]
  # pVals
  
}
#ee <- sapply(55:60, singlePredictorStats, inputMat=fittedPredictionGrid )
runIndex <- 1:dim(fittedPredictionGrid)[1]
detectCores()
nCores <- 4
cl <- makeCluster(nCores, type='FORK')
storeTimePar <- system.time( outModelStats <- parSapply(cl, runIndex, singlePredictorStats, inputMat=fittedPredictionGrid ) )
stopCluster(cl)
print( storeTimePar )
#sum( outModelStats[1,] - outModelLoop[13,c(1:dim(outModelStats)[2])] < 0.0001 )
#outModelStats <- array(0,c( dim(fittedPredictionGrid)[1], 9 ))
storeAllPred <- cbind( t(outModelLoop[ seq(1,17), ]), outModelStats )

#### extract surround size ####
print('surround size...')
if (flagSurround==1) { 
  print('get FWHM...')
  progress <- 0.05
  FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
  xFWHM <- seq(-25,25,0.0125)
  for (k in 1:dim(storeAllPred)[1]) {
    
    parameters <- storeAllPred[k,]
    
    if ( sum(parameters!=0) ) {
      a <- dnorm( xFWHM, 0, parameters[3] )
      b <- dnorm( xFWHM, 0, parameters[4] )
      r <- scaleData(a,1,0) - scaleData(b,1,0)*parameters[7]
      
      maxIdx <- which.max(r)
      rHalf <- r[ maxIdx:length(xFWHM) ]
      xHalf <- xFWHM[ maxIdx:length(xFWHM) ]
      halfMax <- max(rHalf)/2
      xMid <- xHalf[ rHalf<=halfMax  ]  
      FWHMcenter <- min( xMid )*2 
      if ( sum(r<0)>0 ) { #if there is a detectable surround, compute surround size
        rMin <- which.min(rHalf)
        FWHMsurround <- xHalf[rMin[1]]*2
      }
      if ( sum(r<0)==0 ) { #otherwise set it to zero
        FWHMsurround <- 0
      }
      FWHM[k,] <- c(FWHMcenter,FWHMsurround)
    }
    
    if (k>dim(storeAllPred)[1]*progress) {
      cat(paste('*',""))
      progress <- progress + 0.05
    }
  }
}
if (flagSurround==0) { 
  FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
  FWHM[,1] <- storeAllPred[,3]
  FWHM[,2] <- 1000
}

storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 18:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(30, dim(tsTransposedAll)[2] ) )

print('save linear step...')
polCoords <- cart2pol( storeAllPred[,c(2,1)] ) 
storeAllPredOut[,selIdxVoxel] <- t( cbind( storeAllPred, polCoords, FWHM ) )
fileTs <- sprintf('%s_PredixtedTs.nii.gz',outSuffix) 
fileParams <- sprintf('%s_params.nii.gz',outSuffix)

rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( '__tt_expectedTs.nii.gz', brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_expectedTs.nii.gz', orientValue[1,1], fileTs )
system( instr)
system( 'rm __tt_expectedTs.nii.gz' )

rParameters <- array( storeAllPredOut, c( dim(storeAllPredOut)[1], dim(ts$brk)[1:3] ) )
rParameters <- aperm( rParameters, c(2,3,4,1) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', epiFile )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( '__tt_parameters.nii.gz', brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_parameters.nii.gz', orientValue[1,1], fileParams )
system( instr)
system( 'rm __tt_parameters.nii.gz' )

labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multPar','centerVonMis','kVonMis','xGain','yGain','sigmaGain',
            'varExp','intercept','b_eyeFix','b_eyeMoving','b_direction','varExp_cont','t_intercept','t_eyeFix','t_eyeMoving','t_direction',
            'p_intercept','p_eyeFix','p_eyeMoving','p_direction','theta','radius','FWHM_center','surround_size')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}


