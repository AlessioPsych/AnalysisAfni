args <- commandArgs(T)
print( args )

## prepare the stimuli with the portrait as well

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')


#args <- c('maxVarEye.nii.gz', 'meanTs_eye_topUp_res.nii', 'output_del', '5', '3','0.166','3','0')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
#instr <- '3dresample -dxyz 4 4 4 -orient RAI -rmode Lin -prefix barsTs_res.nii.gz -inset meanTsBars_topUp.nii'
#system( instr )
#instr <- '3dTstat -mean -prefix bars_mean_res.nii.gz barsTs_res.nii.gz'
#system( instr )
#instr <- '3dAutomask -prefix bars_mask_res.nii.gz bars_mean_res.nii.gz'
#system( instr )
#instr <- '3dresample -dxyz 4 4 4 -orient RAI -rmode Lin -prefix eyeTs_res.nii.gz -inset meanTsEye_topUp.nii'
#system( instr )




source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( neuRosim )
library( parallel )

epiFile <- args[1]
inputFile <- args[2]
outSuffix <- args[3]
fineFit <- as.numeric( args[4] )
polortArg <- as.numeric( args[5] )
samplingTime <- as.numeric( args[6] )
stimType <- as.numeric( args[7] )
flagSurround <-  as.numeric( args[8] )
#fitGain <- as.numeric( args[9] )

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
print('get prfStimuli.RData...')
#load( file='prfStimuli.RData' )
#setwd('/analyse/Project0226/scanner_train - Copy/orientationgrating_points_fit_long_glas_large_get_stimuli')
#arrayStim <- scan( 'eyeMovingStim.txt' )
#setwd('/analyse/Project0226/scanner_train - Copy/orientationgrating_points_fit_long_glas_large_fixation_get_stimuli_border')
#arrayStim <- scan( 'eyeFixStim_border.txt' )

setwd('/analyse/Project0226/dataSummary')
if (stimType==1) { 
  arrayStim <- scan( 'eyeMovingStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==2) { 
  arrayStim <- scan( 'eyeFixStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==3) { 
  arrayStim <- scan( 'eyeFixStim_border.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==4) { 
  arrayStim <- scan( 'prfStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1860,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==5) { 
  arrayStim <- scan( 'eyeFixStim_border_disappear.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}


#stimMat <- aperm( array( arrayStim, c(128,620,96) ), c( 1, 3, 2 ) )
#stimMat <- aperm( array( arrayStim, c(200,1860,150) ), c( 3, 1, 2 ) ) # prf
#stimMatInterp2 <- array( 0,c(60,60,1510) )
#xStim <- seq(0,200,length.out=200)	
#yStim <- seq(0,150,length.out=150)	
#xStimInt <- seq(0,200,length.out=60)	
#yStimInt <- seq(0,150,length.out=60)	
#for (interpCounter in 1:dim(stimMat)[2] ) {
#	imageLoop <- stimMat[,,interpCounter]
#	stimMatInterp2[,,interpCounter] <- interp2( xStim, yStim, imageLoop, xStimInt, yStimInt, method=c('linear') )
#}

#stimMatFlip <- stimMat[ ,dim(stimMat)[2]:1, ]


stimMatFlip <- aperm( stimMat[ dim(stimMat)[1]:1,, ], c(2,1,3) )

x11( width=3, height=3 )
for ( snap in 1:dim(stimMat)[3] ) {
  image( stimMatFlip[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01)
}
stimSeq <- stimMatFlip
x <- seq(-10,10,length.out = dim(stimSeq)[1] )
y <- seq(-5.5,5.5,length.out = dim(stimSeq)[2] )
outMesh <- meshgrid(y, x)
outMesh$X <- scaleData( outMesh$X, 1, 0 )
outMesh$Y <- scaleData( outMesh$Y, 1, 0 )

rm( stimMat )

# stimImageTemp <- array( 0, c( dim(stimMatFlip)[1]+100, dim(stimMatFlip)[2]+100, dim(stimMatFlip)[3] ) )
# x11( width=3, height=3 )
# for ( snap in 1:dim(stimMat)[3] ) {
#   imageTemp <- stimMatFlip[,,snap]
#   imageTemp <- cbind( matrix(1, nrow = dim(imageTemp)[1], 50 ), imageTemp, matrix(1, nrow = dim(imageTemp)[1], 50 ) )
#   imageTemp <- rbind( matrix(1, ncol = dim(imageTemp)[2], 50 ), imageTemp, matrix(1, ncol = dim(imageTemp)[2], 50 ) )
#   image( imageTemp, axes=FALSE ); par(new=TRUE); Sys.sleep(0.01)
#   stimImageTemp[,,snap] <- imageTemp
# }
# stimSeq <- stimImageTemp
# x <- seq(-12,12,length.out = dim(stimSeq)[1] )
# y <- seq(-10,10,length.out = dim(stimSeq)[2] )
# rm( stimMat )

#this part of the code builds a matrix with all the possible prediction tested, for both models at this stage

if (fineFit==-2) {
  xElements <- 4
  yElements <- 4
  sigmaArrayPositiveElements <- 4
  multParElements <- 3
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
}
if (fineFit==-1) {
  xElements <- 7
  yElements <- 7
  sigmaArrayPositiveElements <- 5
  multParElements <- 4
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
}

if (fineFit==0) {
  xElements <- 4
  yElements <- 4
  sigmaArrayPositiveElements <- 4
  multParElements <- 3
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  xGainElements <- 3
  yGainElements <- 3
}
if (fineFit==1) {
  xElements <- 8
  yElements <- 8
  sigmaArrayPositiveElements <- 9
  multParElements <- 4
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  xGainElements <- 5
  yGainElements <- 5
}

if (fineFit==2) {
  xElements <- 4
  yElements <- 4
  sigmaArrayPositiveElements <- 4
  multParElements <- 3
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  xGainPosElements <- 2
  yGainPosElements <- 2
  sizeGainElements <- 2
}
if (fineFit==3) {
  xElements <- 7
  yElements <- 7
  sigmaArrayPositiveElements <- 5
  multParElements <- 4
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  xGainPosElements <- 4
  yGainPosElements <- 4
  sizeGainElements <- 4
}
if (fineFit==4) { #oblique prfs coarse
  xElements <- 6
  yElements <- 6
  sigmaArrayPositiveElements <- 4
  multParElements <- 3
  hrfDelayOnsetElements <- 3
  hrfDelayUnderShootElements <- 3
  sigmaArrayPositiveElements_2 <- 4
  thetaElements <- 4
}
if (fineFit==5) { #oblique prfs fine
  xElements <- 14
  yElements <- 14
  sigmaArrayPositiveElements <- 8
  multParElements <- 4
  hrfDelayOnsetElements <- 2
  hrfDelayUnderShootElements <- 2
  sigmaArrayPositiveElements_2 <- 8
  thetaElements <- 6
}

print('build prediction...')
xPosFit <- seq( -9, 9, length.out=xElements )
yPosFit <- seq( -4.5, 4.5, length.out=yElements )
sigmaArrayPositive <- seq( 0.25, 7, length.out=sigmaArrayPositiveElements )
#xPosFit <- seq( -2, 2, length.out=xElements )
#yPosFit <- seq( -2, 2, length.out=yElements )
#sigmaArrayPositive <- seq( 0.25, 2, length.out=sigmaArrayPositiveElements )
if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
par_hrf_a1 <- seq( 6, 10, length.out=hrfDelayOnsetElements )
par_hrf_a2 <- seq( 12, 16, length.out=hrfDelayUnderShootElements )
if (flagSurround==1) { multPar <- seq(0,0.8, length.out = multParElements) }
if (flagSurround==0) { multPar <- 0 }

if ( fineFit==-1 | fineFit==-2  ) { # no gain, just prf basically, coarse (-1) and fine (-2)
  xGain <- 0
  yGain <- 0
  addPredictor <- 0
  predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, xGain, yGain, addPredictor )
}
if ( fineFit==0 | fineFit==1 ) { # if coarse fit or fine fit with plane gain field
  xGain <- seq(-1,1,length.out = xGainElements)
  yGain <- seq(-1,1,length.out = yGainElements)
  addPredictor <- 0
  predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, xGain, yGain, addPredictor )
}
if ( fineFit==2 | fineFit==3 ) { #if coarse fit or fine fit with gaussian plane gain field
  xGain <- seq(-9,9,length.out = xGainPosElements)
  yGain <- seq(-4.5,4.5,length.out = yGainPosElements)
  sigmaArrayGain <- seq( 2, 8, length.out=sizeGainElements )
  predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, xGain, yGain, sigmaArrayGain )
}
if (fineFit==4 | fineFit==5 ) {
  sigmaArrayPositive_02 <- seq( 0.25, 7, length.out=sigmaArrayPositiveElements_2 )
  if (flagSurround==1) { sigmaArrayNegative_02 <- sigmaArrayPositive }
  if (flagSurround==0) { sigmaArrayNegative_02 <- 1000 }
  thetaParameter <- seq( 0, 170, length.out = thetaElements )
  predictionGridTemp01 <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, sigmaArrayPositive_02, sigmaArrayNegative_02, thetaParameter )
  keepPredictionIdxTemp01 <- (predictionGridTemp01[ ,3] < predictionGridTemp01[ ,4]) & #positive sigmas01 < than negative sigmas01
    (predictionGridTemp01[ ,8] < predictionGridTemp01[ ,9]) & #positive sigmas02 < than negative sigmas02
    (predictionGridTemp01[ ,3] <= predictionGridTemp01[ ,8]) #positive sigmas01 < than positive sigmas02
  predictionGridTemp02 <- predictionGridTemp01[keepPredictionIdxTemp01,]

  thetaFactor <- as.factor( predictionGridTemp02[,10] )
  thetaFlag <- thetaFactor==levels(thetaFactor)[1]
  keepPredictionIdxTemp02_01 <- thetaFlag & (predictionGridTemp02[,3] == predictionGridTemp02[,8]) #in case of identical positive sigmas remove all the thetas, prf is isotropic in this case
  keepPredictionIdxTemp02_02 <- (predictionGridTemp02[,3] != predictionGridTemp02[,8]) #also keep all the non-isotropic prfs
  predictionGridTemp02 <- predictionGridTemp02[  keepPredictionIdxTemp02_01 | keepPredictionIdxTemp02_02, ]
}
if (fineFit<4) {
  keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4] 
  predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]
}
if (fineFit==4 | fineFit==5) {
  predictionGridGlobal <- predictionGridTemp02
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
for (modelFitCounter in c(1,2,3,4,388,389,390,391) ) { #1:length(runIndexPredictions)
  
  print( sprintf( 'iteration %1.0f of %1.0f, start...', modelFitCounter, length(runIndexPredictions)  ) )
  
  # here I select the portion of all the predictions that are going to be tested later
  predictionGrid <- predictionGridGlobal[ limitsPredictionMatrix[modelFitCounter,1]:limitsPredictionMatrix[modelFitCounter,2], ]
  
  print( sprintf( 'generate predictions in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  
  #### prepare to generate the predictions ####
  stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
  incrementalCounter <- dim(predictionGrid)[1] * 0.05
  timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
  mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(ts$brk)[4] )
  
  #### function to generate the predictions ####
  generatePrediction <- function( indexPrediction, inputPredictionGrid ) {
    prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    
    #prfPar <- c(0, 0, 1.25, 1000, 6, 12, 0, 0.75, 1000, 135) #constraint: sigma pos 1 << sigma pos 2
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
    
    a <- dnorm( x, prfPar[1], prfPar[3] )
    b <- dnorm( y, prfPar[2], prfPar[3] )
    imgCenter <- scaleData( tcrossprod(a,b), 1, 0 )
    
    a <- dnorm( x, prfPar[1], prfPar[4] )
    b <- dnorm( y, prfPar[2], prfPar[4] )
    imgSurround <- scaleData( tcrossprod(a,b), 1, 0)*prfPar[7]

    if (fineFit==-1 | fineFit==-2) { #no gain
      r <- imgCenter - imgSurround
    }
    if (fineFit==0 | fineFit==1 ) { #plane gain
      imgGain <- outMesh$X*prfPar[8] + outMesh$Y*prfPar[9]
      r <- (imgCenter - imgSurround) + imgGain
    }
    if (fineFit==2 | fineFit==3) { #gaussian gain
      a <- dnorm( x, prfPar[8], prfPar[10] )
      b <- dnorm( y, prfPar[9], prfPar[10] )
      imgGain <- scaleData( tcrossprod(a,b), 1, 0 )
      r <- (imgCenter - imgSurround) + imgGain
    }
    if (fineFit==4 | fineFit==5) {
      meshOut = meshgrid(y, x);
      X1 <- meshOut$X
      X2 <- meshOut$Y
      sigma1 <- prfPar[8];
      sigma2 <- prfPar[3];
      Theta <- prfPar[10];
      mu01 <- prfPar[2];
      mu02 <- prfPar[1];
      
      a <- ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
      b <- -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
      c <- ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));
      
      A <- 1;
      rPos <- A*exp(-(a*(X1 - mu01)^2 + 2*b*(X1 - mu01)*(X2 - mu02) + c*(X2 - mu02)^2))
      
      sigma1 <- prfPar[9];
      sigma2 <- prfPar[4];
      Theta <- prfPar[10];
      mu01 <- prfPar[2];
      mu02 <- prfPar[1];
      
      a <- ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
      b <- -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
      c <- ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));
      
      A <- 1;
      rNeg <- A*exp(-(a*(X1 - mu01)^2 + 2*b*(X1 - mu01)*(X2 - mu02) + c*(X2 - mu02)^2))
      r <- rPos-rNeg*prfPar[7]
    }
    
    r[ r<quantile(r,0.1) ] <- 0
    r <- scaleData(r,1,0)
    rMat <- array(r)
    
    
    #trMat <- t(stimSeqMat)
    #system.time( predictionLoop <- as.numeric( trMat%*%rMat ) )
    predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
    pConv <- conv( predictionLoop, hrf )
    pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
    tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('nearest') )
    returnPrediction <- round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 ) #### scale predictions betweeo 0 and 1 and round them to 5 digits
    #par(mfrow=c(2,2))
    #plot( r[,70], type='l' )
    #image( r, col=gray.colors(500) )
    #plot( returnPrediction, type='l' )
    #image( imgGain, col=gray.colors(500)  )
    
    return( returnPrediction ) 
  }
  
  #### generate predictions in parallel ####
  detectCores()
  nCores <- 3
  cl <- makeCluster(nCores, type='FORK')
  storeTimePar <- system.time( tsPredictionTransposed <- parSapply(cl, 1:dim(predictionGrid)[1], generatePrediction, inputPredictionGrid=predictionGrid ) )
  stopCluster(cl)
  #storeTimeSerial <- system.time( sapply(1:100, generatePrediction, inputPredictionGrid=predictionGrid ) )
  print( storeTimePar )
  #print( storeTimeSerial )
  
  #### clean up ts predictions and prediction grid ####
  tsPrediction <- t( tsPredictionTransposed )
  controlPredictions <- apply( tsPrediction, 1, sum ) != 0 
  tsPrediction <- tsPrediction[controlPredictions,]
  predictionGrid <- predictionGrid[controlPredictions,]
  print( dim( predictionGrid ) )
  
  
  #### linear fitting in parallel ####
  print( sprintf( 'linear fitting in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  # indexVol <- meanEpi$brk[,,,1]
  # indexArray <- array( indexVol, prod( dim( indexVol ) ) )
  # tsTransposedAll <- t( tsArray )
  # selIdxVoxel <- which( indexArray == 1 )
  # tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series
  
  totalVoxelsIterations <- ceil( length(selIdxVoxel) / 300 )
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
      dMat <- cbind( tsPrediction[nIndex,] )      
      dMat01 <- cbind( rep(1,length(dMat)), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope)
      #a <- solve( qr(dMat01), selTsVoxel)
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot) #r squares
      return( rbind( r2, a, expectedTs ) )
    }
    outVoxelList <- lapply( 1:dim(tsPrediction)[1], runLinMod  ) #apply the function for all predictions, get a list as output
    outVoxel3D <- abind( outVoxelList, along=3 ) #reshape the list in a 3dmatrix with fields: output X voxels tested X predictions
    betaPositiveMatrix <- outVoxel3D[3,,] > 0
    r2Matrix <- outVoxel3D[1,,]	
    extractData <- function(nSelectedVoxels) {
      indexBetaZero <- betaPositiveMatrix[nSelectedVoxels,]
      if ( sum(indexBetaZero)>0 ) {
        indexBetaZero <- which( indexBetaZero )
        indexVarExp <- which.max( r2Matrix[nSelectedVoxels,indexBetaZero] )
        return( as.numeric( c( predictionGrid[indexBetaZero[indexVarExp],], outVoxel3D[,nSelectedVoxels,indexBetaZero[indexVarExp]] ) ) ) 
      }
      if ( sum(indexBetaZero)==0 ) {
        return( rep(0, dim(predictionGrid)[2]+3+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 3=r2,beta intercept, beta slope
      }
    }	
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }
  #system.time( aaa <- lapply( 1:2, voxelModel ) )
  detectCores()
  nCores <- 3
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
  
  # for ( nElements in 1:length(outModel) ) {
  #   if (nElements==1) {
  #     outMatrix <- outModel[[nElements]]
  #   }
  #   if (nElements>1) {
  #     outMatrix <- cbind( outMatrix, outModel[[nElements]] )
  #   }
  # }
  
  outModel <- outMatrix
  if (modelFitCounter==1) { outModelLoop <- outModel }
  if (modelFitCounter>1) { 
    selectedCols <- outModel[11,] > outModelLoop[11,] #where r2 is store (after 4 prf parameters and 2 hrf parameters, surround parameter and 3 gain field parameters), there remain the intercet and the slope parameter
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
}

#outModelLoop <- outMatrix #fix it from here
storeAllPred <- t(outModelLoop[ seq(1,13), ])
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
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 14:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(17, dim(tsTransposedAll)[2] ) )

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

if (fineFit==4 | fineFit==5 ) {
  labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multPar','sigmaPos_1','sigmaNeg_1','theta','varExp','intercept','slope','theta','radius','fwhmCenter','surroundSize')
  for (k in 1:length(labels)) {
    instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
    print( instr )
    system( instr )
  }
} 
if (fineFit!=4 & fineFit!=5 ) {
  labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multPar','gainX','gainY','sigmaGain','varExp','intercept','slope','theta','radius','fwhmCenter','surroundSize')
  for (k in 1:length(labels)) {
    instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
    print( instr )
    system( instr )
  }
} 


