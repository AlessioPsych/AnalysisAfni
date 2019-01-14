args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')
#args <- c('meanTs_eyeMovement_topUp_res_mask.nii', 'meanTs_eyeMovement_topUp_res.nii', '__delme_output_mixed','/home/alessiof/abin', '1', '4','1','0.166','1')

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
flagSurround <- as.numeric( args[5] )
fitIntercept <- as.numeric( args[7] )
polortArg <- as.numeric( args[6] )
samplingTime <- as.numeric( args[8] )
stimType <- as.numeric( args[9] )
fineFit <- 0

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
#if (stimType==1) { 
arrayStim_eyeMovingStim <- scan( 'eyeMovingStim.txt' )
stimMat_eyeMovingStim <- aperm( array( arrayStim_eyeMovingStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
#}
#if (stimType==2) { 
#  arrayStim <- scan( 'eyeFixStim.txt' )
#  setwd(mainDir)
#  stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
#}
#if (stimType==3) { 
arrayStim_eyeFixStim_border <- scan( 'eyeFixStim_border.txt' )
stimMat_eyeFixStim_border <- aperm( array( arrayStim_eyeFixStim_border, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
setwd(mainDir)
#}
#if (stimType==4) { 
#  arrayStim <- scan( 'prfStim.txt' )
#  setwd(mainDir)
#  stimMat <- aperm( array( arrayStim, c(200,1860,150) ), c( 3, 1, 2 ) ) # eye movement
#}

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

#stimMatFlip_eyeMovingStim <- stimMat_eyeMovingStim[ ,dim(stimMat_eyeMovingStim)[2]:1, ]
#stimMatFlip_eyeFixStim_border <- stimMat_eyeFixStim_border[ ,dim(stimMat_eyeFixStim_border)[2]:1, ]
stimMatFlip_eyeMovingStim <- aperm( stimMat_eyeMovingStim[ dim(stimMat_eyeMovingStim)[1]:1,, ], c(2,1,3) )
stimMatFlip_eyeFixStim_border <- aperm( stimMat_eyeFixStim_border[ dim(stimMat_eyeFixStim_border)[1]:1,, ], c(2,1,3) )

x11( width=3.5, height=3.5 )
for ( snap in 1:dim(stimMat_eyeMovingStim)[3] ) {
  image( stimMatFlip_eyeMovingStim[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01) 
}
x11( width=3.5, height=3.5 )
for ( snap in 1:dim(stimMat_eyeFixStim_border)[3] ) {
  image( stimMatFlip_eyeFixStim_border[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01) 
  #image( stimSeq_eyeFixStim_border[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.001) 
}

stimSeq_eyeMovingStim <- stimMatFlip_eyeMovingStim
stimSeq_eyeFixStim_border <- stimMatFlip_eyeFixStim_border
x <- seq(-6,6,length.out = dim(stimSeq_eyeMovingStim)[1] )
y <- seq(-8,8,length.out = dim(stimSeq_eyeMovingStim)[2] )
rm( list=c('stimMatFlip_eyeFixStim_border','stimMatFlip_eyeMovingStim','stimMat_eyeFixStim_border','stimMat_eyeMovingStim','arrayStim_eyeFixStim_border','arrayStim_eyeMovingStim') )

#### get the standard hrf ####
print('get hrf...')
hrf <- canonicalHRF( seq(0,30,samplingTime) )

# #this part of the code builds a matrix with all the possible prediction tested, for both models at this stage
# addSpace <- abs( min(x) )*0.5 
# print('build prediction...')
# xPosFit <- seq( min(x)-addSpace, max(x)+addSpace, length.out=4 )
# yPosFit <- seq( min(y)-addSpace, max(y)+addSpace, length.out=4 )
# sigmaArrayPositive <- seq( 0.25, 7, length.out=3 )
# if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
# if (flagSurround==0) { sigmaArrayNegative <- 1000 }
# predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative )
# keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
# predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]
# 
# # this part builds a matrix with the starting and ending prediction index to fit for each loop in the next for loop
# limitsPrediction <- round( seq(1,dim(predictionGridGlobal)[1], length.out=50) ) #split between 100 (arbitrary number) chuncks for iterative fitting
# limitsPrediction[1] <- 1
# limitsPrediction[ length(limitsPrediction) ] <- dim(predictionGridGlobal)[1]
# limitsPredictionMatrix <- array( 0, c( (length(limitsPrediction)-1) , 2 ) )
# limitsPredictionMatrix[,1] <- limitsPrediction[1:(length(limitsPrediction)-1)]
# limitsPredictionMatrix[,2] <- limitsPrediction[2:(length(limitsPrediction))]
# limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] <- limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
# runIndexPredictions <- seq( 1:dim(limitsPredictionMatrix)[1] )

modelType <- stimType
multVector <- c(0.5,1,1.2)#c(0.1, 0.5, 1, 1.25)
multVectorGrid <- expand.grid( multVector, multVector, multVector )
#modelFitCounter <- 1
for ( modelFitCounter in 1:dim( multVectorGrid )[1] ) {
  
  print( sprintf( 'iteration %1.0f of %1.0f, start...', modelFitCounter, dim( multVectorGrid )[1]  ) )
  
  #### prepare to generate the predictions ####
  print( sprintf( 'iteration %1.0f of %1.0f, building predictions...', modelFitCounter, dim( multVectorGrid )[1]  ) )
  for (stimuliCondition in 1:2) {
    if (stimuliCondition==1) { stimSeq <- stimSeq_eyeMovingStim }
    if (stimuliCondition==2) { stimSeq <- stimSeq_eyeFixStim_border }
    
    # here I select the portion of all the predictions that are going to be tested later
    #predictionGrid <- predictionGridGlobal[ limitsPredictionMatrix[modelFitCounter,1]:limitsPredictionMatrix[modelFitCounter,2], ]
    
    addSpace <- abs( min(x) )*0.05
    print('build prediction...')
    xPosFit <- seq( min(x)-addSpace, max(x)+addSpace, length.out=2 ) * multVectorGrid[ modelFitCounter, 1 ]
    yPosFit <- seq( min(y)-addSpace, max(y)+addSpace, length.out=2 ) * multVectorGrid[ modelFitCounter, 2 ]
    sigmaArrayPositive <- seq( 0.25, 6, length.out=4 ) * multVectorGrid[ modelFitCounter, 3 ]

    if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
    if (flagSurround==0) { sigmaArrayNegative <- 1000 }
    predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative )
    keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
    predictionGrid <- predictionGridTemp[ keepPredictionIdx, ]
    
    #tsPrediction <- array( 0, c( dim(predictionGrid)[1], dim( ts$brk )[4] ) )
    stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
    incrementalCounter <- dim(predictionGrid)[1] * 0.05
    timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
    mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(ts$brk)[4] )
    
    #### function to generate the predictions ####
    generatePrediction <- function( indexPrediction, inputPredictionGrid ) {
      prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
      
      #prfPar <- c(-2,-4,1,1.1)
      a <- dnorm( x, prfPar[1], prfPar[3] )
      b <- dnorm( y, prfPar[2], prfPar[3] )
      #a <- normFun( x, prfPar[1], prfPar[3] )
      #b <- normFun( y, prfPar[1], prfPar[3] )
      imgCenter <- tcrossprod(a,b)
      #image( imgCenter )
      
      a <- dnorm( x, prfPar[1], prfPar[4] )
      b <- dnorm( y, prfPar[2], prfPar[4] )
      #a <- normFun( x, prfPar[1], prfPar[4] )
      #b <- normFun( y, prfPar[1], prfPar[4] )
      imgSurround <- tcrossprod(a,b)
      
      if (flagSurround==1) { r <- imgCenter - imgSurround }
      if (flagSurround==0) { r <- imgCenter }
      rMat <- array(r)
      
      #trMat <- t(stimSeqMat)
      #system.time( predictionLoop <- as.numeric( trMat%*%rMat ) )
      predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
      pConv <- conv( predictionLoop, hrf )
      pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
      tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('linear') )
      #tsPrediction[nPrediction,] <- tsPredictionMriInterp / max( tsPredictionMriInterp )
      #return( round( tsPredictionMriInterp / max( tsPredictionMriInterp ), 5 ) ) #### scale predictions to 1 and round them to 5 digits
      return( round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 ) ) #### scale predictions betweeo 0 and 1 and round them to 5 digits
    }
    
    #### generate predictions in parallel ####
    library(parallel)
    detectCores()
    nCores <- 4
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
    if (stimuliCondition==1) { 
      tsPrediction_eyeMovingStim <- tsPrediction
      predictionGrid_eyeMovingStim <- predictionGrid
      rm('tsPrediction','predictionGrid')
    }
    if (stimuliCondition==2) { 
      tsPrediction_eyeFixStim_border <- tsPrediction
      predictionGrid_eyeFixStim_border <- predictionGrid
      rm('tsPrediction','predictionGrid')
    }
  }
  
  #### linear fitting in parallel ####
  print( sprintf( 'iteration %1.0f of %1.0f, linear fit...', modelFitCounter, dim( multVectorGrid )[1]  ) )
  indexVol <- meanEpi$brk[,,,1]
  indexArray <- array( indexVol, prod( dim( indexVol ) ) )
  tsTransposedAll <- t( tsArray )
  selIdxVoxel <- which( indexArray == 1 )
  tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series
  limitsSelVoxels <- round( seq(1,length(selIdxVoxel), length.out=100) ) #split between 50 (arbitrary number) chuncks for parallel fitting
  limitsSelVoxels[1] <- 1
  limitsSelVoxels[ length(limitsSelVoxels) ] <- length(selIdxVoxel)
  limitsSelVoxelsMatrix <- array( 0, c( (length(limitsSelVoxels)-1) , 2 ) )
  limitsSelVoxelsMatrix[,1] <- limitsSelVoxels[1:(length(limitsSelVoxels)-1)]
  limitsSelVoxelsMatrix[,2] <- limitsSelVoxels[2:(length(limitsSelVoxels))]
  limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] <- limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
  runIndex <- seq( 1:dim(limitsSelVoxelsMatrix)[1] )
  predictionIndexesMatrix <- expand.grid( 1:dim(predictionGrid_eyeMovingStim)[1], 1:dim(predictionGrid_eyeFixStim_border)[1] )
  voxelModel <- function( passIdx ) { #this fits the model on serveral voxels at a time, see limitsSelVoxelsMatrix
    selTsVoxel <- tsTransposedSel[ , limitsSelVoxelsMatrix[passIdx,1]:limitsSelVoxelsMatrix[passIdx,2] ]
    selTsVoxelMean <- apply( selTsVoxel, 2, mean ) #voxels average
    ssTot <- apply( (selTsVoxel-selTsVoxelMean)^2, 2, sum) #voxels total sum of squares
    runLinMod <- function( nIndex ) { #get best fit for every prediction (nIndex = counter of predictions)
      index_eyeMoving <- predictionIndexesMatrix[nIndex,1] #this is the index that links predictionIndexesMatrix with predictionGrid_eyeMovingStim
      index_eyeFixStim_border <- predictionIndexesMatrix[nIndex,2] #this is the index that links predictionIndexesMatrix with predictionGrid_eyeFixStim_border
      if (modelType==1) { # mixed model
        dMat <- cbind( tsPrediction_eyeMovingStim[index_eyeMoving,], tsPrediction_eyeFixStim_border[index_eyeFixStim_border,] )      
        dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
        a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (3: intercept and slope _eyeMovingStim and slope _eyeFixStim_border )
      }
      if (modelType==2) { # eyeFix model with border
        dMat <- cbind( tsPrediction_eyeFixStim_border[index_eyeFixStim_border,] )      
        #dMat <- cbind( tsPrediction_eyeMovingStim[index_eyeMoving,] )      
        dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
        a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope and slope _eyeFixStim_border )
      }
      if (modelType==3) { # eye moving model
        dMat <- cbind( tsPrediction_eyeMovingStim[index_eyeMoving,] )      
        #dMat <- cbind( tsPrediction_eyeMovingStim[index_eyeMoving,] )      
        dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
        a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope and slope _eyeMovingStim )
      }
      #a <- solve( qr(dMat01), selTsVoxel)
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot) #r squares
      if (modelType==2 | modelType==3) { a <- rbind( a, rep(999,dim(a)[2]) ) }
      eyeMovingArray <- rep( index_eyeMoving, length(r2) )
      eyeFixStimArray <- rep( index_eyeFixStim_border, length(r2) )
      return( rbind( r2, a, eyeMovingArray, eyeFixStimArray, expectedTs ) ) #fields: r2, 
      # 3 beta coeffs (intercept slope _eyeMovingStim and slope _eyeFixStim_border), 
      # index predictor eye moving (da predictionIndexesMatrix), 
      # index predictor eyeFixStim_border (da predictionIndexesMatrix), 
      # expected ts
    }
    outVoxelList <- lapply( 1:dim(predictionIndexesMatrix)[1], runLinMod  ) #apply the function for all predictions, get a list as output
    outVoxel3D <- abind( outVoxelList, along=3 ) #reshape the list in a 3dmatrix with fields: output X voxels tested X predictions (from predictionIndexesMatrix)
    betaPositiveMatrix <- outVoxel3D[3,,] > 0 & outVoxel3D[4,,] > 0 #fields slope _eyeMovingStim and slope _eyeFixStim_border must be positive, output = matrix with voxelsXpredictions (from predictionIndexesMatrix)
    r2Matrix <- outVoxel3D[1,,]	# extracts r2, output = matrix with voxelsXpredictions (from predictionIndexesMatrix)
    index_eyeMovingMatrix <- outVoxel3D[5,,] # output = matrix with voxelsXpredictions (from predictionIndexesMatrix)
    index_eyeFixStim_borderMatrix <- outVoxel3D[6,,] # output = matrix with voxelsXpredictions (from predictionIndexesMatrix)
    extractData <- function(nSelectedVoxels) { #function that, for each voxel tested, looks for positive betas and extracts max(r2) and relevant parameters
      indexBetaZero <- betaPositiveMatrix[nSelectedVoxels,]
      if ( sum(indexBetaZero)>0 ) {
        indexBetaZero <- which( indexBetaZero )
        indexVarExp <- which.max( r2Matrix[nSelectedVoxels,indexBetaZero] )
        eyeMovingPredictorIndex <- index_eyeMovingMatrix[nSelectedVoxels, indexBetaZero[indexVarExp] ] #to be found in predictionGrid_eyeMovingStim
        eyeFixBorderPredictorIndex <- index_eyeFixStim_borderMatrix[nSelectedVoxels, indexBetaZero[indexVarExp] ] #to be found in predictionGrid_eyeFixStim_border
        prediction_eyeMoving <- predictionGrid_eyeMovingStim[eyeMovingPredictorIndex, ]
        prediction_eyeFixBorder <- predictionGrid_eyeFixStim_border[eyeFixBorderPredictorIndex, ]
        return( as.numeric( c( prediction_eyeMoving, prediction_eyeFixBorder, outVoxel3D[,nSelectedVoxels,indexBetaZero[indexVarExp]] ) ) ) 
      }
      if ( sum(indexBetaZero)==0 ) {
        return( rep(0, dim(predictionGrid_eyeMovingStim)[2]+dim(predictionGrid_eyeFixStim_border)[2]+6+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 6=r2, beta intercept, beta slope eyeMovingStim, beta slope eyeFixStim_border, index predictor eye moving (da predictionIndexesMatrix), index predictor eye fix border (da predictionIndexesMatrix)
      } 
    }	
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }
  #system.time( aaa <- lapply( 1:3, voxelModel ) ) #to test out comment otherwise
  library(parallel)
  detectCores()
  nCores <- 4
  cl <- makeCluster(nCores, type='FORK')
  storeTimePar <- system.time( outModel <- parLapply(cl, runIndex, voxelModel ) )
  stopCluster(cl)
  print( storeTimePar )
  
  for ( nElements in 1:length(outModel) ) { #this concatenates al outputs (49 loops ran in parallel, normally, into 1 matrix with dimensions: outputModel X nVoxel (global))
    if (nElements==1) {
      outMatrix <- outModel[[nElements]]
    }
    if (nElements>1) {
      outMatrix <- cbind( outMatrix, outModel[[nElements]] )
    }
  }
  
  outModel <- outMatrix #this updates the results, taking r2 values (outModel[9,]) across the loops and taking the max
  if (modelFitCounter==1) { outModelLoop <- outModel }
  if (modelFitCounter>1) { 
    selectedCols <- outModel[9,] > outModelLoop[9,]
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
}


#outModelLoop <- outMatrix
storeAllPred <- t(outModelLoop[ seq(1,14), ]) #takes all the output except the expected time serie
#### extract surround size #### fix this bor both stimuli tested and combined
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
      r <- a - b
      
      maxIdx <- which.max(r)
      rHalf <- r[ maxIdx:length(xFWHM) ]
      xHalf <- xFWHM[ maxIdx:length(xFWHM) ]
      halfMax <- max(rHalf)/2
      xMid <- xHalf[ rHalf<=halfMax  ]  
      FWHMcenter <- min( xMid )*2 
      rMin <- which.min(rHalf)
      FWHMsurround <- xHalf[rMin[1]]*2
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
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 15:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(18, dim(tsTransposedAll)[2] ) )

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

labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','x-PosFix','y-PosFix','sigmaPosFix','sigmaNegFix','varExp','intercept','slopeEyeMoving','slopeEyeFix','indexEyeMoving','indexEyeFix','theta','radius','fwhmCenter','surroundSize')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

