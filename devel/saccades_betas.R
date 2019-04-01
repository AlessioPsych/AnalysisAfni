#args <- commandArgs(T)
#print( args )

#models eye direction and saccade onset, no visual stimuli

rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

args <- c('greyMask.nii.gz', 'meanTs_eye_topUp_res.nii', 'sacc_test_08_pos', '0.166','1')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir ) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
source( sprintf('%s/vonMisesDist.R', generalPurposeDir ) )

library( pracma )
library( abind )
library( neuRosim )
library( parallel )
library( matrixcalc )

epiFile <- args[1]
inputFile <- args[2]
outSuffix <- args[3]
samplingTime <- as.numeric( args[4] )
polortArg <- as.numeric( args[5] )

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

# load design eye fix to get the baseline
setwd('/analyse/Project0226/dataSummary')
eyeDesignFile <- read.table( file='eyeMovingDesign.txt', header = FALSE, as.is=TRUE )
eyeDesignBaseline <- ifelse( eyeDesignFile[,3]>15, 0, 1 )
baselines <- eyeDesignBaseline
directionPol <- eyeDesignFile[,c(1,2)]
directionPol[,1] <- deg2rad( round( directionPol[,1] , 3 ) )
directionArrayCart <- pol2cart( cbind( directionPol[,1], directionPol[,2] ) ) # this is necessary to convert from matlab to R polar coords
directionArrayCart[,1] <- directionArrayCart[,1]
directionArrayCart[,2] <- directionArrayCart[,2]*-1

par(mfrow=c(2,2)) # plot first target (up-right) and all targets wrt fixation
plotidx <- 1:100
plot( directionArrayCart[plotidx,1]*baselines[plotidx], directionArrayCart[plotidx,2]*baselines[plotidx], xlim=c(-3,3), ylim=c(-2,2) ) 
plot( directionArrayCart[,1]*baselines, directionArrayCart[,2]*baselines*-1, xlim=c(-3,3), ylim=c(-2,2)  ) 
polarCoords <- cart2pol( cbind( directionArrayCart[,1], directionArrayCart[,2] ) )
polar( polarCoords[plotidx,1]*baselines[plotidx], polarCoords[plotidx,2]*baselines[plotidx], type='p'  )
polar( polarCoords[,1]*baselines, polarCoords[,2]*baselines, type='p'  )

stimuliSequenceTr <- samplingTime
setwd(mainDir)

# design
#hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=6, a2=10, b1=0.9, b2=0.9, c=0.5), verbose=FALSE )
#plot( hrf~seq(0,30,samplingTime) )
hrfDelayOnsetElements <- 2
hrfDelayUnderShootElements <- 2
hrfb1Elements <- 3
hrfb2Elements <- 3
hrfcElements <- 3
meanElements <- 10
sdElements <- 10
hrf_a1 <- seq( 6, 10, length.out=hrfDelayOnsetElements )
hrf_a2 <- seq( 12, 16, length.out=hrfDelayUnderShootElements )
hrf_b1 <- seq( 0.7, 1, length.out=hrfb1Elements  ) 
hrf_b2 <- seq( 0.7, 1, length.out=hrfb2Elements  ) 
hrf_c <- seq( 0.1, 0.6, length.out=hrfcElements  ) 
centerTuning <- seq( -pi*0.99, pi*0.99, length.out=meanElements ) 
sdTuning <- seq( 0.1, 3, length.out=sdElements ) 

predictionGridGlobal <- expand.grid( hrf_a1, hrf_a2, hrf_b1, hrf_b2, hrf_c, centerTuning, sdTuning )
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

for (modelFitCounter in 1:length(runIndexPredictions) ) { #length(runIndexPredictions) #modelFitCounter <- 1
  
  ## test vonMises distribution ##
  #par(mfrow=c(1,3))
  #x <- seq(-pi,pi,0.1)
  #mu <- pi*0.1
  #k <- 2
  #vonValues <- vonMisesDist( x, mu, k, 3, 0.6, 0.1 )
  #plot( vonValues~x )
  
  print( sprintf( 'iteration %1.0f of %1.0f, start...', modelFitCounter, length(runIndexPredictions)  ) )
  
  # here I select the portion of all the predictions that are going to be tested later
  predictionGrid <- predictionGridGlobal[ limitsPredictionMatrix[modelFitCounter,1]:limitsPredictionMatrix[modelFitCounter,2], ]
  print( sprintf( 'generate predictions in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  
  ### eye direction predictions ####
  generatePrediction_direction <- function(indexPrediction, inputPredictionGrid) {

    #indexPrediction <- 15
    #inputPredictionGrid <- predictionGrid

    #predPar <- rep(0,12)
    #prfPar[c(5,6,8,9,10,11)] <- c(6,12,-pi*0.9,6,0.5,0.8)
    #plot( abs(directionPol[,1]*baselines)>0, type='l' )
    
    predPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    #predPar <- c( 6, 12, 0.9, 0.9, 0.35, 2.8, 4)
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=predPar[1], a2=predPar[2], b1=predPar[3], b2=predPar[4], c=predPar[5]), verbose=FALSE )

    #polar( polarCoords[,1]*baselines, polarCoords[,2]*baselines, type='p'  )
    
    muRad <- predPar[6]
    kappaPar <- predPar[7]
    #sdRadNeg <- predPar[8]
    directionArray <- polarCoords[,1]
    #expectedTs01 <- dnorm( directionArray, meanRad, sdRad ) #plot( dnorm( seq( -pi, pi, 0.05 ), meanRad, sdRad ) )
    #expectedTs01 <- sin( directionArray + meanRad ) #plot( dnorm( seq( -pi, pi, 0.05 ), meanRad, sdRad ) )
    #expectedTs01 <- dnorm( directionArray, meanRad, sdRad )-dnorm( directionArray, meanRad, sdRadNeg ) #plot( dnorm( seq( -pi, pi, 0.05 ), meanRad, sdRad )-dnorm( seq( -pi, pi, 0.05 ), meanRad, sdRadNeg ) )
    expectedTs01 <- suppressWarnings( vonMisesDist( directionArray, muRad, kappaPar,  1 ) ) # direction read throught the current VonMises distribution
    expectedTs02 <- expectedTs01*baselines  #code saccades according to the vonMises distribution, every 200 ms (stimuliSequenceTr)

    convTs <- conv( expectedTs02, hrf )
    convTsTrim <- convTs[ 1:length( expectedTs02 ) ]
    timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
    iTs <- interp1( timeArray, convTsTrim, seq(0,max(timeArray),length.out = dim(ts$brk)[4] ) )
    scaleConvTsTrim <- scaleData( round( convTsTrim, 5 ), 1, 0 )
    scaleIts <- scaleData( round( iTs, 5 ), 1, 0 )
    return( scaleIts )
    
    #convTs_saccOn <- conv( baselines, hrf )
    #convTsTrim_saccOn <- convTs_saccOn[ 1:length( firstStepTs_orig_baseline ) ]
    #timeArray <- seq( 0, length(convTsTrim_saccOn)*stimuliSequenceTr, length.out = length(convTsTrim_saccOn) )#  *stimuliSequenceTr/1000
    #iTs_saccOn <- interp1( timeArray, convTsTrim_saccOn, seq(0,max(timeArray),length.out = dim(ts$brk)[4] ) )
    #scaleConvTsTrim_saccOn <- scaleData( round( convTsTrim_saccOn, 5 ), 1, 0 )
    #scaleIts_saccOn <- scaleData( round( iTs_saccOn, 5 ), 1, 0 )
    
    #if (returnFlag==1) { return( scaleIts ) }
    #if (returnFlag==0) { return( scaleIts_saccOn ) }

    ##to plot
    #par(mfrow=c(3,3))
    #piSpace <- seq(-pi*0.95, +pi*0.95, length.out = 200)
    #plot( baselines, type='l' )
    #plot( expectedTs01, type='l' )
    #plot( expectedTs02, type='l' )
    #plot( dnorm( piSpace, meanRad, sdRad ) ~ piSpace )
    #plot( sin( piSpace+meanRad ) ~ piSpace )
    #plot( convTsTrim~timeArray )
    #lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), iTs, col='red' )
    #plot( scaleConvTsTrim~timeArray )
    #lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), scaleIts, col='red' )
    #plot( scaleConvTsTrim_saccOn~timeArray )
    #lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), scaleIts_saccOn, col='red' )
    
  }
  eyePredictions_direction <- sapply( 1:dim(predictionGrid)[1], generatePrediction_direction, inputPredictionGrid=predictionGrid )
  #eyePredictions_saccOn <- sapply( 1:dim(predictionGrid)[1], generatePrediction_direction, inputPredictionGrid=predictionGrid, returnFlag=0 )
  

  #### clean up ts predictions and prediction grid ####
  tsPrediction <- t( eyePredictions_direction )
  controlPredictions <- apply( tsPrediction, 1, sum ) != 0 & !is.na( apply( tsPrediction, 1, sum ) )
  #tsPrediction_saccOn <- t( eyePredictions_saccOn )
  #controlPredictions <- apply( tsPrediction, 1, sum ) != 0 & apply( tsPrediction_saccOn, 1, sum ) != 0
  tsPrediction_direction <- tsPrediction[controlPredictions,]
  #tsPrediction_saccOn <- tsPrediction_saccOn[controlPredictions,]
  predictionGrid <- predictionGrid[controlPredictions,]
  print( dim( predictionGrid ) )
  
  #### linear fitting in parallel ####
  print( sprintf( 'linear fitting in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  
  if (length(selIdxVoxel) <= 200 ) {
    totalVoxelsIterations <- 3
  }
  if (length(selIdxVoxel) > 200 ) {
    totalVoxelsIterations <- ceil( length(selIdxVoxel)/100 )
  }
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
      
      dMat <- cbind( tsPrediction_direction[nIndex,] )  #visual, movement and direction predictors, multiplied to avoid collinearity issues    
      
      dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (4: intercept, beta eyeFix, beta eyeMoving, beta direction
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot) #r squares
      
      return( rbind( r2, a, expectedTs ) )
      
    }
    outVoxelList <- lapply( 1:dim(tsPrediction_direction)[1], runLinMod  ) #apply the function for all predictions, get a list as output
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
          return( rep(0, dim(predictionGrid)[2]+1+2+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 1+3=r2,+intercept, slope direction
      }
    }	
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }
  #system.time( aaa <- lapply( 1:20, voxelModel ) )
  detectCores()
  nCores <- 4
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
    selectedCols <- outModel[dim(predictionGrid)[2]+1,] > outModelLoop[dim(predictionGrid)[2]+1,] #where r2 is store (after model parameters, there are: r2, beta intercept, beta slope direction 
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
} #length(runIndexPredictions)


hist( outModelLoop[8, outModelLoop[8,]>0 ] )
# # get single parameter contribution
# print('get single parameter contribution...')
# fittedPredictionGrid <- t( outModelLoop[seq(1:12),] )
# singlePredictorStats <- function( idxLine, inputMat ) {
#   
#   #idxLine <- 57
#   #inputMat <- fittedPredictionGrid
#   
#   if ( sum(inputMat[idxLine,])!=0 ) {
#     
#     directionPrediction <- eyeDirectionPrediction(idxLine,inputMat)
#     gainPrediction <- generatePrediction_moving(idxLine,inputMat)
#     eyeCenterPrediction <- generatePrediction_visual(idxLine,inputMat)
#     tsToFit <- tsTransposedSel[,idxLine]
#     
#     ssTotSingleVoxel <- sum ( (tsToFit-mean(tsToFit))^2 )  #voxel total sum of squares
#     dMat <- cbind( directionPrediction, gainPrediction, eyeCenterPrediction )  #visual, movement and direction predictors    
#     dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
#     a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,tsToFit) ) #beta coefficients (4: intercept, beta eyeFix, beta eyeMoving, beta direction
#     expectedTs <- crossprod( t(dMat01), a ) #expected ts
#     residualsSquared <- ( tsToFit - expectedTs )^2 
#     ssRes <- sum( residualsSquared )
#     r2 <- 1-(ssRes/ssTotSingleVoxel) #r squares
#     
#     dfTot <- dim(dMat01)[1]-dim(dMat01)[2]
#     seCoeff <- sqrt( diag( ssRes/dfTot*matrix.power( ( t(dMat01) %*% dMat01 ) , -1 ) ) ) # standard error based without lm
#     tStats <- a/seCoeff
#     pVals <- 2 * pt(-abs(tStats), dfTot)
#   }
#   
#   if ( sum(inputMat[idxLine,])==0 ) {
#     r2 <- 0
#     tStats <- rep(0,4)
#     pVals <- rep(0,4)
#   }
#   
#   return( c(r2, tStats, pVals) )
#   
#   ## to test
#   # nVox <- 1
#   # r2[nVox]
#   # rbind( coef( lm( selTsVoxel[,nVox] ~ dMat ) ), a[,nVox] ) #these two lines should be VERY similar
#   # fit = lm( selTsVoxel[,nVox] ~ dMat )
#   # se <- sqrt(diag(vcov(fit)))
#   # se # standard error based on the lm function in r, slow
#   # seCoeff[nVox] # standard error based on the analytical formula
#   # # S^2{b} = MSE * (X^T * X)^-1 # the analytical formula for standard error in multiple linear regression
#   # # t stats
#   # summary(fit)$coefficients[,3]
#   # tStats
#   # # pvals
#   # summary(fit)$coefficients[,4]
#   # pVals
#   
# }
# runIndex <- 1:dim(fittedPredictionGrid)[1]
# detectCores()
# nCores <- 4
# cl <- makeCluster(nCores, type='FORK')
# storeTimePar <- system.time( outModelStats <- parSapply(cl, runIndex, singlePredictorStats, inputMat=fittedPredictionGrid ) )
# stopCluster(cl)
# print( storeTimePar )
# 
# storeAllPred <- cbind( t(outModelLoop[ seq(1,15), ]), t( outModelStats ) )

# #### extract surround size ####
# print('surround size...')
# if (flagSurround==1) { 
#   print('get FWHM...')
#   progress <- 0.05
#   FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
#   xFWHM <- seq(-25,25,0.0125)
#   for (k in 1:dim(storeAllPred)[1]) {
#     
#     parameters <- storeAllPred[k,]
#     
#     if ( sum(parameters!=0) ) {
#       a <- dnorm( xFWHM, 0, parameters[3] )
#       b <- dnorm( xFWHM, 0, parameters[4] )
#       r <- scaleData(a,1,0) - scaleData(b,1,0)*parameters[7]
#       
#       maxIdx <- which.max(r)
#       rHalf <- r[ maxIdx:length(xFWHM) ]
#       xHalf <- xFWHM[ maxIdx:length(xFWHM) ]
#       halfMax <- max(rHalf)/2
#       xMid <- xHalf[ rHalf<=halfMax  ]  
#       FWHMcenter <- min( xMid )*2 
#       if ( sum(r<0)>0 ) { #if there is a detectable surround, compute surround size
#         rMin <- which.min(rHalf)
#         FWHMsurround <- xHalf[rMin[1]]*2
#       }
#       if ( sum(r<0)==0 ) { #otherwise set it to zero
#         FWHMsurround <- 0
#       }
#       FWHM[k,] <- c(FWHMcenter,FWHMsurround)
#     }
#     
#     if (k>dim(storeAllPred)[1]*progress) {
#       cat(paste('*',""))
#       progress <- progress + 0.05
#     }
#   }
# }
# if (flagSurround==0) { 
#   FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
#   FWHM[,1] <- storeAllPred[,3]
#   FWHM[,2] <- 1000
# }

#storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
#storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 18:dim(outModelLoop)[1], ]	
#storeAllPredOut <- array( 0, c(30, dim(tsTransposedAll)[2] ) )

storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 11:dim(outModelLoop)[1], ]	
storeAllPred <- cbind( t(outModelLoop[ seq(1,10), ]) ) # 7 parameters + R2 + intercept and slope direction
storeAllPredOut <- array( 0, c(10, dim(tsTransposedAll)[2] ) )
storeAllPredOut[,selIdxVoxel] <- t( storeAllPred )

print('save linear step...')
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

# labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multPar','centerVonMis','kVonMis','xGain','yGain','sigmaGain',
#             'varExp','intercept','b_eyeFix','b_eyeMoving','b_direction','varExp_cont','t_intercept','t_eyeFix','t_eyeMoving','t_direction',
#             'p_intercept','p_eyeFix','p_eyeMoving','p_direction','theta','radius','FWHM_center','surround_size')
# for (k in 1:length(labels)) {
#   instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
#   print( instr )
#   system( instr )
# }

labels <- c('hrf_a1','hrf_a2','hrf_b1','hrf_b1','hrf_c','meanPi','sdPi',
            'varExp','intercept','b_birection','b_saccOn')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

#xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, centerVon, kVon, saccSuppression, xPosFit_gain, yPosFit_gain, sigmaArrayPositive_gain )

