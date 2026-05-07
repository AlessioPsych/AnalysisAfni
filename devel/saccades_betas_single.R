#args <- commandArgs(T)
#print( args )

#models each sacc direction in order to take the betas and fit a model

rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

args <- c('greyMask.nii.gz', 'meanTs_eye_topUp_res.nii', 'sacc_test_single', '0.166','1')

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
hrfDelayOnsetElements <- 4
hrfDelayUnderShootElements <- 4
hrfb1Elements <- 4
hrfb2Elements <- 4
hrfcElements <- 4
hrf_a1 <- seq( 6, 10, length.out=hrfDelayOnsetElements )
hrf_a2 <- seq( 12, 16, length.out=hrfDelayUnderShootElements )
hrf_b1 <- seq( 0.6, 1.1, length.out=hrfb1Elements  ) 
hrf_b2 <- seq( 0.6, 1.1, length.out=hrfb2Elements  ) 
hrf_c <- seq( 0.1, 0.6, length.out=hrfcElements  ) 

predictionGridGlobal <- expand.grid( hrf_a1, hrf_a2, hrf_b1, hrf_b2, hrf_c )
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
  
  print( sprintf( 'iteration %1.0f of %1.0f, start...', modelFitCounter, length(runIndexPredictions)  ) )
  
  # here I select the portion of all the predictions that are going to be tested later
  predictionGrid <- predictionGridGlobal[ limitsPredictionMatrix[modelFitCounter,1]:limitsPredictionMatrix[modelFitCounter,2], ]
  print( sprintf( 'generate predictions in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  
  ### eye direction predictions ####
  generatePrediction_direction <- function(indexPrediction, inputPredictionGrid, returnFlag) {
    
    #indexPrediction <- 15
    #inputPredictionGrid <- predictionGrid
    
    #predPar <- rep(0,12)
    #prfPar[c(5,6,8,9,10,11)] <- c(6,12,-pi*0.9,6,0.5,0.8)
    #plot( abs(directionPol[,1]*baselines)>0, type='l' )
    
    predPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=predPar[1], a2=predPar[2], b1=predPar[3], b2=predPar[4], c=predPar[5]), verbose=FALSE )
    
    directionFactor <- as.factor( round( polarCoords[,1]*baselines, 3 ) )
    rawPredictors <- array( 0, c( dim(directionPol)[1], length( levels( directionFactor ) )-1 ) )
    predictorCounter <- 1
    for ( dirCounter in 1:length(levels( directionFactor )) ) {
      if ( levels( directionFactor )[ dirCounter ] != '0' ) {
        rawPredictors[,predictorCounter] <- directionFactor==levels( directionFactor )[ dirCounter ]
        predictorCounter <- predictorCounter+1
      }
    }
    polar( polarCoords[,1]*baselines, polarCoords[,2]*baselines, type='p'  )
    
    convPredictors <- array( 0, c( dim(ts$brk)[4], length( levels( directionFactor ) )-1 ) )
    for ( predCount in 1:dim(rawPredictors)[2] ) {
      convTs <- conv( rawPredictors[,predCount], hrf )
      convTsTrim <- convTs[ 1:length( rawPredictors[,predCount] ) ]
      timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
      iTs <- interp1( timeArray, convTsTrim, seq(0,max(timeArray),length.out = dim(ts$brk)[4] ) )
      scaleIts <- scaleData( round( iTs, 4 ), 1, 0 )
      convPredictors[,predCount] <- scaleIts
    }

    return( convPredictors )
    
  }
  
  print( sprintf( 'generating predictions in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  system.time( eyePredictions_direction <- lapply( 1:dim(predictionGrid)[1], generatePrediction_direction, inputPredictionGrid=predictionGrid ) )
  
  # #### clean up ts predictions and prediction grid ####
  # tsPrediction <- t( eyePredictions_direction )
  # tsPrediction_saccOn <- t( eyePredictions_saccOn )
  # controlPredictions <- apply( tsPrediction, 1, sum ) != 0 & apply( tsPrediction_saccOn, 1, sum ) != 0
  # tsPrediction_direction <- tsPrediction[controlPredictions,]
  # tsPrediction_saccOn <- tsPrediction_saccOn[controlPredictions,]
  # predictionGrid <- predictionGrid[controlPredictions,]
  # print( dim( predictionGrid ) )
  
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
      
      dMat <- cbind( eyePredictions_direction[[ nIndex ]] )  #visual, movement and direction predictors, multiplied to avoid collinearity issues    
      
      dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (4: intercept, beta eyeFix, beta eyeMoving, beta direction
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot) #r squares
      
      return( rbind( r2, a, expectedTs ) )
      
    }
    outVoxelList <- lapply( 1:length(eyePredictions_direction), runLinMod  ) #apply the function for all predictions, get a list as output
    outVoxel3D <- abind( outVoxelList, along=3 ) #reshape the list in a 3dmatrix with fields: output X voxels tested X predictions
    
    #betaPositiveMatrix <- outVoxel3D[3,,] > 0 & outVoxel3D[4,,] > 0
    
    r2Matrix <- outVoxel3D[1,,]	
    extractData <- function(nSelectedVoxels) {
      indexVarExp <- which.max( r2Matrix[nSelectedVoxels,] )
      if ( sum(indexVarExp)>0 ) {
        return( as.numeric( c( predictionGrid[indexVarExp,], outVoxel3D[,nSelectedVoxels,indexVarExp] ) ) ) 
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
    selectedCols <- outModel[dim(predictionGrid)[2]+1,] > outModelLoop[dim(predictionGrid)[2]+1,] #where r2 is store (after 13 model parameters, there are: r2, beta intercept, beta slope visual, beta slope movement, beta slope direction 
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
} 
storeAllPred <- cbind( t(outModelLoop[ seq(1,17), ]) ) # 9 parameters + R2 + intercept and slope1 and slope 2

## voxel tuning model
print( sprintf( 'voxel tuning model ...', modelFitCounter, length(runIndexPredictions)  ) )
#xRad <- seq(-pi,pi,length.out = 12)
#kPar <- seq(0.05,6,length.out = 14)
#voxTuningDesign <- expand.grid(xRad,kPar)
meanMod <- seq(-pi,pi,length.out = 12)
sdMod <- seq(0.25,4,length.out = 12)
voxTuningDesign <- expand.grid(meanMod,sdMod)
whichHighR2 <- which( storeAllPred[,6] > 0.25 )
directionFactor <- as.factor( round( polarCoords[,1]*baselines, 3 ) )
directionFactorLevels <- as.numeric( as.character( levels( directionFactor ) ) )
directionFactorLevelsNoZero <- directionFactorLevels[ abs( directionFactorLevels )>0.0001 ]
voxelTuning <- function( idx ) {
  voxelPredictions <- function( idxPrediction, selVoxelBetas ) {
    #idxPrediction <- 50
    #vDist <- vonMisesDist( directionFactorLevelsNoZero, voxTuningDesign[idxPrediction,1], voxTuningDesign[idxPrediction,2], 1 )
    vDist <- scaleData( dnorm( directionFactorLevelsNoZero, voxTuningDesign[idxPrediction,1], voxTuningDesign[idxPrediction,2] ), 1, 0 )
    dMat01 <- cbind( rep(1,length(vDist)), vDist ) #column of ones and column of predictor
    a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selVoxelBetas) ) #beta coefficients (2: intercept, beta )
    expectedBetas <- crossprod( t(dMat01), a ) #expected ts
    residualsSquared <- ( selVoxelBetas - expectedBetas )^2 
    ssRes <- sum( residualsSquared )
    ssTot <- sum( ( selVoxelBetas-mean( selVoxelBetas ) )^2 )
    r2 <- 1-(ssRes/ssTot) #r squared
    return( c(r2,a,expectedBetas,selVoxelBetas) )
  }
  idxVoxel <- whichHighR2[idx]
  selBetas <- storeAllPred[ idxVoxel, 8:dim(storeAllPred)[2] ]
  voxelModelOutput <- sapply( 1:dim(voxTuningDesign)[1], voxelPredictions, selVoxelBetas=selBetas )
  r2Output <- voxelModelOutput[1,]
  betaOutput <- voxelModelOutput[3,]
  fitSelected <- which( betaOutput > 0 )
  if (sum(fitSelected)>0) {
    r2max <- which.max( r2Output[ fitSelected ] )
    selectedR2 <- fitSelected[ r2max ]
    winningmodelParameters <- as.numeric( voxTuningDesign[ selectedR2, ] )
    winningExpectedBetas <- voxelModelOutput[ 4:13, selectedR2  ]
    winningObservedBetas <- voxelModelOutput[ 14:dim(voxelModelOutput)[1], selectedR2  ]
    winningR2 <- round( r2Output[ selectedR2 ], 3 )
    #vDist <- vonMisesDist( directionFactorLevelsNoZero, winningmodelParameters[1], winningmodelParameters[2], 1 )
    #plot( storeAllPred[ idxVoxel, 8:dim(storeAllPred)[2] ], col='black', pch=16 )
    #lines( winningExpectedBetas, col='red' )
    #plot(vDist)
    #selBetas==winningObservedBetas
  }
  if (sum(fitSelected)==0) {
    winningmodelParameters <- rep(0,2)
    winningExpectedBetas <- rep(0,10)
    winningObservedBetas <- rep(0,10)
    winningR2 <- 0
  }
  return( c( winningR2, winningmodelParameters, winningExpectedBetas, winningObservedBetas ) )
}
detectCores()
nCores <- 4
cl <- makeCluster(nCores, type='FORK')
runIndex <- 1:length( whichHighR2 )
storeTimePar <- system.time( outTuning <- parSapply(cl, runIndex, voxelTuning ) )
stopCluster(cl)
print( storeTimePar )


## save output
storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 18:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(17, dim(tsTransposedAll)[2] ) )
storeAllPredOut[,selIdxVoxel] <- t( storeAllPred )
storeAllPredOutModel <- array( 0, c(3, dim(tsTransposedAll)[2] ) )
storeAllPredOutModel[ , selIdxVoxel[whichHighR2] ] <- outTuning[1:3,]
storeAllPredOutExpectedBetas <- array( 0, c( 10, dim(tsTransposedAll)[2] ) )
storeAllPredOutExpectedBetas[ , selIdxVoxel[whichHighR2] ] <- outTuning[4:13,]

###
print('save linear step...')
fileTs <- sprintf('%s_PredixtedTs.nii.gz',outSuffix) 
fileParams <- sprintf('%s_params.nii.gz',outSuffix)
fileModelParams <- sprintf('%s_modelParams.nii.gz',outSuffix) 
fileModelPredictedBetas <- sprintf('%s_PredixtedBetas.nii.gz',outSuffix)

###
rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( '__tt_expectedTs.nii.gz', brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_expectedTs.nii.gz', orientValue[1,1], fileTs )
system( instr)
system( 'rm __tt_expectedTs.nii.gz' )

###
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
labels <- c('hrf_a1','hrf_a2','hrf_b1','hrf_b1','hrf_c','varExp')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

###
instr <- sprintf( '3dTcat -prefix %s_observedBetas.nii.gz %s_params.nii.gz[7..16]', outSuffix, outSuffix )
system( instr )

###
rParameters <- array( storeAllPredOutModel, c( dim(storeAllPredOutModel)[1], dim(ts$brk)[1:3] ) )
rParameters <- aperm( rParameters, c(2,3,4,1) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', epiFile )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( '__tt_parameters.nii.gz', brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_parameters.nii.gz', orientValue[1,1], fileModelParams )
system( instr)
system( 'rm __tt_parameters.nii.gz' )
labels <- c('varExp','centerVon','kVon')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileModelParams )
  print( instr )
  system( instr )
}

###
rParameters <- array( storeAllPredOutExpectedBetas, c( dim(storeAllPredOutExpectedBetas)[1], dim(ts$brk)[1:3] ) )
rParameters <- aperm( rParameters, c(2,3,4,1) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', epiFile )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( '__tt_parameters.nii.gz', brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_parameters.nii.gz', orientValue[1,1], fileModelPredictedBetas )
system( instr)
system( 'rm __tt_parameters.nii.gz' )
