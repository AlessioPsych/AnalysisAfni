

args <- commandArgs(T)
print( args )

#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#args <- c('meanTs_bars_res_mask.nii','meanTs_bars_res.nii','outTest_delme.nii.gz','/packages/afni/17.0.13','1')

## prepare the stimuli with the portrait as well

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')


#args <- c('meanTs_bars_res_mask.nii', 'meanTs_bars_res.nii', 'output_test_delme_bars','/home/alessiof/abin', '1', '1','1','0.166','3')
#args <- c('meanTs_eyeMovement_topUp_res_mask.nii', 'meanTs_eyeMovement_topUp_res.nii', 'output_test_delme','/home/alessiof/abin', '1', '1','1','0.166','1','border_fit_parallel_nohrf_params.nii.gz','border_fit_parallel_nohrf_PredixtedTs.nii.gz')

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
samplingTime <- 0.166#as.numeric( args[8] )
stimType <- as.numeric( args[8] )
fineFit <- 0
paramsFile <- args[9]
predictedTsFile <- args[10]

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

# already fitted model
tsPredicted <- read.AFNI( predictedTsFile )
tsArrayPredicted <- array( tsPredicted$brk, c( prod( dim( tsPredicted$brk )[1:3] ), dim(tsPredicted$brk)[4] ) )
paramsPredicted <- read.AFNI( paramsFile )
paramsArrayPredicted <- array( paramsPredicted$brk, c( prod( dim( paramsPredicted$brk )[1:3] ), dim(paramsPredicted$brk)[4] ) )

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
  stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==2) { 
  arrayStim <- scan( 'eyeFixStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==3) { 
  arrayStim <- scan( 'eyeFixStim_border.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==4) { 
  arrayStim <- scan( 'prfStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(200,1860,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==5) { 
  arrayStim <- scan( 'eyeFixStim_border_disappear.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==6) { 
  arrayStim <- scan( 'eyeMovingStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
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
if (stimType==6) { 
  stimMatFlip <- stimMatFlip[ dim(stimMatFlip)[1]:1,dim(stimMatFlip)[2]:1, ]
}

x11( width=3, height=3 )
for ( snap in 1:dim(stimMat)[3] ) {
  image( stimMatFlip[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01) 
}
stimSeq <- stimMatFlip
x <- seq(-8,8,length.out = dim(stimSeq)[1] )
y <- seq(-6,6,length.out = dim(stimSeq)[2] )

#### get the standard hrf ####
#print('get hrf...')
#hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=9, a2=15, b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
#plot(hrf~seq(0,30,samplingTime))

#this part of the code builds a matrix with all the possible prediction tested, for both models at this stage
addSpace <- abs( min(x) )*0.1
print('build prediction...')
xPosFit <- seq( min(x)+addSpace, max(x)-addSpace, length.out=6 )
yPosFit <- seq( min(y)+addSpace, max(y)-addSpace, length.out=6 )
sigmaArrayPositive <- seq( 0.25, 6, length.out=6 )
if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative )
keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]

# this part builds a matrix with the starting and ending prediction index to fit for each loop in the next for loop
limitsPrediction <- round( seq(1,dim(predictionGridGlobal)[1], length.out=50) ) #split between 100 (arbitrary number) chuncks for iterative fitting
limitsPrediction[1] <- 1
limitsPrediction[ length(limitsPrediction) ] <- dim(predictionGridGlobal)[1]
limitsPredictionMatrix <- array( 0, c( (length(limitsPrediction)-1) , 2 ) )
limitsPredictionMatrix[,1] <- limitsPrediction[1:(length(limitsPrediction)-1)]
limitsPredictionMatrix[,2] <- limitsPrediction[2:(length(limitsPrediction))]
limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] <- limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
runIndexPredictions <- seq( 1:dim(limitsPredictionMatrix)[1] )

#fit visual respone and hrf
#modelFitCounter <- 1
for (modelFitCounter in 1:length(runIndexPredictions)) {
  
  print( sprintf( 'iteration %1.0f of %1.0f, start...', modelFitCounter, length(runIndexPredictions)  ) )
  # here I select the portion of all the predictions that are going to be tested later
  # I need to work on a voxel by voxel basis
  predictionGrid <- predictionGridGlobal[ limitsPredictionMatrix[modelFitCounter,1]:limitsPredictionMatrix[modelFitCounter,2], ]

  # define stimuli and timing
  stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
  incrementalCounter <- dim(predictionGrid)[1] * 0.05
  timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
  mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(ts$brk)[4] )
  
  ## linear fitting in parallel ##
  print( sprintf( 'linear fitting in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  indexVol <- meanEpi$brk[,,,1]
  indexArray <- array( indexVol, prod( dim( indexVol ) ) )
  r2ArraySelection <- paramsArrayPredicted[,7]
  hrf_a1Selection <- paramsArrayPredicted[,5]
  hrf_a2Selection <- paramsArrayPredicted[,6]
  tsTransposedAll <- t( tsArray )
  predictionTransposedAll <- t( tsArrayPredicted )
  selIdxVoxel <- which( indexArray == 1 ) # selection based on intensity and variance explained
  tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series
  tsPredictionTransposedSel <- predictionTransposedAll[,selIdxVoxel] #all selected predictions
  a1Selected <- hrf_a1Selection[ selIdxVoxel ]
  a2Selected <- hrf_a2Selection[ selIdxVoxel ]
  
  print( sprintf( 'generate predictions in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  ### function to generate the predictions for each single voxels ###
  generatePrediction <- function( indexPrediction, inputPredictionGrid, a1_par, a2_par ) {
    prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=a1_par, a2=a2_par, b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
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
  #generate all the predictions for a single voxels (needs to be voxel specific for the hrf parameters)
  #generatePrediction( predictionGrid, a1Selected[1], a2Selected[1] )
  # a1_par and a2_par are now fixed (6 and 12, respectively, for timing purposes)
  singleVoxel_predTs <- sapply( 1:dim(predictionGrid)[1], generatePrediction, inputPredictionGrid=predictionGrid, a1_par=a1Selected[1], a2_par=a2Selected[1] )
  
  voxelModel <- function( voxelIdx ) { #this fits the model on each selected voxel at a time
    selTsVoxel <- tsTransposedSel[,voxelIdx]
    selTsVoxelMean <- mean( selTsVoxel ) #voxel average
    ssTot <- sum( (selTsVoxel-selTsVoxelMean)^2 ) #voxels total sum of squares
    visual_predicted <- tsPredictionTransposedSel[,voxelIdx]
    #run the following to visualize the visual prediction, for testing purposes:
    #plot( selTsVoxel ); lines( visual_predicted, col='red' ); summary(lm(selTsVoxel~visual_predicted))$r.squared; r2ArraySelection[selIdxVoxel[voxelIdx]] #these last two numbers  must be (almost) idetical
    #summary(lm(selTsVoxel~visual_predicted)) #why is the parameter always 1?? because I have already estimated it in the previous fit, and this represents the best predicted ts, so the estimate will be very very close to 1 
    
    # function to run the model with the visual prediction + the gain field pred
    #to test: nIndex <- 1; visualPrediction <- visual_predicted
    runLinMod <- function(nIndex, visualPrediction) { #get best fit for every prediction (nIndex = counter of predictions)
      dMat <- cbind( singleVoxel_predTs[,nIndex], visualPrediction )      
      dMat01 <- cbind( rep(1,dim(dMat)[1]), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (3: intercept and slope gain, slope visualPred)
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      #run the following to visualize the visual prediction, and the new visual + gain prediction for testing purposes:
      #plot( selTsVoxel ); lines( visual_predicted, col='red' ); lines( expectedTs, col='blue' );
      
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- sum( residualsSquared )
      r2 <- 1-(ssRes/ssTot) #r squares
      return( rbind( r2, a, expectedTs ) )
    }
    
    # run the model on the singlevoxels for each gain field prediction
    singleVoxel_output <- sapply( 1:dim(singleVoxel_predTs)[2], runLinMod, visualPrediction=visual_predicted )

    #select positive betas for both visual and gain field predictors
    selectedBetas <- singleVoxel_output[3,] > 0 & singleVoxel_output[4,] > 0
    
    #select highest var exp:
    selectedBetasIndex <- which( selectedBetas )
    indexVarExp <- which.max( singleVoxel_output[1,selectedBetasIndex]  )
    indexSelectedPrediction <- selectedBetasIndex[indexVarExp]
    
    if ( sum( selectedBetas )>0 ) {
      newTs <- singleVoxel_output[5:dim(singleVoxel_output)[1],indexSelectedPrediction]
      #modVis <- lm( selTsVoxel~visual_predicted )
      #modVisGain <- lm( selTsVoxel~visual_predicted+newTs )
      #modTest <- anova(modVis,modVisGain)
      modTest <- list()
      modTest$F[2] <- 0
      modTest$`Pr(>F)`[2] <- 0
      singleVoxel_selectedOutput <- c( as.numeric( predictionGrid[indexSelectedPrediction,] ),modTest$F[2],modTest$`Pr(>F)`[2], singleVoxel_output[,indexSelectedPrediction] ) #array: 4 prediction, f test comparison, p value comparison,  r2, 3 beta coefficients (intercept and slope gain, slope visualPred) and the time series itself
    }
    if ( sum( selectedBetas )==0 ) {
      singleVoxel_selectedOutput <- rep( 0, dim(predictionGrid)[2]+2+dim(singleVoxel_output)[1] )
    }
    #run the following to visualize the visual prediction, and the best visual + gain prediction for testing purposes:
    #plot( selTsVoxel ); lines( visual_predicted, col='red' ); lines( singleVoxel_selectedOutput[9:length(singleVoxel_selectedOutput)], col='blue' );
    
    return( singleVoxel_selectedOutput )
  }
  
  #system.time( aaa <- sapply( 1:10, voxelModel ) )
  print( sprintf( 'fitting in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
  library(parallel)
  detectCores()
  nCores <- 4
  cl <- makeCluster(nCores, type='FORK')
  runIndexVoxel <- 1:length(selIdxVoxel)
  storeTimePar <- system.time( outModel <- parSapply(cl, runIndexVoxel, voxelModel ) )
  stopCluster(cl)
  print( storeTimePar )
  
  if (modelFitCounter==1) { outModelLoop <- outModel }
  if (modelFitCounter>1) { 
    selectedCols <- outModel[7,] > outModelLoop[7,] #where r2 is store (after 6 grid parameters, x,y,sigmaPos,sigmaNeg, f test comparison and p value comparison)
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
}

# add a eye fovea parameter???

#outModelLoop <- outMatrix
storeAllPred <- t(outModelLoop[ seq(1,10), ])

## extract surround size ##

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

#
storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 11:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(14+dim(paramsArrayPredicted)[2], dim(tsTransposedAll)[2] ) )


print('save linear step...')
polCoords <- cart2pol( storeAllPred[,c(2,1)] ) 
storeAllPredOut[1:14,selIdxVoxel] <- t( cbind( storeAllPred, polCoords, FWHM ) )
storeAllPredOut[15:dim(storeAllPredOut)[1],selIdxVoxel] <- t( paramsArrayPredicted[ selIdxVoxel,  ] )
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

labels <- c('x-Pos-eye','y-Pos-eye','sigmaPos-eye','sigmaNeg-eye','f-test','pval-f-test','varExp-eye','intercept-eye','slope_eye','slope_vis','theta-eye','radius-eye','fwhmCenter-eye','surroundSize-eye',
            'x-Pos-vis','y-Pos-vis','sigmaPos-vis','sigmaNeg-vis','hrf_a1-vis','hrf_a2-vis','varExp-vis','intercept-vis','slope-only-vis','theta-vis','radius-vis','fwhmCenter-vis','surroundSize-vis')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

