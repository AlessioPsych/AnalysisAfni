args <- commandArgs(T)
print( args )

## prepare the stimuli with the portrait as well

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GVW19')

#args <- c('greyMaskPrf10.nii.gz', 'meanTs_eyeMovement_topUp_detrend_res.nii', 'eyeNoBorder', '1', '0', '0.166', '2', '1')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
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

setwd('/analyse/Project0226/dataSummary')
if (stimType==1) { 
  arrayStim <- scan( 'eyeMovingStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==2) { 
  arrayStim <- scan( 'eyeFixStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) 
}
if (stimType==3) { 
  arrayStim <- scan( 'eyeFixStim_border_kat.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) )
}
if (stimType==4) { 
  arrayStim <- scan( 'prfStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1860,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==5) { 
  arrayStim <- scan( 'eyeFixStim_border_kat_disapp.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==6) { 
  arrayStim <- scan( 'eyeFixStim_border_occluded.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==7) { 
  arrayStim <- scan( 'prfStim_borders_fixation_kat.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1860,135) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==8) { 
  arrayStim <- scan( 'prfStim_borders_fixation_kat_occluded.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1860,135) ), c( 3, 1, 2 ) ) # eye movement
}



stimMatFlip <- aperm( stimMat[ dim(stimMat)[1]:1,, ], c(2,1,3) )

x11( width=3, height=3 )
for ( snap in 1:dim(stimMat)[3] ) {
  image( stimMatFlip[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01)
}
stimSeq <- stimMatFlip
x <- seq(-10,10,length.out = dim(stimSeq)[1] )
y <- seq(-5.5,5.5,length.out = dim(stimSeq)[2] )
rm( stimMat )

#this part of the code builds a matrix with all the possible prediction tested, for both models at this stage
if (fineFit==0) {
  xElements <- 6
  yElements <- 6
  sigmaArrayPositiveElements <- 5
  multParElements <- 4
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
}
if (fineFit==1) {
  xElements <- 12
  yElements <- 12
  sigmaArrayPositiveElements <- 9
  multParElements <- 5
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
}
addSpace <- abs( min(x) )*0.1
print('build prediction...')
xPosFit <- seq( -12.5, 12.5, length.out=xElements )
yPosFit <- seq( -8, 8, length.out=yElements )
sigmaArrayPositive <- seq( 0.25, 9, length.out=sigmaArrayPositiveElements )
#xPosFit <- seq( -2, 2, length.out=xElements )
#yPosFit <- seq( -2, 2, length.out=yElements )
#sigmaArrayPositive <- seq( 0.25, 2, length.out=sigmaArrayPositiveElements )
if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
par_hrf_a1 <- seq( 6, 9, length.out=hrfDelayOnsetElements )
par_hrf_a2 <- seq( 12, 15, length.out=hrfDelayUnderShootElements )
if (flagSurround==1) { multPar <- seq(0,0.8, length.out = multParElements) }
if (flagSurround==0) { multPar <- 0 }
predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar )
keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]

# this part builds a matrix with the starting and ending prediction index to fit for each loop in the next for loop
if (dim( predictionGridGlobal )[1] <= 250 ) {
  totalIterations <- 3
}
if (dim( predictionGridGlobal )[1] > 250 ) {
  totalIterations <- ceil( dim( predictionGridGlobal )[1]/200 )
}
limitsPrediction <- round( seq(1,dim(predictionGridGlobal)[1], length.out=totalIterations ) ) 
limitsPrediction[1] <- 1
limitsPrediction[ length(limitsPrediction) ] <- dim(predictionGridGlobal)[1]
limitsPredictionMatrix <- array( 0, c( (length(limitsPrediction)-1) , 2 ) )
limitsPredictionMatrix[,1] <- limitsPrediction[1:(length(limitsPrediction)-1)]
limitsPredictionMatrix[,2] <- limitsPrediction[2:(length(limitsPrediction))]
limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] <- limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
runIndexPredictions <- seq( 1:dim(limitsPredictionMatrix)[1] )

indexVol <- meanEpi$brk[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
selIdxVoxel <- which( indexArray == 1 )
tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series

#multVector <- c(1)#c(0.8, 0.9, 1, 1.1, 1.2)
#modelFitCounter <- 1
for (modelFitCounter in 1:length(runIndexPredictions)) { #
  
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
    #to test: 
    #inputPredictionGrid <- predictionGrid; indexPrediction <- 2
    prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
    
    #prfPar <- c(2,2,1.5,2,6,12,0.3)
    hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
    
    a <- dnorm( x, prfPar[1], prfPar[3] )
    b <- dnorm( y, prfPar[2], prfPar[3] )
    a[ x<qnorm( 0.01, mean=prfPar[1], sd=prfPar[3] ) | x>qnorm( 0.99, mean=prfPar[1], sd=prfPar[3] ) ] <- 0 #cut the tails x
    b[ y<qnorm( 0.01, mean=prfPar[2], sd=prfPar[3] ) | y>qnorm( 0.99, mean=prfPar[2], sd=prfPar[3] ) ] <- 0 #cut the tails y 
    imgCenter <- scaleData( tcrossprod(a,b), 1, 0 )
    #cut the tails, to plot:
    #par(mfrow=c(2,1))
    #plot( a ~ x, bty='n', las=1 ); abline( v=qnorm( 0.01, mean=prfPar[1], sd=prfPar[3] ), col='red', lty=2, lwd=1  ); abline( v=qnorm( 0.99, mean=prfPar[1], sd=prfPar[3] ), col='red', lty=2, lwd=1  )
    #plot( b ~ y, bty='n', las=1 ); abline( v=qnorm( 0.01, mean=prfPar[2], sd=prfPar[3] ), col='red', lty=2, lwd=1  ); abline( v=qnorm( 0.99, mean=prfPar[2], sd=prfPar[3] ), col='red', lty=2, lwd=1  )
    
    a <- dnorm( x, prfPar[1], prfPar[4] )
    b <- dnorm( y, prfPar[2], prfPar[4] )
    a[ x<qnorm( 0.01, mean=prfPar[1], sd=prfPar[4] ) | x>qnorm( 0.99, mean=prfPar[1], sd=prfPar[4] ) ] <- 0 #cut the tails x
    b[ y<qnorm( 0.01, mean=prfPar[2], sd=prfPar[4] ) | y>qnorm( 0.99, mean=prfPar[2], sd=prfPar[4] ) ] <- 0 #cut the tails y 
    imgSurround <- scaleData( tcrossprod(a,b), 1, 0)*prfPar[7]
    
    r <- imgCenter - imgSurround
    #image( r, col=grey.colors(1000) )
    rMat <- array(r)
    
    if ( sum( is.na(rMat) ) == 0 ) {
      #trMat <- t(stimSeqMat)
      #system.time( predictionLoop <- as.numeric( trMat%*%rMat ) )
      predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
      pConv <- conv( predictionLoop, hrf )
      pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
      tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('linear') )
      returnPrediction <- round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 )
      #par(mfrow=c(1,3))
      #plot( r[,70], type='l' )
      #image( r, col=gray.colors(500) )
      #plot( returnPrediction, type='l' )
    }
    if ( sum( is.na(rMat) ) > 0 ) {
      returnPrediction <- rep(0,length(mriTime))
    }
    
    return( returnPrediction ) #### scale predictions betweeo 0 and 1 and round them to 5 digits
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
  controlPredictions <- apply( tsPrediction, 1, sum ) != 0 & apply( is.na(tsPrediction), 1, sum ) == 0
  tsPrediction <- tsPrediction[controlPredictions,]
  predictionGrid <- predictionGrid[controlPredictions,]
  print( dim( predictionGrid ) )
  
  if (dim( predictionGrid )[1]==0) { outMatrix <- array(0, c(dim(tsTransposedSel)[1]+10,dim(tsTransposedSel)[2])) }
  if (dim( predictionGrid )[1]!=0) { 
    
    #### linear fitting in parallel ####
    print( sprintf( 'linear fitting in iteration %1.0f of %1.0f ...', modelFitCounter, length(runIndexPredictions)  ) )
    # indexVol <- meanEpi$brk[,,,1]
    # indexArray <- array( indexVol, prod( dim( indexVol ) ) )
    # tsTransposedAll <- t( tsArray )
    # selIdxVoxel <- which( indexArray == 1 )
    # tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series
    
    totalVoxelsIterations <- ceil( length(selIdxVoxel) / 200 )
    limitsSelVoxels <- round( seq(1,length(selIdxVoxel), length.out=totalVoxelsIterations) ) #split between chuncks for parallel fitting (about 200 voxels each)
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
    #system.time( aaa <- lapply( 4:5, voxelModel ) )
    library(parallel)
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
  }  
  
  outModel <- outMatrix
  if (modelFitCounter==1) { outModelLoop <- outModel }
  if (modelFitCounter>1) { 
    selectedCols <- outModel[8,] > outModelLoop[8,] #where r2 is store (after 4 prf parameters and 2 hrf parameters)
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
}

#outModelLoop <- outMatrix #fix it from here
storeAllPred <- t(outModelLoop[ seq(1,10), ])
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
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 11:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(14, dim(tsTransposedAll)[2] ) )

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

labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multpar','varExp','intercept','slope','theta','radius','fwhmCenter','surroundSize')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

