

#args <- commandArgs(T)
#print( args )

#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')
#args <- c('bars_mask_res.nii.gz','barsTs_res.nii.gz','outTest_delme.nii.gz','/packages/afni/17.0.13','1')

## prepare the stimuli with the portrait as well

rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

mainDir <- getwd()
args <- c('meanTs_eyeMovement_topUp_res_mask.nii', 'meanTs_eyeMovement_topUp_res.nii', 'output_test_eye_test_fix_moving','/home/alessiof/abin', '1', '4','1','0.166')
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
arrayStim <- scan( 'eyeMovingStim.txt' )

setwd(mainDir)
#stimMat <- aperm( array( arrayStim, c(128,620,96) ), c( 1, 3, 2 ) )
stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
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
#stimMat <- stimMatInterp2
stimMatFlip <- stimMat[ ,dim(stimMat)[2]:1, ]
image( stimMatFlip[,,556], axes=FALSE )
#image( stimMat[,,56], axes=FALSE )
x11( width=4, height=4 )
for ( snap in 1:dim(stimMat)[3] ) {
  image( stimMat[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01) 
}
stimSeq <- stimMatFlip
x <- seq(-7,7,length.out = dim(stimSeq)[1] )
y <- seq(-10,10,length.out = dim(stimSeq)[2] )

print('get hrf...')
hrf <- canonicalHRF( seq(0,30,samplingTime) )

addSpace <- abs( min(x) )*0.5
print('build prediction...')
xPosFit <- seq( min(x)-addSpace, max(x)+addSpace, length.out=6 )
yPosFit <- seq( min(y)-addSpace, max(y)+addSpace, length.out=6 )
sigmaArrayPositive <- seq( 0.15, 6, length.out=6 )
if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative )
keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
predictionGrid <- predictionGridTemp[ keepPredictionIdx, ]

tsPrediction <- array( 0, c( dim(predictionGrid)[1], dim( ts$brk )[4] ) )
stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
incrementalCounter <- dim(predictionGrid)[1] * 0.05
timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(ts$brk)[4] )

for (nPrediction in 1:dim(predictionGrid)[1] ) {
  
  prfPar <- as.numeric( predictionGrid[nPrediction,] )
  
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
  rMat <- round( array(r), 5 )

  #trMat <- t(stimSeqMat)
  #system.time( predictionLoop <- as.numeric( trMat%*%rMat ) )
  predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) )
  pConv <- conv( predictionLoop, hrf )
  pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
  tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('linear') )
  tsPrediction[nPrediction,] <- tsPredictionMriInterp / max( tsPredictionMriInterp )
  if ( nPrediction > incrementalCounter ) {
    incrementStep <- dim(predictionGrid)[1] * 0.05
    incrementalCounter <- incrementalCounter + incrementStep
    percentageDone <- incrementalCounter / dim(predictionGrid)[1] * 100
    cat(paste('*',""))
  }
}

controlPredictions <- apply( tsPrediction, 1, sum ) != 0 

tsPrediction <- tsPrediction[controlPredictions,]
predictionGrid <- predictionGrid[controlPredictions,]


generatePrediction <- function( indexPrediction ) {
  prfPar <- as.numeric( predictionGrid[indexPrediction,] )
  
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
  predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) )
  pConv <- conv( predictionLoop, hrf )
  pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
  tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('linear') )
  #tsPrediction[nPrediction,] <- tsPredictionMriInterp / max( tsPredictionMriInterp )
  return( round( tsPredictionMriInterp / max( tsPredictionMriInterp ), 5 ) )
}
library(parallel)
detectCores()
nCores <- 8
cl <- makeCluster(nCores, type='FORK')
storeTimePar <- system.time( outParallel <- parSapply(cl, 1:dim(predictionGrid)[1], generatePrediction ) )
stopCluster(cl)
storeTimeSerial <- system.time( sapply(1:dim(predictionGrid)[1], generatePrediction ) )
print( storeTimePar )
print( storeTimeSerial )



### downsample predictions to epi TR
#timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
#mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(ts$brk)[4] )
#tsPredictionMriInterp <- array(0, c( dim(predictionGrid)[1], length(mriTime) ) )
#for ( nPrediction in 1:dim(predictionGrid)[1] ) {
#  tsPredictionMriInterp[nPrediction,] <- interp1( x=timeStimuli, y=tsPrediction[nPrediction,], xi=mriTime, method=c('linear') )
#  tsPredictionMriInterp[nPrediction,] <- scaleData( tsPredictionMriInterp[nPrediction,], 1, 0)
#}
#tsPrediction <- tsPredictionMriInterp

#scale predictions to 1

#tsPrediction <- tsPrediction / apply( tsPrediction, 1, max )

print('linear fitting...')
indexVol <- meanEpi$brk[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
tsTransposedOriginal <- t( tsOriginalArray )
tsOriginalAverage <- apply( tsTransposedOriginal, 2, mean )
nVoxels <- dim( tsTransposedAll )[2]
#fitPart <- ceiling( seq( 1, nVoxels, length.out = 401 ) )
fitPart <- ceiling( seq( 1, nVoxels, by = 50 ) )
fitPart <- fitPart[ 1:(length(fitPart)-2) ]
storeAllPred <- array(0,c(nVoxels,7))
storeAllExpectedTs <- array(0,dim(tsTransposedAll))

for (dataPart in 1:(length(fitPart)-1) ) {
  
  cat( paste( c('Iteration:',dataPart, 'out of:',length(fitPart)-1, "" ) ) )
  print('...')
  tsTransposedFull <- tsTransposedAll[ , seq( fitPart[dataPart], fitPart[dataPart+1] ) ]
  arrayFull <- indexArray[ seq( fitPart[dataPart], fitPart[dataPart+1] ) ]
  #tsSd <- apply( tsTransposedFull, 2, sd  )
  
  idxFit <- arrayFull == 1
  if (sum(idxFit)>1) {
    tsTransposed <- tsTransposedFull[,idxFit]
    
    #tsTransposedMean <- apply( tsTransposed, 2, mean )
    #ssTot <- apply( (tsTransposed-tsTransposedMean)^2, 2, sum)
    
    if (sum(idxFit)>1) { 
      tsTransposedMean <- apply( tsTransposed, 2, mean ) 
      ssTot <- apply( (tsTransposed-tsTransposedMean)^2, 2, sum)
    }
    if (sum(idxFit)==1) { 
      tsTransposedMean <- mean( tsTransposed ) 
      ssTot <- sum((tsTransposed-tsTransposedMean)^2)
    }
    #tsMeans <- apply( tsTransposed, 2, mean  )
    #tsMeanScaled <- scale( tsTransposed, center=tsMeans, scale=FALSE  ) 
    #ssTot <- apply( tsMeanScaled^2, 2, sum)
    progress <- 0.05
    
    for (k in 1:dim(tsPrediction)[1] ) {
      
      dMat <- cbind( tsPrediction[k,] )
      
      if (fitIntercept==1) { dMat01 <- cbind( rep(1,length(dMat)), dMat ) }
      if (fitIntercept==0) { dMat01 <- dMat }
      #a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,tsTransposed) )
      a <- solve( qr(dMat01), tsTransposed)
      
      ### remove the intercept fit, it is better ############
      
      #system.time( a2 <- solve(t(dMat01) %*% dMat01, t(dMat01) %*% tsTransposed) )
      #system.time( a3 <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,tsTransposed) ))
      #t( ginv(dMat01) )
      #system.time( a <- solve( qr(dMat01), tsTransposed) )
      #system.time( a <- lm( tsTransposed ~ dMat01 ) )
      
      #expectedTs01 <- dMat01%*%a
      expectedTs <- crossprod( t(dMat01), a )
      residualsSquared <- ( tsTransposed - expectedTs )^2
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot)
      
      if (fitIntercept==0) { a <- rbind( rep(1000, dim(a)[2] ), a ) }
      
      if (k==1) {
        if (sum(idxFit)>1) { 
          storePred <- matrix( rep( as.numeric( predictionGrid[k,] ), dim(tsTransposed)[2] ), nrow = dim(tsTransposed)[2], ncol=dim(predictionGrid)[2], byrow=TRUE  )
        }
        if (sum(idxFit)==1) { 
          storePred <- matrix( rep( as.numeric( predictionGrid[k,] ), 1 ), nrow = 1, ncol=dim(predictionGrid)[2], byrow=TRUE  )
        }
        #storePred <- matrix( rep( as.numeric( predictionGrid[k,] ), dim(tsTransposed)[2] ), nrow = dim(tsTransposed)[2], ncol=dim(predictionGrid)[2], byrow=TRUE  )
        storePred <- cbind( storePred, r2, t( a ) )
        storePred[ storePred[ , dim(storePred)[2] ] < 0, c( (dim(storePred)[2]-1):dim(storePred)[2] ) ] <- 0 
        storeFit <- expectedTs
      }
      if (k>1) {
        updateIdx <- r2>storePred[,5] & a[2,]>0
        if ( sum( updateIdx ) > 0 ) {
          storePredUpdated <- matrix( rep( as.numeric( predictionGrid[k,] ), sum( updateIdx ) ), nrow = sum( updateIdx ), ncol=dim(predictionGrid)[2], byrow=TRUE  )
          storePredUpdated <- cbind( storePredUpdated, r2[updateIdx], t( a[,updateIdx] ) )
          storePred[ updateIdx, ] <- storePredUpdated
          storeFit[, updateIdx ] <- expectedTs[,updateIdx]
        }
      }
      if (k>dim(tsPrediction)[1]*progress) {
        cat(paste('*',""))
        progress <- progress + 0.05
      }
    }
    
    storePredFilled <- array( 0, c( dim(tsTransposedFull)[2], dim(storePred)[2] )  )
    storeFitFilled <- array( 0, dim(tsTransposedFull) )
    storePredFilled[idxFit,] <- storePred
    storeFitFilled[,idxFit] <- storeFit
    storeAllPred[ fitPart[dataPart]:fitPart[dataPart+1], ] <- storePredFilled
    storeAllExpectedTs[, fitPart[dataPart]:fitPart[dataPart+1] ] <- storeFitFilled
  }
  print('...')
}



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





#take it from here, function to  optimize the fit, maybe a couple of times

print('fine fit...')
fineFitPredictor <- function( paramsArray ) {
    if (flagSurround==1) { 
      xCoord <- paramsArray[1]
      yCoord <- paramsArray[2]
      sizeCenter <- paramsArray[3]
      sizeSurround <- paramsArray[4]
      #intercept <- paramsArray[5]
      #slope <- paramsArray[6]
      a <- dnorm( x, xCoord, sizeCenter )
      b <- dnorm( y, yCoord, sizeCenter )
      imgCenter <- tcrossprod(a,b)
      a <- dnorm( x, xCoord, sizeSurround )
      b <- dnorm( y, yCoord, sizeSurround )
      imgSurround <- tcrossprod(a,b)
      r <- imgCenter - imgSurround
    }
    if (flagSurround==0) { 
      xCoord <- paramsArray[1]
      yCoord <- paramsArray[2]
      sizeCenter <- paramsArray[3]
      #intercept <- paramsArray[4]
      #slope <- paramsArray[5]
      a <- dnorm( x, xCoord, sizeCenter )
      b <- dnorm( y, yCoord, sizeCenter )
      imgCenter <- tcrossprod(a,b)
      r <- imgCenter
    }
    rMat <- array(r)
    predictionLoop <- as.numeric( crossprod( stimSeqMat, rMat ) )
    pConv <- conv( predictionLoop, hrf )
    pConv <- pConv[ 1 : dim(stimSeq)[3] ]
    ### downsampling step ###
    timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out=dim(stimSeq)[3] )
    mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out=dim(ts$brk)[4] )
    pConv <- interp1( x=timeStimuli, y=pConv, xi=mriTime, method=c('linear') )
    ### make max to 1 ###
    pConv <- pConv / max(pConv)
    #pConvTrim <- intercept + pConv * slope
    return(pConv)  
}


test



xMod <- c(-0.2,0,0.2)
yMod <- xMod
sigmaMod <- xMod
surroundMod <- xMod 
multiplicativeGrid <- expand.grid(xMod, yMod, sigmaMod, surroundMod )
rep(  ) 





# from here, check the position of the slope in the storePred matrix

print('save linear step...')
polCoords <- cart2pol( storeAllPred[,c(1,2)] )
storeAllPredOut <- cbind( storeAllPred, polCoords, FWHM )
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

rParameters <- array( storeAllPredOut, c( dim(ts$brk)[1:3], dim(storeAllPredOut)[2] ) )
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

labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','varExp','intercept','slope','theta','radius','fwhmCenter','surroundSize')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}



print('non-linear fit...')
nLinNls <- function( paramsArray ) {
  if (fitIntercept==1) {
    if (flagSurround==1) { 
      xCoord <- paramsArray[1]
      yCoord <- paramsArray[2]
      sizeCenter <- paramsArray[3]
      sizeSurround <- paramsArray[4]
      intercept <- paramsArray[5]
      slope <- paramsArray[6]
      a <- dnorm( x, xCoord, sizeCenter )
      b <- dnorm( y, yCoord, sizeCenter )
      imgCenter <- tcrossprod(a,b)
      a <- dnorm( x, xCoord, sizeSurround )
      b <- dnorm( y, yCoord, sizeSurround )
      imgSurround <- tcrossprod(a,b)
      r <- imgCenter - imgSurround
    }
    if (flagSurround==0) { 
      xCoord <- paramsArray[1]
      yCoord <- paramsArray[2]
      sizeCenter <- paramsArray[3]
      intercept <- paramsArray[4]
      slope <- paramsArray[5]
      a <- dnorm( x, xCoord, sizeCenter )
      b <- dnorm( y, yCoord, sizeCenter )
      #imgCenter <- a%*%t(b)
      imgCenter <- tcrossprod(a,b)
      r <- imgCenter
    }
    rMat <- array(r)
    predictionLoop <- as.numeric( crossprod( stimSeqMat, rMat ) )
    pConv <- conv( predictionLoop, hrf )
    pConv <- pConv[ 1 : dim(stimSeq)[3] ]
    ### downsampling step ###
    timeStimuli <- seq( 0, dim(stimSeq)[3]*0.5-0.5, by=0.5 )
    mriTime <- seq( 0, dim(stimSeq)[3]*0.5-0.5, by=2 )
    pConv <- interp1( x=timeStimuli, y=pConv, xi=mriTime, method=c('linear') )
    ### downsampling step ###
    pConv <- pConv / max(pConv)
    #pConvTrim <- intercept + pConv[ 1 : dim(stimSeqMat)[2] ] * slope
    pConvTrim <- intercept + pConv * slope
    return(pConvTrim)  
  }
  if (fitIntercept==0) {
    if (flagSurround==1) { 
      xCoord <- paramsArray[1]
      yCoord <- paramsArray[2]
      sizeCenter <- paramsArray[3]
      sizeSurround <- paramsArray[4]
      #intercept <- paramsArray[5]
      slope <- paramsArray[5]
      a <- dnorm( x, xCoord, sizeCenter )
      b <- dnorm( y, yCoord, sizeCenter )
      imgCenter <- tcrossprod(a,b)
      a <- dnorm( x, xCoord, sizeSurround )
      b <- dnorm( y, yCoord, sizeSurround )
      imgSurround <- tcrossprod(a,b)
      r <- imgCenter - imgSurround
    }
    if (flagSurround==0) { 
      xCoord <- paramsArray[1]
      yCoord <- paramsArray[2]
      sizeCenter <- paramsArray[3]
      #intercept <- paramsArray[4]
      slope <- paramsArray[4]
      a <- dnorm( x, xCoord, sizeCenter )
      b <- dnorm( y, yCoord, sizeCenter )
      #imgCenter <- a%*%t(b)
      imgCenter <- tcrossprod(a,b)
      r <- imgCenter
    }
    rMat <- array(r)
    predictionLoop <- as.numeric( crossprod( stimSeqMat, rMat ) )
    pConv <- conv( predictionLoop, hrf )
    pConv <- pConv[ 1 : dim(stimSeq)[3] ]
    ### downsampling step ###
    timeStimuli <- seq( 0, dim(stimSeq)[3]*0.5-0.5, by=0.5 )
    mriTime <- seq( 0, dim(stimSeq)[3]*0.5-0.5, by=2 )
    pConv <- interp1( x=timeStimuli, y=pConv, xi=mriTime, method=c('linear') )
    ### downsampling step ###
    
    pConv <- pConv / max(pConv)
    pConvTrim <- intercept + pConv * slope
    #pConvTrim <- pConv[ 1 : dim(stimSeqMat)[2] ] * slope
    return(pConvTrim)  
  }
}

# to debug and visualize # from here
sortR2 <- sort( storeAllPred[,5], index.return=TRUE, decreasing = TRUE )
X11( width=12, height=8 )
par(mfrow=c(3,3))
for (nSubplot in 1:9) {
  selIndexEpi <- sortR2$ix[nSubplot]
  observedData <- tsTransposedAll[,selIndexEpi]
  paramsArray <- storeAllPred[selIndexEpi,]
  if (fitIntercept==1) {
    if (flagSurround==1) { paramsArray <- paramsArray[c(1,2,3,4,6,7)] } # x, y, sigmaPos, sigmaNeg, intercept, slope
    if (flagSurround==0) { paramsArray <- paramsArray[c(1,2,3,6,7)] } # x, y, sigmaPos, intercept, slope
  }
  if (fitIntercept==0) {
    if (flagSurround==1) { paramsArray <- paramsArray[c(1,2,3,4,7)] } # x, y, sigmaPos, sigmaNeg, slope
    if (flagSurround==0) { paramsArray <- paramsArray[c(1,2,3,7)] } # x, y, sigmaPos, slope
  }
  hrfFun <- hrf
  plot( observedData, type='p', pch=3, bty='n', las=1, xlab='TR number', ylab='detrended BOLD', col='gray50' )
  lines( nLinNls(paramsArray), col='blue', lwd=2 )
}
dev.copy2pdf( file=sprintf( '%s_linear_estimates_examples.pdf', outSuffix ), width=12, height=8 ) 
dev.off()


#paramsNlinTest <- t(outParallel)[96,]
#paramsNlinTest <- paramsNlinTest[c(1,2,3,4,6,7)]
#idxNonLin <- idx_to_run_sapply[96]
#paramsArray <- storeAllPred[idx_to_run_sapply[96],]
#paramsArray <- paramsArray[c(1,2,3,4,6,7)]
#observedData <- tsTransposedAll[,idx_to_run_sapply[96]]
#plot( observedData, type='p', pch=8 )
#lines( nLinNls(paramsNlinTest), col='red' )
#lines( nLinNls(paramsArray), col='blue' )



nLinFit <- function( startParams, observedData ) {
  pConvTrim <- nLinNls( startParams )
  out <- sum( ( pConvTrim - observedData )^2 ) 
  return(out)
}

#nLinFit( paramsArray, observedData )

fToComp <- function( startParams, dataToFit ) {
  options(warn=-1)
  #optOut01 <- optim( par=startParams, fn=nLinFit,
  #                   xInput = xInput, yInput = yInput, hrfInput = hrfInput, stimSeqInput=stimSeqInput, constParInput=constPar,
  #                   method='CG', control=list(ndeps=c( 0.02, 0.02, 0.02, 0.05, 200, 200 ) , reltol=1e-1, maxit=10 )
  
  #if (startParams[5] < 0) {
  #    lowBound <- c( startParams[1]*0.5, startParams[2]*0.5, startParams[3]*0.98, startParams[4]*0.5, startParams[5]*1.5, startParams[6]*0.5 )
  #    upBound <- c( startParams[1]*1.5, startParams[3]*0.95, startParams[3]*1.5, startParams[4]*1.5, startParams[5]*0.5, startParams[6]*1.5 )
  #}
  #if (startParams[5] > 0) {
  #    lowBound <- c( startParams[1]*0.5, startParams[2]*0.5, startParams[3]*0.98, startParams[4]*0.5, startParams[5]*0.5, startParams[6]*0.5 )
  #    upBound <- c( startParams[1]*1.5, startParams[3]*0.95, startParams[3]*1.5, startParams[4]*1.5, startParams[5]*1.5, startParams[6]*1.5 )
  #}
  
  # loP <- constPar
  # 
  # if (startParams[5] < 0) {
  #     lowBound <- c( startParams[1]*0.5, 0.2,                0.5, startParams[4]*0.5, startParams[5]*1.5, startParams[6]*0.5 )
  #     upBound <- c( startParams[1]*1.5, 7, 8, startParams[4]*1.5, startParams[5]*0.5, startParams[6]*1.5 )
  # }
  # 
  # if (startParams[5] > 0) {
  #     lowBound <- c( startParams[1]*0.5, 0.2, startParams[3]*0.98, startParams[4]*0.5, startParams[5]*0.5, startParams[6]*0.5 )
  #     upBound <- c( startParams[1]*1.5, startParams[3]*0.95, startParams[3]*1.5, startParams[4]*1.5, startParams[5]*1.5, startParams[6]*1.5 )
  # }
  # 
  # startingParams <- startParams
  # startingParams[1] <- startingParams[1]+startingParams[1]*loP
  # startingParams[2] <- startingParams[2]+startingParams[2]*loP
  # startingParams[3] <- startingParams[3]+startingParams[3]*loP
  
  #print( startParams )
  #print( startingParams )
  
  #lowBound <- startParams - abs( startParams )*0.15
  #upBound <- startParams + abs( startParams )*0.15
  #upBound[3] <- lowBound[4]*0.95
  
  
  #optOut01 <- optim( par=startParams, fn=nLinFit, observedData=dataToFit, lower = lowBound , upper = upBound,
  #                   method='L-BFGS-B', control=list(ndeps=c( 0.03, 0.03, 0.03, 0.06, 100, 100 ), reltol=1e-3, maxit=5 )
  if (fitIntercept==1) { 
    if (flagSurround==1) { # x, y, sigmaPos, sigmaNeg, intercept, slope
      optOut01 <- optim( par=startParams, fn=nLinFit, observedData=dataToFit,
                         method='CG', control=list(ndeps=c( 0.03, 0.03, 0.03, 0.06, 100, 100 ), reltol=1e-3, maxit=12 )
      )
    }
    if (flagSurround==0) { # x, y, sigmaPos, intercept, slope
      optOut01 <- optim( par=startParams, fn=nLinFit, observedData=dataToFit,
                         method='CG', control=list(ndeps=c( 0.03, 0.03, 0.03, 100, 100 ), reltol=1e-3, maxit=12 )
      )
    }
  }
  if (fitIntercept==0) {
    if (flagSurround==1) { # x, y, sigmaPos, sigmaNeg, slope
      optOut01 <- optim( par=startParams, fn=nLinFit, observedData=dataToFit,
                         method='CG', control=list(ndeps=c( 0.03, 0.03, 0.03, 0.06, 100 ), reltol=1e-3, maxit=12 )
      )
    }
    if (flagSurround==0) { # x, y, sigmaPos, slope
      optOut01 <- optim( par=startParams, fn=nLinFit, observedData=dataToFit,
                         method='CG', control=list(ndeps=c( 0.03, 0.03, 0.03, 100 ), reltol=1e-3, maxit=12 )
      )
    }
  }
  
  return( optOut01$par )
  options(warn=0)
}

#fToComp( paramsArray, observedData )

runFit <- function( idx_to_run, idxSelected ) {
  selIdx <- idx_to_run[idxSelected]
  paramsArray <- storeAllPred[ selIdx, ]
  if (fitIntercept==1) {
    if (flagSurround==1) { 
      paramsArray <- paramsArray[c(1,2,3,4,6,7)] # x, y, sigmaPos, sigmaNeg, intercept, slope
      observedData <- tsTransposedAll[,selIdx]
      uuu <- tryCatch( fToComp( paramsArray, observedData ),
                       error = function(e) 'error' )
    }
    if (flagSurround==0) { 
      paramsArray <- paramsArray[c(1,2,3,6,7)] # x, y, sigmaPos, intercept, slope
      observedData <- tsTransposedAll[,selIdx]
      uuu <- tryCatch( fToComp( paramsArray, observedData ),
                       error = function(e) 'error' )
      uuuOut <- rep(0,6)
      uuuOut[1:3] <- uuu[1:3]
      uuuOut[4] <- 1000
      uuuOut[5:6] <- uuu[4:5]
      uuu <- uuuOut
    }
  }
  if (fitIntercept==0) {
    if (flagSurround==1) { 
      paramsArray <- paramsArray[c(1,2,3,4,7)] # x, y, sigmaPos, sigmaNeg, slope
      observedData <- tsTransposedAll[,selIdx]
      uuu <- tryCatch( fToComp( paramsArray, observedData ),
                       error = function(e) 'error' )
      uuuOut <- rep(0,6)
      uuuOut[1:4] <- uuu[1:4]
      uuuOut[5] <- 1000
      uuuOut[6] <- uuu[5]
      uuu <- uuuOut
    }
    if (flagSurround==0) { 
      paramsArray <- paramsArray[c(1,2,3,7)] # x, y, sigmaPos, slope
      observedData <- tsTransposedAll[,selIdx]
      uuu <- tryCatch( fToComp( paramsArray, observedData ),
                       error = function(e) 'error' )
      uuuOut <- rep(0,6)
      uuuOut[1:3] <- uuu[1:3]
      uuuOut[4:5] <- c(1000,1000)
      uuuOut[6] <- uuu[4]
      uuu <- uuuOut
    }
  }
  if (uuu=='error') {
    errorCatch <- rep(0,6)
    dataOrig <- nLinNls( paramsArray )
    r2Old <- summary( lm( observedData ~ dataOrig ) )$r.squared
    return( c( storeAllPred[ selIdx, ], storeAllPred[ selIdx, 5], r2Old ) )
  }
  if (uuu!='error') {
    dataOrig <- nLinNls( paramsArray )
    r2Old <- summary( lm( observedData ~ dataOrig ) )$r.squared
    dataNew <- nLinNls( uuu )
    r2New <- summary( lm( observedData ~ dataNew ) )$r.squared
    return( c( uuu[1:4], storeAllPred[ selIdx, 5], uuu[5:6], r2Old, r2New ) )
  }
}

#keepPredictors <- apply( storeAllPredOut == 0, 1, sum )
#idx_to_run_sapply <- which( indexArray == 1 & storeAllPredOut[,5] > 0.10 & keepPredictors==0 )

idx_to_run_sapply <- which( storeAllPredOut[,5] > 0.10 )

#idx_to_run <- idx_to_run_sapply
#idxSelected <- 2
#runFit(idx_to_run_sapply,2)
#selIdx <- idx_to_run_sapply[idxSelected]
#summary( lm(storeAllExpectedTs[,selIdx]~observedData) )$r.squared

no_cores <- 4
cl <- makeCluster(no_cores, type='FORK')
storeTime <- system.time( outParallel <- parSapply( cl, c(1:length(idx_to_run_sapply)), runFit, idx_to_run=idx_to_run_sapply ) )
#storeTime <- system.time( outParallel <- parSapply( cl, c(1:1000), runFit, idx_to_run=idx_to_run_sapply ) )
stopCluster(cl)
print( storeTime )

if (flagSurround==1) { 
  progress <- 0.05
  FWHMNlin <- array(0, c( dim(storeAllPred)[1], 2 ) )
  for (k in 1:dim(outParallel)[2]) {
    
    parameters <- outParallel[,k]
    idxOverall <- idx_to_run_sapply[k]
    
    if ( sum(parameters)!=0 ) {
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
      FWHMNlin[idxOverall,] <- c(FWHMcenter,FWHMsurround)
      
      paramsNlinTest <- t(outParallel)[k,]
      if (fitIntercept==1) {
        paramsNlinTest <- paramsNlinTest[c(1,2,3,4,6,7)] # x, y, sigmaPos, sigmaNeg, intercept, slope
      }
      if (fitIntercept==0) {
        paramsNlinTest <- paramsNlinTest[c(1,2,3,4,7)] # x, y, sigmaPos, sigmaNeg, slope
      }
      idxNonLin <- idx_to_run_sapply[k]
      tsNlin <- nLinNls(paramsNlinTest)
      storeAllExpectedTs[,idxNonLin] <- tsNlin
    }
    
    if (k>dim(outParallel)[2]*progress) {
      cat(paste('*',""))
      progress <- progress + 0.05
    }
    
  }
}
if (flagSurround==0) {
  FWHMNlin <- array(0, c( dim(storeAllPred)[1], 2 ) )
  
  #parameters <- outParallel[,k]
  #idxOverall <- idx_to_run_sapply[k]
  
  FWHMNlin[idx_to_run_sapply,] <- outParallel[3,]
  FWHMNlin[2,] <- 1000
  
  for (k in 1:dim(outParallel)[2]) {
    parameters <- outParallel[,k]
    idxOverall <- idx_to_run_sapply[k]
    
    if ( sum(parameters)!=0 ) {
      
      paramsNlinTest <- t(outParallel)[k,]
      if (fitIntercept==1) {
        paramsNlinTest <- paramsNlinTest[c(1,2,3,6,7)] # x, y, sigmaPos, intercept, slope
      }
      if (fitIntercept==0) {
        paramsNlinTest <- paramsNlinTest[c(1,2,3,7)] # x, y, sigmaPos, slope
      }
      idxNonLin <- idx_to_run_sapply[k]
      tsNlin <- nLinNls(paramsNlinTest)
      storeAllExpectedTs[,idxNonLin] <- tsNlin
    }
  }
}

#test the two participants

emptyNlinOut <- array( 0, c( dim(storeAllPred)[1], dim(outParallel)[1]  ) )
emptyNlinOut[ idx_to_run_sapply, ] <-  t( outParallel )

print('save non-linear step...')
polCoords <- cart2pol( storeAllPred[,c(1,2)] )
polCoordsNlin <- cart2pol( emptyNlinOut[,c(1,2)] )
percentSignalChange <- emptyNlinOut[,7]/tsOriginalAverage*100 #tsOriginalAverage computed on line 115
storeAllPredOut <- cbind( storeAllPred, FWHM, emptyNlinOut, polCoords, polCoordsNlin, FWHMNlin, percentSignalChange )
fileTs <- sprintf('%s_PredixtedTs_nLin.nii.gz',outSuffix) 
fileParams <- sprintf('%s_params_nlin.nii.gz',outSuffix)

rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( '__tt_expectedTs.nii.gz', brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_expectedTs.nii.gz', orientValue[1,1], fileTs )
system( instr)
system( 'rm __tt_expectedTs.nii.gz' )

rParameters <- array( storeAllPredOut, c( dim(ts$brk)[1:3], dim(storeAllPredOut)[2] ) )
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

labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','varExp','intercept','slope','fwhmCenter','surroundSize','x-Pos_n_lin', 'y-Pos_n_lin',
            'sigmaPos_n_lin', 'sigmaNeg_n_lin','varExpRecomputed','intercept_n_lin', 'slope_n_lin', 'varExp_lin', 'varExp_n_lin', 'theta', 'ecc', 'theta_nlin', 'ecc_nlin', 'fwhmCenter_n_lin', 'surroundSize_n_lin', 'BOLDsignalChange' )
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}


sortR2 <- sort( storeAllPredOut[,17], index.return=TRUE, decreasing = TRUE )
X11( width=12, height=8 )
par(mfrow=c(3,3))
for (nSubplot in 1:9) {
  selIndexEpi <- sortR2$ix[nSubplot]
  observedData <- tsTransposedAll[,selIndexEpi]
  paramsArray <- storeAllPredOut[selIndexEpi,]
  if (fitIntercept==1) {
    if (flagSurround==1) { paramsArray <- paramsArray[c(10,11,12,13,15,16)] } # x, y, sigmaPos, sigmaNeg, intercept, slope
    if (flagSurround==0) { paramsArray <- paramsArray[c(10,11,12,15,16)] } # x, y, sigmaPos, intercept, slope
  }
  if (fitIntercept==0) {
    if (flagSurround==1) { paramsArray <- paramsArray[c(10,11,12,13,16)] } # x, y, sigmaPos, sigmaNeg, slope
    if (flagSurround==0) { paramsArray <- paramsArray[c(10,11,12,16)] } # x, y, sigmaPos, slope
  }
  hrfFun <- hrf
  plot( observedData, type='p', pch=3, bty='n', las=1, xlab='TR number', ylab='detrended BOLD', col='gray50' )
  lines( nLinNls(paramsArray), col='blue', lwd=2 )
}
dev.copy2pdf( file=sprintf( '%s_non_linear_estimates_examples.pdf', outSuffix ), width=12, height=8 ) 
dev.off()


