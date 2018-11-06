#args <- commandArgs(T)
#print( args )

#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V3813Alessio_copy/prfModel')
#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V3814Serge')
#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V4162Ben')
rm(list=ls())
setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

args <- c('meanTsEye_res_mask.nii.gz','meanTsEye_res_detrend.nii.gz','output_meanTsEye_res_detrend','/home/alessiof/abin','directionFile.1D' ,'/home/alessiof/Programs/dlDevelopment/analysisAfni/generalPurpose')

epiFile <- args[1]
detrendFile <- args[2]
outSuffix <- args[3]

source( sprintf('%s/AFNIio.R', args[4] ) )
source( sprintf('%s/scaleData.R', args[6] ) )
library(neuRosim)
library(circular)

# get canonical hrf
print('get hrf...')
#hrf <- canonicalHRF( seq(0,25,0.2), list(a1=4.5,a2=7,b1=0.9,b2=0.9,c=0.35) )
hrf <- canonicalHRF( seq(0,25,0.2) )
plot( hrf~seq(0,25,0.2), bty='n', las=1 )

# direction file
directionFile <- read.table( file=args[5], header = FALSE, as.is=TRUE )
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


# load the data
print('load data...')
ts <- read.AFNI( detrendFile )
tsArray <- array( ts$brk, c( prod( dim( ts$brk )[1:3] ), dim(ts$brk)[4] ) )
meanEpi <- read.AFNI( epiFile )

fitSurround <- 0
debug <- 0

print('build prediction...')
if (fitSurround==0) {
  centerArray <- seq( -0.98*pi, 0.98*pi, length.out=20 ) 
  sigmaArrayPositive <- seq( 0.1, 6, length.out=20 ) 
  logParVonMises <- FALSE
  sigmaArrayNegative <- -1000
  multArray <- c(1)
  predictionGridTemp <- expand.grid( centerArray, sigmaArrayPositive, sigmaArrayNegative, multArray )
  keepPredictionIdx <- predictionGridTemp[ ,2] > predictionGridTemp[ ,3]
}
if (fitSurround==1) {
  centerArray <- seq( -0.95*pi, 0.95*pi, length.out=15 ) 
  sigmaArrayPositive <- seq( 0.1, 0.95, length.out=15 ) 
  sigmaArrayNegative <- seq( 0.01, 0.90, length.out=15 ) 
  multArray <- c(1)
  predictionGridTemp <- expand.grid( centerArray, sigmaArrayPositive, sigmaArrayNegative, multArray )
  keepPredictionIdx <- predictionGridTemp[ ,2] > predictionGridTemp[ ,3]
}
if (fitSurround==2) {
  centerArray <- seq( -0.95*pi, 0.95*pi, length.out=30 ) 
  sigmaArrayPositive <- seq( 0.05, 0.95, length.out=20 ) 
  sigmaArrayNegative <- seq( 0.1, 0.9, length.out=8 ) 
  multArray <- c(1)
  predictionGridTemp <- expand.grid( centerArray, sigmaArrayPositive, sigmaArrayNegative, multArray )
  keepPredictionIdx <- rep(TRUE,dim(predictionGridTemp)[1])
}
predictionGrid <- predictionGridTemp[ keepPredictionIdx, ]
tsPrediction <- array( 0, c( dim(predictionGrid)[1], length( seq( 1, dim(ts$brk)[4] ) ) ) )
progress <- 0.05
for (nPrediction in 1:dim(predictionGrid)[1] ) {
  
  if (nPrediction>dim(predictionGrid)[1]*progress) {
    cat(paste('*',""))
    progress <- progress + 0.05
  }
  
  if (fitSurround==0) {
    prfPar <- as.numeric( predictionGrid[nPrediction,] )
    muRad <- prfPar[1]
    kappaPar <- prfPar[2]
    stimuliSequenceTr <- 200 #in milliseconds
    sampleRepetitions <- directionFile[,3] / stimuliSequenceTr
    baselines <- ifelse( sampleRepetitions>10, 1, 0 )
    directionArray <- directionPol[,1]
    directionArray_expanded_radiants <- rep( directionArray, sampleRepetitions )
    baselines_expanded <- rep( baselines, sampleRepetitions )
    baselines_expanded <- rep(1,length(baselines_expanded))
    baselines_expanded[ saccStart ] <- 0
    
    firstStepTs_orig <- suppressWarnings( dvonmises( directionArray_expanded_radiants, muRad, kappaPar, log=logParVonMises  ) )
    firstStepTs_orig_baseline <- abs(baselines_expanded-1) * firstStepTs_orig
    
    convTs <- conv( firstStepTs_orig_baseline, hrf )
    convTsTrim <- convTs[1:length(firstStepTs_orig_baseline)]
    timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr/1000, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
    iTs <- interp1( timeArray, convTsTrim, seq(0,max(timeArray),length.out = dim(ts$brk)[4] ) )
    tsPrediction[nPrediction,] <- iTs
    
    if (debug==1) {
      par(mfrow=c(2,4))
      muRad01 <- -pi*0.5
      kappaPar01 <- 10
      muRad02 <- -pi*0.5
      kappaPar02 <- 0.1
      xPlot <- seq(-pi,pi,length.out = 200)
      plot( suppressWarnings( dvonmises( xPlot, muRad01, kappaPar01  ) ) ~ xPlot, type='l', bty='n', ylab='density', xlab='angle' )
      plot( suppressWarnings( dvonmises( xPlot, muRad02, kappaPar02  ) ) ~ xPlot, type='l', bty='n', ylab='density', xlab='angle' )
      plot( scaleData( suppressWarnings( dvonmises( xPlot, muRad01, kappaPar01  ) ), 1, 0 ) ~ xPlot, type='l', bty='n', ylab='density', xlab='angle' )
      plot( scaleData( suppressWarnings( dvonmises( xPlot, muRad02, kappaPar02  ) ), 1, 0 ) ~ xPlot, type='l', bty='n', ylab='density', xlab='angle' )
      muGauss01 <- -2
      sigmaGauss01 <- 1
      muGauss02 <- -2
      sigmaGauss02 <- 2.5
      xPlotGaussian <- seq(-5,5,length.out = 200)
      plot( dnorm( xPlotGaussian, muGauss01, sigmaGauss01 ) ~ xPlotGaussian, type='l', bty='n', ylab='density', xlab='dva'  ) 
      plot( dnorm( xPlotGaussian, muGauss02, sigmaGauss02 ) ~ xPlotGaussian, type='l', bty='n', ylab='density', xlab='dva'  ) 
      plot( scaleData( dnorm( xPlotGaussian, muGauss01, sigmaGauss01 ), 1, 0 ) ~ xPlotGaussian, type='l', bty='n', ylab='density', xlab='dva'  ) 
      plot( scaleData( dnorm( xPlotGaussian, muGauss02, sigmaGauss02 ), 1, 0 ) ~ xPlotGaussian, type='l', bty='n', ylab='density', xlab='dva'  ) 
      
      par(mfrow=c(3,2))
      piSpace <- seq(-pi*0.95, +pi*0.95, length.out = 200)
      plot( baselines_expanded * directionArray_expanded_radiants, type='l' )
      plot( suppressWarnings( dvonmises( piSpace, muRad, kappaPar, log=logParVonMises  ) ) ~ piSpace )
      plot( firstStepTs_orig, type='l' )
      #plot( firstStepTs_orig_inverted, type='l' )
      plot( firstStepTs_orig_baseline, type='l' )
      plot( convTsTrim~timeArray )
      lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), iTs, col='red' )
    }
  }
  if (fitSurround==1) {
    prfPar <- as.numeric( predictionGrid[nPrediction,] )
    muRad <- prfPar[1]
    rhoParPositive <- prfPar[2]
    rhoParNegative <- prfPar[3]
    stimuliSequenceTr <- 200 #in milliseconds
    sampleRepetitions <- directionFile[,3] / stimuliSequenceTr
    baselines <- ifelse( sampleRepetitions>10, 1, 0 )
    directionArray <- directionPol[,1]
    directionArray_expanded_radiants <- rep( directionArray, sampleRepetitions )
    baselines_expanded <- rep( baselines, sampleRepetitions )
    baselines_expanded <- rep(1,length(baselines_expanded))
    baselines_expanded[ saccStart ] <- 0
    
    piSpace <- seq(-pi*0.95, +pi*0.95, length.out = 200)
    positivePrf <- suppressWarnings( dwrappednormal( piSpace, muRad, rhoParPositive  ) )
    negativePrf <- -1*suppressWarnings( dwrappednormal( piSpace, muRad, rhoParNegative  ) )
    sumPrf <- positivePrf + negativePrf
    sumPrfFit <- smooth.spline( piSpace, sumPrf )
    firstStepTs_orig_prf <- predict( sumPrfFit, directionArray_expanded_radiants )
    firstStepTs_orig_prf_baseline <- abs(baselines_expanded-1) * firstStepTs_orig_prf$y
      
    #firstStepTs_orig_center <- suppressWarnings( dwrappednormal( directionArray_expanded_radiants, muRad, rhoPar01  ) )
    #firstStepTs_orig_surround <- suppressWarnings( dwrappednormal( directionArray_expanded_radiants, muRad, rhoPar02  ) )
    #firstStepTs_orig_inverted_baseline <- abs(baselines_expanded-1) * firstStepTs_orig_inverted
    
    convTs <- conv( firstStepTs_orig_prf_baseline, hrf )
    convTsTrim <- convTs[1:length(firstStepTs_orig_prf_baseline)]
    timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr/1000, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
    iTs <- interp1( timeArray, convTsTrim, seq(0,max(timeArray),length.out = dim(ts$brk)[4] ) )
    tsPrediction[nPrediction,] <- iTs
    
    if (debug==1) {
      par(mfrow=c(3,2))
      plot( positivePrf ~ piSpace, type='l', col='red', ylim=c(min(negativePrf),max(positivePrf)) ) 
      lines( piSpace, negativePrf, col='blue' )
      plot( sumPrf ~ piSpace, type='l', col='darkorange', lwd=4 ) 
      lines( sumPrfFit$x, sumPrfFit$y, lty=2, col='black' )
      plot( firstStepTs_orig_prf$y, type='l' )
      plot( firstStepTs_orig_prf_baseline, type='l' )
      plot( convTsTrim~timeArray )
      lines( seq(0,max(timeArray),length.out = dim(ts$brk)[4] ), iTs, col='red' )
    }      
  }
  if (fitSurround==2) {
    prfPar <- as.numeric( predictionGrid[nPrediction,] )
    muRad <- prfPar[1]
    rhoPar <- prfPar[2]
    suppressionRatio <- prfPar[3]
    logParVonMises <- TRUE
    stimuliSequenceTr <- 200 #in milliseconds
    sampleRepetitions <- directionFile[,3] / stimuliSequenceTr
    baselines <- ifelse( sampleRepetitions>10, 1, 0 )
    directionArray <- directionPol[,1]
    directionArray_expanded_radiants <- rep( directionArray, sampleRepetitions )
    baselines_expanded <- rep( baselines, sampleRepetitions )
    baselines_expanded <- rep(1,length(baselines_expanded))
    baselines_expanded[ saccStart ] <- 0

    piSpace <- seq(-pi*0.95, +pi*0.95, length.out = 200)
    positivePrf <- suppressWarnings( dwrappednormal( piSpace, muRad, rhoPar  ) )
    suppressionConstant <- diff( range( positivePrf ) )*suppressionRatio
    positivePrf <- positivePrf - suppressionConstant
    positivePrfFit <- smooth.spline( piSpace, positivePrf )
    firstStepTs_orig_prf <- predict( positivePrfFit, directionArray_expanded_radiants )
    firstStepTs_orig_prf_baseline <- abs(baselines_expanded-1) * firstStepTs_orig_prf$y
    
    convTs <- conv( firstStepTs_orig_prf_baseline, hrf )
    convTsTrim <- convTs[1:length(firstStepTs_orig_prf_baseline)]
    timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr/1000, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
    timeArrayInterp <- seq(0,max(timeArray),length.out = dim(ts$brk)[4] )
    iTs <- interp1( timeArray, convTsTrim, timeArrayInterp )
    tsPrediction[nPrediction,] <- iTs
    
    if (debug==1) {
      par(mfrow=c(3,2))
      plot( positivePrf ~ piSpace, type='l', col='red' ) 
      lines( positivePrfFit$x, positivePrfFit$y, lty=2, col='black' )
      plot( firstStepTs_orig_prf$y, type='l' )
      plot( firstStepTs_orig_prf_baseline, type='l' )
      plot( convTsTrim~timeArray )
      lines( timeArrayInterp, iTs, col='red' )
    }
  }
  if (fitSurround==3) {
    prfPar <- as.numeric( predictionGrid[nPrediction,] )
    muRad <- -pi*0.8 #prfPar[1]
    rhoPar <- 0.8 #prfPar[2]
    suppressionRatio <- 0 #prfPar[3]
    motorUncertainty <- 0.8
    logParVonMises <- TRUE
    stimuliSequenceTr <- 200 #in milliseconds
    sampleRepetitions <- directionFile[,3] / stimuliSequenceTr
    baselines <- ifelse( sampleRepetitions>10, 1, 0 )
    directionArray <- directionPol[,1]
    directionArray_expanded_radiants <- rep( directionArray, sampleRepetitions )
    baselines_expanded <- rep( baselines, sampleRepetitions )
    baselines_expanded <- rep(1,length(baselines_expanded))
    baselines_expanded[ saccStart ] <- 0
    
    piSpace <- seq(-pi*0.95, +pi*0.95, length.out = 200)
    #positivePrf <- suppressWarnings( dwrappednormal( piSpace, muRad, rhoPar  ) )
    #suppressionConstant <- diff( range( positivePrf ) )*suppressionRatio
    #positivePrf <- positivePrf - suppressionConstant
    #positivePrfFit <- smooth.spline( piSpace, positivePrf )
    
    positivePrf <- suppressWarnings( dwrappednormal( piSpace, muRad, rhoPar  ) )
    storeMotorFields <- array( 0, c( length(directionArray_expanded_radiants), length(piSpace)   ) )
    storeTs <- rep( 0, length(directionArray_expanded_radiants), 1   )
    for (motorCounter in 1:length(directionArray_expanded_radiants)) {
      storeMotorFields[motorCounter, ] <- suppressWarnings( dwrappednormal( piSpace, directionArray_expanded_radiants[motorCounter], motorUncertainty  ) )
      storeTs[motorCounter] <- sum( positivePrf*storeMotorFields[motorCounter, ] )
    }
    
    firstStepTs_orig_prf_baseline <- abs(baselines_expanded-1) * storeTs
    
    convTs <- conv( firstStepTs_orig_prf_baseline, hrf )
    convTsTrim <- convTs[1:length(firstStepTs_orig_prf_baseline)]
    timeArray <- seq( 0, length(convTsTrim)*stimuliSequenceTr/1000, length.out = length(convTsTrim) )#  *stimuliSequenceTr/1000
    timeArrayInterp <- seq(0,max(timeArray),length.out = dim(ts$brk)[4] )
    iTs <- interp1( timeArray, convTsTrim, timeArrayInterp )
    plot(iTs, type='l')
    tsPrediction[nPrediction,] <- iTs
    
    if (debug==1) {
      par(mfrow=c(3,2))
      plot( positivePrf ~ piSpace, type='l', col='red' ) 
      lines( positivePrfFit$x, positivePrfFit$y, lty=2, col='black' )
      plot( firstStepTs_orig_prf$y, type='l' )
      plot( firstStepTs_orig_prf_baseline, type='l' )
      plot( convTsTrim~timeArray )
      lines( timeArrayInterp, iTs, col='red' )
    }
  }  
}
plot(tsPrediction[25,],type='l')

#source( sprintf('/data1/projects/myelin/analysisAfni/powerSpectra.R' ) )
#plot(pConvTrim~seq(1,length(pConvTrim)),type='l',col='black')
#lines(iTs~seq(1,length(pConvTrim),by=2),type='l',col='red')
#outFft <- powerSpectra(iTs,2)
#plot( outFft$amp~outFft$freq, type='l' )

print('fitting...')
indexVol <- meanEpi$brk[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
nVoxels <- dim( tsTransposedAll )[2]
fitPart <- ceiling( seq( 1, nVoxels, length.out = 80 ) )
storeAllPred <- array(0,c(nVoxels,7))
storeAllExpectedTs <- array(0,dim(tsTransposedAll))

for (dataPart in 1:(length(fitPart)-1) ) {
  
  cat( paste( c('Iteration:',dataPart, 'out of:',length(fitPart), "" ) ) )
  print('...')
  tsTransposedFull <- tsTransposedAll[ , seq( fitPart[dataPart], fitPart[dataPart+1] ) ]
  arrayFull <- indexArray[ seq( fitPart[dataPart], fitPart[dataPart+1] ) ]
  #tsSd <- apply( tsTransposedFull, 2, sd  )
  
  idxFit <- arrayFull == 1
  if (sum(idxFit)>1) {
    tsTransposed <- tsTransposedFull[,idxFit]
    tsMeans <- apply( tsTransposed, 2, mean  )
    #tsMeanScaled <- scale( tsTransposed, center=tsMeans, scale=FALSE  ) 
    #ssTot <- apply( tsMeanScaled^2, 2, sum)
    ssTot <- apply( (tsTransposed-tsMeans)^2, 2, sum)
    progress <- 0.05
    
    for (k in 1:dim(tsPrediction)[1] ) {
      
      dMat <- cbind( tsPrediction[k,] )
      dMat01 <- cbind( rep(1,length(dMat)), dMat )
      #a <- solve( qr(dMat01), tsTransposed)
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,tsTransposed) )
      
      #expectedTs <- dMat01%*%a
      expectedTs <- crossprod( t(dMat01), a )
      residualsSquared <- ( tsTransposed - expectedTs )^2
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot)
      
      if (k==1) {
        storePred <- matrix( rep( as.numeric( predictionGrid[k,] ), dim(tsTransposed)[2] ), nrow = dim(tsTransposed)[2], ncol=dim(predictionGrid)[2], byrow=TRUE  )
        storePred <- cbind( storePred, r2, t( a ) )
        storePred[ storePred[ , dim(storePred)[2] ] < 0 , c( (dim(storePred)[2]-2):dim(storePred)[2] ) ] <- 0 
        storeFit <- expectedTs
      }
      if (k>1) {
        updateIdx <- r2 > storePred[,5] & a[2,] > 0
        #updateIdx <- r2 > storePred[,5]
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

progress <- 0.05
FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
#xFWHM <- seq(-7,7,0.0015)
#for (k in 1:dim(storeAllPred)[1]) {
#  if ( indexArray[k] == 1 ) {
#    parameters <- storeAllPred[k,]
#    a <- dnorm( xFWHM, 0, parameters[2] )
#    b <- dnorm( xFWHM, 0, parameters[3] )
#    r <- ( a - b*parameters[4] )

#    maxIdx <- which.max(r)
#    rHalf <- r[ maxIdx:length(xFWHM) ]
#    xHalf <- xFWHM[ maxIdx:length(xFWHM) ]
#    halfMax <- max(rHalf)/2
#    xMid <- xHalf[ rHalf<=halfMax  ]  
#    FWHMcenter <- min( xMid )*2 
#    rMin <- which.min(rHalf)
#    FWHMsurround <- xHalf[rMin[1]]*2
#    FWHM[k,] <- c(FWHMcenter,FWHMsurround)
#  }
#  if (k>dim(storeAllPred)[1]*progress) {
#    cat(paste('*',""))
#    progress <- progress + 0.05
#  }
#}

storeAllPredOut <- cbind( storeAllPred, FWHM )

fileParams <- sprintf('%s_params.nii.gz',outSuffix)
fileTs <- sprintf('%s_Ts.nii.gz',outSuffix)

rParameters <- array( storeAllPredOut, c( dim(ts$brk)[1:3], dim(storeAllPredOut)[2] ) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', epiFile )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( '__tt_parameters.nii.gz', brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')
instr <- sprintf( 'mv __tt_parameters.nii.gz %s', fileParams )
system( instr)

rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( '__tt_expectedTs.nii.gz', brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)
instr <- sprintf( 'mv __tt_expectedTs.nii.gz %s', fileTs )
system( instr)

#predictionGridTemp <- expand.grid( centerArray, sigmaArrayPositive, sigmaArrayNegative,multArray)
labels <- c('ecc','sigmaPos','sigmaNeg','const','var_exp','intercept','slope','FWHM_center','FWHM_surr')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

#instr <- sprintf( '3dNwarpApply -master coregistration/targetVolume.nii.gz -source output01_params.nii.gz -nwarp \u0027coregistration/tMat.1D topUpDir/warpTop_PLUS_WARP+orig\u0027 -prefix results/prfModel_R.nii.gz -interp NN' );
#system(instr)

#instr <- sprintf('3dAllineate -final NN -master coregistration2/targetVolume_res.nii.gz -source prfModel/outParamsAll.nii.gz -1Dmatrix_apply coregistration2/tMat.1D -prefix results/params_res_all.nii.gz')
#system(instr)

#instr <- sprintf('3dAllineate -master results/anatomy_res.nii.gz -source outParamsAll_V2.nii.gz -1Dmatrix_apply coregistrationV2/tMat.1D -prefix results/params_res_V2.nii.gz -final NN')
#system(instr)

#instr <- sprintf('3dAllineate -master results/anatomy_res.nii.gz -source outParamsAll_80.nii.gz -1Dmatrix_apply coregistration/tMat.1D -prefix results/params_res_80.nii.gz -final NN')
#system(instr)

# 
# #ecc <- storeAllPredOut[,1] 
# #sigma01 <- storeAllPredOut[,8]
# #r2 <- storeAllPredOut[,5]
# #eccFilt <- ecc[ r2>0.25 ]
# #sigma01Filt <- sigma01[ r2 > 0.25 ]
# #plot( sigma01Filt ~ eccFilt, pch='.' )
# #summary(lm( sigma01Filt ~ eccFilt ))
# #abline( lm( sigma01Filt ~ eccFilt ) )
# 
# # nlin fit
# options(warn=-1)
# nLinNls <- function( center, sigmaPositive, sigmaNegative, multFact, intercept, slope ) {
#   a <- dnorm( x, center, sigmaPositive )
#   b <- dnorm( x, center, sigmaNegative )
#   r <- ( a - b*multFact )
#   rMat <- matrix( rep( r, dim(stimSeq)[1] ), ncol=length(x), byrow=T )
#   p <- apply( rMat * stimSeq, 1, sum )
#   pConv <- conv( p, hrf )
#   pConvTrim <- intercept + pConv[ 1 : length(p) ] * slope
#   return(pConvTrim)
# }
# nLin <- function( parameters, xInput, hrfInput, stimSeqInput ) {
#   a <- dnorm( xInput, parameters[1], parameters[2] )
#   b <- dnorm( xInput, parameters[1], parameters[3] )
#   r <- ( a - b*parameters[4] )
#   rMat <- matrix( rep( r, dim(stimSeqInput)[1] ), ncol=length(x), byrow=T )
#   p <- apply( rMat * stimSeqInput, 1, sum )
#   pConv <- conv( p, hrf )
#   pConvTrim <- parameters[5] + pConv[ 1 : length(p) ] * parameters[6]
#   return(pConvTrim)
# }
# nLinFit <- function( parameters, xInput, yInput, hrfInput, stimSeqInput ) {
#   pConvTrim <- nLin( parameters, xInput, hrfInput, stimSeqInput )
#   out <- sum(( pConvTrim - yInput )^2) 
#   return(out)
# }
# 
# #center <- parameters[1]
# #sigmaPositive <- parameters[2]
# #sigmaNegative <- parameters[3]
# #multFact <- parameters[4]
# #intercept <- parameters[5]
# #slope <- parameters[6]
# #nLinNls( center, sigmaPositive, sigmaNegative, multFact, intercept, slope )
# nLin( parameters, x, hrf, stimSeq )
# nLinPars <- array(0, c( dim(storeAllPred)[1], 8 ) )
# keepPredictors <- apply( storeAllPredOut == 0, 1, sum )
# iii <- which( indexArray == 1 & storeAllPredOut[,5] > 0.15 & keepPredictors==0 )
# options(warn=-1)
# for (k in 1:40) { #length(iii)
#   parameters <- storeAllPred[ iii[k], c(1,2,3,4,6,7) ]
#   tsToFit <- tsTransposedAll[ ,iii[k] ] 
#   if (parameters[5] < 0) {
#     lowBound <- c( parameters[1]*0.5, parameters[2]*0.5, parameters[3]*0.98, parameters[4]*0.5, parameters[5]*1.5, parameters[6]*0.5 )
#     upBound <- c( parameters[1]*1.5, parameters[3]*0.95, parameters[3]*1.5, parameters[4]*1.5, parameters[5]*0.5, parameters[6]*1.5 )
#   }
#   if (parameters[5] > 0) {
#     lowBound <- c( parameters[1]*0.5, parameters[2]*0.5, parameters[3]*0.98, parameters[4]*0.5, parameters[5]*0.5, parameters[6]*0.5 )
#     upBound <- c( parameters[1]*1.5, parameters[3]*0.95, parameters[3]*1.5, parameters[4]*1.5, parameters[5]*1.5, parameters[6]*1.5 )
#   }
#   #optOut01 <- optim( par=parameters, fn=nLinFit,
#   #                 xInput = x, yInput = tsToFit, hrfInput = hrf, stimSeqInput=stimSeq,
#   #                 lower=lowBound, upper=upBound, method="L-BFGS-B", control=list( factr=1e-9, maxit=100))
#   
#   centerS <- parameters[1]
#   sigmaPositiveS <- parameters[2]
#   sigmaNegativeS <- parameters[3]
#   multFactS <- parameters[4]
#   interceptS <- parameters[5]
#   slopeS <- parameters[6]
#   outNls <- tryCatch( nls( tsToFit ~ nLinNls( center, sigmaPositive, sigmaNegative, multFact, intercept, slope ),
#                 start = list( center=centerS, sigmaPositive=sigmaPositiveS, 
#                      sigmaNegative=sigmaNegativeS, multFact=multFactS,
#                      intercept=interceptS, slope=slopeS),
#                 lower = lowBound,
#                 upper = upBound,
#                 algorithm = 'port',
#                 nls.control(maxiter = 10, tol = 1e-02, minFactor = 1/1024,
#                             printEval = FALSE, warnOnly = TRUE) ), 
#                 error = function(e) 'error' );
#   if (outNls!='error') {
#     r2Nlin <- summary( lm( nLin(   coef( outNls ), x, hrf, stimSeq ) ~ tsToFit ) )$r.squared
#   }
#     #optOut01 <- optim( par=parameters, fn=nLinFit,
#   #                   xInput = x, yInput = tsToFit, hrfInput = hrf, stimSeqInput=stimSeq,
#   #                   method='CG', control=list(ndeps=c( 0.03, 0.03, 0.03, 0.1, 200, 200 ) , reltol=1e-2, maxit=10 ) )
#   #r2Nlin <- summary( lm( nLin( optOut01$par, x, hrf, stimSeq ) ~ tsToFit ) )$r.squared
#   nLinPars[ iii[k], ] <- c( optOut01$par, r2Nlin, storeAllPred[ iii[k],5] )
#   
# }
# options(warn=0)
# nLinPars[iii,]
# ####warnings off
# #plot( tsToFit, type='l' )
# #lines( nLin( optOut$par, x, hrf, stimSeq ), col='red')
# #optOut$par
# #lines( nLin( optOut01$par, x, hrf, stimSeq ), col='green')
# #optOut01$par
# #lines( nLin( parameters, x, hrf, stimSeq ), col='blue')
# #parameters
# 
# #summary( lm( nLin( optOut$par, x, hrf, stimSeq ) ~ tsToFit ) )$r.squared
# #summary( lm( nLin( parameters, x, hrf, stimSeq ) ~ tsToFit ) )$r.squared
# 
# 
# 
# 
# 
# # add non-linear fit step with optim, only for thresholded data
# 
# #idxPlot <- storeAllPredOut[,5] > 0.15
# #plot( storeAllPredOut[idxPlot,8] ~ storeAllPredOut[idxPlot,1], pch='.' )
# #abline( lm( storeAllPredOut[idxPlot,8] ~ storeAllPredOut[idxPlot,1] ) )
# #summary( lm( storeAllPredOut[idxPlot,8] ~ storeAllPredOut[idxPlot,1] ) )
# 
# #hist( storePred[,4] )
# 
# 
# 
