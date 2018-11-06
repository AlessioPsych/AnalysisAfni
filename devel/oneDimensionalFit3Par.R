args <- commandArgs(T)
print( args )

setwd('/home/fracasso/data/HighRes/AF_HighRes_04112016')
args <- c('epiMask.nii.gz','meanTsCat_det.nii.gz','par03','/packages/afni/17.0.13')

epiFile <- args[1]
detrendFile <- args[2]
outSuffix <- args[3]

source( sprintf('%s/AFNIio.R', args[4] ) )
library(neuRosim)
library(signal)

# get canonical hrf
print('get hrf...')
hrf <- canonicalHRF( seq(0,32,4), list(a1=6,a2=12,b1=0.9,b2=0.9,c=0.35) )
plot( hrf )

# build stimuli sequence
print('build stimuli...')
x <- seq(-10,10,0.01)
bWidth <- 0.5
bStep <- 0.25
nSteps <- 12
bStart <- 0.5
stim <- array( 0, c( nSteps, length(x) ) )
for (k in 1:12) {
  if (k == 1) {
    xStart <- bStart
    xEnd <- xStart + bWidth
    xSel <- x>=xStart & x<xEnd
  }
  if (k>1) {
    xStart <- xStart + bStep
    xEnd <- xEnd + bStep
    xSel <- x>=xStart & x<xEnd    
  }
  stim[k,xSel] <- 1 
}
baseline01 <- array( 0, c( 4, length(x) ) )
baseline <- array( 0, c( 5, length(x) ) )
rOrder <- 1:dim(stim)[1]
stimSeq <- rbind( baseline01,  # stack the baseline and the basic stimuli sequences 
                  stim,        # one after the other, revert the basic stimuli sequence  
                  baseline,    # if necessary
                  stim[ rev(rOrder), ],
                  baseline,
                  stim[ rev(rOrder), ],
                  baseline,
                  stim,
                  baseline)

plot( stimSeq[5,] ~ x, type='l', xlim=c(0,2))

# load the data
print('load data...')
ts <- read.AFNI( detrendFile )
tsArray <- array( ts$brk, c( prod( dim( ts$brk )[1:3] ), dim(ts$brk)[4] ) )
meanEpi <- read.AFNI( epiFile )


print('build prediction...')
centerArray <- seq( 0.4, 7, length.out=15 )
sigmaArrayPositive <- seq( 0.25, 4, length.out=15 )
sigmaArrayNegative <- sigmaArrayPositive
predictionGridTemp <- expand.grid( centerArray, sigmaArrayPositive, sigmaArrayNegative )
keepPredictionIdx <- predictionGridTemp[ ,2] < predictionGridTemp[ ,3]
predictionGrid <- predictionGridTemp[ keepPredictionIdx, ]

tsPrediction <- array( 0, c( dim(predictionGrid)[1], dim(stimSeq)[1] ) )
for (nPrediction in 1:dim(predictionGrid)[1] ) {
  prfPar <- as.numeric( predictionGrid[nPrediction,] )
  a <- dnorm( x, prfPar[1], prfPar[2] )
  b <- dnorm( x, prfPar[1], prfPar[3] )
  r <- a - b
  rMat <- matrix( rep( r, dim(stimSeq)[1] ), ncol=length(x), byrow=T )
  p <- apply( rMat * stimSeq, 1, sum )
  pConv <- conv( p, hrf )
  pConvTrim <- pConv[ 1 : length(p) ]
  tsPrediction[nPrediction,] <- pConvTrim
}
plot(tsPrediction[2,],type='l')

print('fitting...')
indexVol <- meanEpi$brk[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
nVoxels <- dim( tsTransposedAll )[2]
fitPart <- ceiling( seq( 1, nVoxels, length.out = 40 ) )
storeAllPred <- array(0,c(nVoxels,6))
storeAllExpectedTs <- array(0,dim(tsTransposedAll))

for (dataPart in 1:(length(fitPart)-1) ) {
  
  cat( paste( c('Iteration:',dataPart, 'out of:',length(fitPart), "" ) ) )
  print('...')
  tsTransposedFull <- tsTransposedAll[ , seq( fitPart[dataPart], fitPart[dataPart+1] ) ]
  arrayFull <- indexArray[ seq( fitPart[dataPart], fitPart[dataPart+1] ) ]
  #tsSd <- apply( tsTransposedFull, 2, sd  )
  
  idxFit <- arrayFull == 1
  tsTransposed <- tsTransposedFull[,idxFit]
  #tsMeans <- apply( tsTransposed, 2, mean  )
  #tsMeanScaled <- scale( tsTransposed, center=tsMeans, scale=FALSE  ) 
  #ssTot <- apply( tsMeans^2, 2, sum)
  
  tsTransposedMean <- apply( tsTransposed, 2, mean )
  ssTot <- apply( (tsTransposed-tsTransposedMean)^2, 2, sum)
  
  
  progress <- 0.05
  
  for (k in 1:dim(tsPrediction)[1] ) {
    
    dMat <- cbind( tsPrediction[k,] )
    dMat01 <- cbind( rep(1,length(dMat)), dMat )
    a <- solve( qr(dMat01), tsTransposed)
    #a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,tsTransposed) )
    
    
    #expectedTs <- dMat01%*%a
    expectedTs <- crossprod( t(dMat01), a )
    
    residualsSquared <- ( tsTransposed - expectedTs )^2
    ssRes <- apply( residualsSquared, 2, sum )
    r2 <- 1-(ssRes/ssTot)
    
    if (k==1) {
      storePred <- matrix( rep( as.numeric( predictionGrid[k,] ), dim(tsTransposed)[2] ), nrow = dim(tsTransposed)[2], ncol=dim(predictionGrid)[2], byrow=TRUE  )
      storePred <- cbind( storePred, r2, t( a ) )
      storePred[ storePred[ , dim(storePred)[2] ] < 0, c( (dim(storePred)[2]-1):dim(storePred)[2] ) ] <- 0 
      storeFit <- expectedTs
    }
    if (k>1) {
      updateIdx <- r2 > storePred[,4] & a[2,] > 0
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
  print('...')
}

progress <- 0.05
FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
xFWHM <- seq(-10,10,0.0125)
for (k in 1:dim(storeAllPred)[1]) {
  parameters <- storeAllPred[k,]
  a <- dnorm( xFWHM, 0, parameters[2] )
  b <- dnorm( xFWHM, 0, parameters[3] )
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
  
  if (k>dim(storeAllPred)[1]*progress) {
    cat(paste('*',""))
    progress <- progress + 0.05
  }
}

storeAllPredOut <- cbind( storeAllPred, FWHM )

fileTs <- sprintf('%s_PredixtedTs.nii.gz',outSuffix) 
fileParams <- sprintf('%s_params.nii.gz',outSuffix)
rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( fileTs, brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)
rParameters <- array( storeAllPredOut, c( dim(ts$brk)[1:3], 8 ) )
write.AFNI( fileParams, brk = rParameters,
            origin = meanEpi$origin, orient = meanEpi$orient,
            defhead = meanEpi$NI_head)

# add non linear fit only for certain PRFS

#idxPlot <- storeAllPredOut[,5] > 0.15
#plot( storeAllPredOut[idxPlot,8] ~ storeAllPredOut[idxPlot,1], pch='.' )
#abline( lm( storeAllPredOut[idxPlot,8] ~ storeAllPredOut[idxPlot,1] ) )
#summary( lm( storeAllPredOut[idxPlot,8] ~ storeAllPredOut[idxPlot,1] ) )

#hist( storePred[,4] )



