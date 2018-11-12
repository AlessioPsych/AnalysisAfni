args <- commandArgs(T)
print( args )

#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V3813Alessio_copy/prfModel')
#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V3814Serge')
#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V4162Ben')
#setwd('/home/fracasso/data/HighRes/AF_HighRes_04112016')

#args <- c('epi_mask.nii.gz','meanTs_det.nii.gz','output_sampling','/packages/afni/17.0.13')

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
centerArray <- seq( 0.4, 7, length.out=20 ) #15
sigmaArrayPositive <- seq( 0.75, 4, length.out=20 ) #15
sigmaArrayNegative <- sigmaArrayPositive
multArray <- c(1)
predictionGridTemp <- expand.grid( centerArray, sigmaArrayPositive, sigmaArrayNegative,multArray)
keepPredictionIdx <- predictionGridTemp[ ,2] < predictionGridTemp[ ,3]
predictionGrid <- predictionGridTemp[ keepPredictionIdx, ]

tsPrediction <- array( 0, c( dim(predictionGrid)[1], dim(stimSeq)[1] ) )
for (nPrediction in 1:dim(predictionGrid)[1] ) {
  prfPar <- as.numeric( predictionGrid[nPrediction,] )
  a <- dnorm( x, prfPar[1], prfPar[2] )
  b <- dnorm( x, prfPar[1], prfPar[3] )
  r <- a - b*prfPar[4]
  rMat <- matrix( rep( r, dim(stimSeq)[1] ), ncol=length(x), byrow=T )
  p <- apply( rMat * stimSeq, 1, sum )
  pConv <- conv( p, hrf )
  pConvTrim <- pConv[ 1 : length(p) ]
  tsPrediction[nPrediction,] <- pConvTrim
}
plot(tsPrediction[6,],type='l')

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
        storePred[ storePred[ , dim(storePred)[2] ] < 0, c( (dim(storePred)[2]-1):dim(storePred)[2] ) ] <- 0 
        storeFit <- expectedTs
      }
      if (k>1) {
        updateIdx <- r2 > storePred[,5] & a[2,] > 0
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
xFWHM <- seq(-7,7,0.0015)
for (k in 1:dim(storeAllPred)[1]) {
  if ( indexArray[k] == 1 ) {
    parameters <- storeAllPred[k,]
    a <- dnorm( xFWHM, 0, parameters[2] )
    b <- dnorm( xFWHM, 0, parameters[3] )
    r <- ( a - b*parameters[4] )
  
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
