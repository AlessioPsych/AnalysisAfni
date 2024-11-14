args <- commandArgs(T)
print( args )
mainDir <- getwd()
setwd( mainDir )

# to debug
#rm( list=ls() ); gc();
#setwd('/analyse/Project0226/dataToSort/eyeMovPoster/ASM16_BarsEyes')
#args <- c('zzz_meanTs_bars_topUp_dt_blur_sm.nii.gz', '4', '0', '0', '0', 'zzz_del_prfSimple')
#mainDir <- getwd()
#setwd( mainDir )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( neuRosim )
library( parallel )

# arrange inputs
inputFile <- args[1]
stimType <- as.numeric( args[2] )
fineFit <-  as.numeric( args[3] )
prfType <-   as.numeric( args[4] )
flagSurround <-   as.numeric( args[5] )
outDirectory <- args[6]

print( sprintf('input file = %s', inputFile ) )
print( sprintf('stim type = %1.0f', stimType ) )
print( sprintf('fineFit = %1.0f', fineFit ) )
print( sprintf('prfType = %1.0f', prfType ) )
print( sprintf('flagSurround = %1.0f', flagSurround ) )
print( sprintf('output directory = %s', outDirectory ) )

# create output directory
system( sprintf('mkdir %s', outDirectory) )

# get TR
instr <- sprintf('3dinfo -tr %s > __tttt.1D', inputFile); system( instr )
trValue <- read.table('__tttt.1D'); system('rm __tttt.1D')
trValue <- as.numeric( trValue )

# get number of dynamics
instr <- sprintf('3dinfo -nv %s > __tttt.1D', inputFile); system( instr )
nDynamics <- read.table('__tttt.1D'); system('rm __tttt.1D')
nDynamics <- as.numeric( nDynamics )

# load stimuli definition
print('get stimulus...')
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

#x11( width=3, height=3 )
#for ( snap in 1:dim(stimMat)[3] ) {
#  image( stimMatFlip[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01)
#}
stimSeq <- stimMatFlip
x <- seq(-10,10,length.out = dim(stimSeq)[1] )
y <- seq(-5.5,5.5,length.out = dim(stimSeq)[2] )
stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
rm( stimMat )
rm( stimMatFlip )
gc()

#this part of the code builds a matrix with all the possible prediction tested, for both models at this stage
if (fineFit==0) {
  xElements <- 22 #22
  yElements <- 24 #24
  sigmaArrayPositiveElements <- 25 #32
  multParElements <- 0 #5
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  s_nLinElements <- 3
  n_cyclesElements <- 1
  thetaElements <- 5
  phaseElements <- 4  
}
if (fineFit==1) {
  xElements <- 12
  yElements <- 12
  sigmaArrayPositiveElements <- 14
  multParElements <- 6
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  s_nLinElements <- 3
  n_cyclesElements <- 1
  thetaElements <- 7
  phaseElements <- 5  
}

print('build prediction...')
xPosFit <- seq( -10, 10, length.out=xElements )
yPosFit <- seq( -5.5, 5.5, length.out=yElements )
#sigmaArrayPositive <- seq( 0.5, 11, length.out=sigmaArrayPositiveElements )
sigmaArrayPositive <- exp( seq( log( 0.2 ), log( 11 ), length.out=sigmaArrayPositiveElements ) )
#sigmaArrayPositive <- exp( seq( log( 0.2 ), log( 8 ), length.out=sigmaArrayPositiveElements ) )
#staticNonLin <- seq( 0.3, 1.8, length.out=s_nLinElements )
staticNonLin <- 1
if (flagSurround==1) { sigmaArrayNegative <- exp( c( log(1), seq( log(1.05), log(3.5), length.out = sigmaArrayPositiveElements) ) ) } # 1.05 2.5
#if (flagSurround==1) { sigmaArrayNegative <- c( 1, seq( 1.05, 3.5, length.out = 15) ) } # 1.05 2.5
#if (flagSurround==1) { sigmaArrayNegative <- seq( 1.05, 4, length.out = 5) }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
par_hrf_a1 <- seq( 6, 9, length.out=hrfDelayOnsetElements )
par_hrf_a2 <- seq( 12, 15, length.out=hrfDelayUnderShootElements )
if (flagSurround==1) { multPar <- 0 }#{ multPar <- seq(0,0.8, length.out = multParElements) }
if (flagSurround==0) { multPar <- 0 }
if (prfType==1){
  #n_cycles <- seq( 0.20, 0.8, length.out = n_cyclesElements )
  n_cycles <- 0.25 #hard-code different values to check effect of SF content
  #angle <- c(5, 45, 90, 135, 175)
  angle <- seq(5, 175, length.out = thetaElements)
  #angle <- seq(0, 180, length.out = thetaElements)
  #angle <- 0 #hard-code different values to check effect of orientation
  phase <- seq( 0.1, 1, length.out = phaseElements)
  #phase <- 0
}
if (prfType==0) { 
  n_cycles <- 0
  angle <- 0
  phase <- 0
}

# generate the global prediction grid
#predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, staticNonLin, n_cycles, angle, phase )
#keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
#predictionGridGlobal <- predictionGridTemp[ keepPredictionIdx, ]
predictionGridGlobal <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative, par_hrf_a1, par_hrf_a2, multPar, staticNonLin, n_cycles, angle, phase )
predictionGridGlobal[,4] <- predictionGridGlobal[,3]*predictionGridGlobal[,4]
print( sprintf('total number of predictions to generate: %1.0f', dim( predictionGridGlobal )[1] ) )

# this part builds a matrix with the starting and ending prediction index to fit for each loop in the next for loop
if (dim( predictionGridGlobal )[1] <= 900 ) {
  totalIterations <- 3
}
if (dim( predictionGridGlobal )[1] > 900 ) {
  totalIterations <- ceil( dim( predictionGridGlobal )[1]/500 )
}

limitsPrediction <- round( seq(1,dim(predictionGridGlobal)[1], length.out=totalIterations ) ) 
limitsPrediction[1] <- 1
limitsPrediction[ length(limitsPrediction) ] <- dim(predictionGridGlobal)[1]
limitsPredictionMatrix <- array( 0, c( (length(limitsPrediction)-1) , 2 ) )
limitsPredictionMatrix[,1] <- limitsPrediction[1:(length(limitsPrediction)-1)]
limitsPredictionMatrix[,2] <- limitsPrediction[2:(length(limitsPrediction))]
limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] <- limitsPredictionMatrix[2:dim(limitsPredictionMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk to be fitted
runIndexPredictions <- seq( 1:dim(limitsPredictionMatrix)[1] )
print( sprintf('total number of iterations: %1.0f', length( runIndexPredictions ) ) )

# MRI time and stimuli time
mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
timeStimuli <- seq( 0, max(mriTime)+as.numeric( trValue[1] )/2, length.out = dim(stimSeq)[3] )
samplingTime <- as.numeric( nDynamics[1]*trValue[1] )/dim(stimSeq)[3]

#### function to generate the predictions ####
generatePrediction <- function( indexPrediction, inputPredictionGrid, samplingTime, mriTime, timeStimuli, prfType, flagSurround, xSpace, ySpace, scaleData, stimSeqMat ) {
  
  library( neuRosim )
  library( pracma )
  x <- xSpace
  y <- ySpace
  prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
  
  #prfPar <- c( 0, 0, 1.5, 1.65, 6, 12, 0, 1, 0.25, 10+180, 0); flagSurround <- 0 # test ts shape without multpar
  
  hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE ); #plot( hrf ~ seq(0,30, samplingTime ) )
  
  if (prfType==0) { # simple PRF or CSS
    a <- dnorm( x, prfPar[1], prfPar[3] )
    b <- dnorm( y, prfPar[2], prfPar[3] )
    a[ x<qnorm( 0.001, mean=prfPar[1], sd=prfPar[3] ) | x>qnorm( 0.999, mean=prfPar[1], sd=prfPar[3] ) ] <- 0 #cut the tails x
    b[ y<qnorm( 0.001, mean=prfPar[2], sd=prfPar[3] ) | y>qnorm( 0.999, mean=prfPar[2], sd=prfPar[3] ) ] <- 0 #cut the tails y 
    imgCenter <- tcrossprod(a,b)
    #cut the tails, to plot:
    #par(mfrow=c(2,1))
    #plot( a ~ x, bty='n', las=1 ); abline( v=qnorm( 0.01, mean=prfPar[1], sd=prfPar[3] ), col='red', lty=2, lwd=1  ); abline( v=qnorm( 0.99, mean=prfPar[1], sd=prfPar[3] ), col='red', lty=2, lwd=1  )
    #plot( b ~ y, bty='n', las=1 ); abline( v=qnorm( 0.01, mean=prfPar[2], sd=prfPar[3] ), col='red', lty=2, lwd=1  ); abline( v=qnorm( 0.99, mean=prfPar[2], sd=prfPar[3] ), col='red', lty=2, lwd=1  )
    
    a <- dnorm( x, prfPar[1], prfPar[4] )
    b <- dnorm( y, prfPar[2], prfPar[4] )
    a[ x<qnorm( 0.001, mean=prfPar[1], sd=prfPar[4] ) | x>qnorm( 0.999, mean=prfPar[1], sd=prfPar[4] ) ] <- 0 #cut the tails x
    b[ y<qnorm( 0.001, mean=prfPar[2], sd=prfPar[4] ) | y>qnorm( 0.999, mean=prfPar[2], sd=prfPar[4] ) ] <- 0 #cut the tails y 
    imgSurround <- tcrossprod(a,b)
    
    if (flagSurround==1)
      r <- imgCenter - imgSurround
    else{
      r <- imgCenter
    }
    #image( r, col=grey.colors(1000) )
    rMat <- array(r)
  }
  
  if (prfType==1) { #gabor filter
    
    a <- length(x)
    b <- length(y)
    x_vecGab <- matrix(rep(x - prfPar[1], b), nrow = a) 
    y_vecGab <- matrix(rep(y - prfPar[2], each = a), nrow = a)
    
    #gabMat <- meshgrid((seq(min(x), max(x), length.out = length(x))),(seq(min(y), max(y), length.out = length(y))))
    #x_vecGab <- gabMat$X
    #y_vecGab <- gabMat$Y
    
    n_cycles <- prfPar[9]
    sigma_x <- prfPar[3]
    sigma_y <- (prfPar[3])/1
    freq <- n_cycles/sigma_x
    angle <- prfPar[10]
    theta <- angle*(pi/180)
    phase <- prfPar[11]
    
    #orientation definition
    x_theta <- x_vecGab * -cos(theta) + y_vecGab * sin(theta)
    #x_theta <- (x_vecGab * -cos(theta) + y_vecGab * sin(theta)) - prfPar[1]
    #y_theta <- (-x_vecGab * cos(theta) + y_vecGab * sin(theta)) - prfPar[2]
    
    #creating a gaussian envelope
    #a <- dnorm(x, prfPar[1], prfPar[3])
    #b <- dnorm(y, prfPar[2], prfPar[3])
    #envelope <- tcrossprod(a, b)
    envelope <- exp(-0.5 *(x_vecGab^2/sigma_x^2 + y_vecGab^2/sigma_y^2 ))
    
    
    grating <- cos((2*pi*freq*x_theta) + phase) #changed code to t() because dimensions weren't compatible ?
    gabor <- envelope*grating
    r <- gabor
    #image(gabor, col = hsv(seq(0.7, 1, 0.01)))
    
    rMat <- array(gabor)
  }
  
  if ( sum( is.na(rMat) ) == 0 ) {
    
    if (prfType==0) { #simple prf or CSS
      predictionLoop <- as.numeric( crossprod( rMat, stimSeqMat ) )^prfPar[8] #### this is the slow step, this is the reason for the parallel computing ####    
    }
    if (prfType==1) { #gabor filter model
      #predictionLoop <- Mod(as.complex(as.numeric( crossprod( rMat, stimSeqMat ) ) )^prfPar[8] ) #### this is the slow step, this is the reason for the parallel computing ####
      predictionLoop <- as.numeric( crossprod( rMat, stimSeqMat ) )^prfPar[8] #### this is the slow step, this is the reason for the parallel computing ####
    }
    
    pConv01 <- conv( predictionLoop, hrf )
    pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]; #plot(pConvTrim01,type='l')
    predictionLoop01Interp <- interp1( x=timeStimuli, y=pConvTrim01, xi=mriTime, method=c('nearest') )
    returnPrediction <- round( scaleData( predictionLoop01Interp, 1, 0 ), 5 ); #plot(returnPrediction,type='l')
    
  }
  
  #par(mfrow=c(2,2))
  #image( r, col=grey.colors(1000) )
  #plot( r[,75] ~ x, bty='n', las=1, type='l' )
  #plot( r[120,] ~ y, bty='n', las=1, type='l' )
  #plot( returnPrediction ~ mriTime, bty='n', las=1, type='l'  )#
  
  if ( sum( is.na(rMat) ) > 0 ) {
    returnPrediction <- rep(0,length(mriTime))
  }
  
  return( returnPrediction ) 
  
} #### scale predictions between 0 and 1 and round them to 5 digits

predictionWrapper <- function( passIdx, limitsMat, inputGrid ) {
  
  print( sprintf('iteration: %1.0f of %1.0f', passIdx, dim( limitsMat )[1] ) )
  selectedRows <- limitsMat[passIdx, ]
  
  # generate predictions:
  nCores <- 12
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( tsPredictions <- parSapply(cl, selectedRows[1]:selectedRows[2], generatePrediction,
                                                          inputPredictionGrid=inputGrid,
                                                          samplingTime=samplingTime,
                                                          mriTime=mriTime,
                                                          timeStimuli=timeStimuli,
                                                          prfType=prfType,
                                                          flagSurround=flagSurround,
                                                          xSpace=x,
                                                          ySpace=y,
                                                          scaleData=scaleData,
                                                          stimSeqMat=stimSeqMat ) )
  stopCluster(cl)
  showConnections()
  rm(cl); gc();
  print( storeTimePar )
  #storeTimePar <- system.time( tsPredictions <- sapply( selectedRows[1]:selectedRows[2], generatePrediction, inputPredictionGrid=inputGrid ) )
  #print( storeTimePar )
  
  
  # concatenate parameters and predictions:
  outMatrix <- rbind( t( predictionGridGlobal[ selectedRows[1]:selectedRows[2], ] ), tsPredictions )
  
  # save parameters and predictions:
  save( outMatrix, file = sprintf('%s/ttt_predictionFile_%1.0f.RData', outDirectory, passIdx ) )
  
  return(1)
}

sapply( 1:dim(limitsPredictionMatrix)[1], predictionWrapper, limitsMat = limitsPredictionMatrix, inputGrid = predictionGridGlobal)

rm( list=ls() ); gc()

# nCores <- 4
# cl <- makeCluster(nCores, type='FORK')
# #storeTimePar <- system.time( tsPredictions <- parSapply(cl, 1:dim(predictionGridGlobal)[1], generatePrediction, inputPredictionGrid=predictionGridGlobal ) )
# #storeTimePar <- system.time( tsPredictions <- parSapply(cl, 1:1000, generatePrediction, inputPredictionGrid=predictionGridGlobal ) )
# #storeTimePar <- system.time( tsPredictions <- pbsapply(1:2000, generatePrediction, inputPredictionGrid=predictionGridGlobal, cl = cl  ) )
# #storeTimePar <- system.time( tsPredictions <- sapply( 1:1000, generatePrediction, inputPredictionGrid=predictionGridGlobal  ) ) 
# stopCluster(cl)
# print( storeTimePar )
