args <- commandArgs(T)
print( args )
mainDir <- getwd()
setwd( mainDir )

# to debug
#rm( list=ls() ); gc();
#setwd('/analyse/Project0226/KastnerModel/data_KastnerModel/AFH28_201119')
#args <- c('KASTMDL_stc_mc_tu_dt_avg.nii.gz','zzz_cw_del/','0','baselineModel')
#mainDir <- getwd()
#setwd( mainDir )
# limit width only for those that are more sensitive for the fitting results (telling apart left and right)
# try the variant with different angles also
# keep orig version and shift fitted phases

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/vonMisesDist.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( neuRosim )
library( parallel )
library( circular )

# arrange inputs
inputFile <- args[1]
outDirectory <- args[2]
fineFit <- as.numeric( args[3] )
modelType <- args[4]

print( sprintf('input file = %s', inputFile ) )
print( sprintf('output directory = %s', outDirectory ) )
print( sprintf('fineFit = %1.0f', fineFit ) )
print( sprintf('modelType = %s', modelType ) )

# create output directory
system( sprintf('mkdir %s', outDirectory) )

# get TR
instr <- sprintf('3dinfo -tr %s > __tttt.1D', inputFile); system( instr )
trValue <- read.table('__tttt.1D'); system('rm __tttt.1D')

# get number of dynamics
instr <- sprintf('3dinfo -nv %s > __tttt.1D', inputFile); system( instr )
nDynamics <- read.table('__tttt.1D'); system('rm __tttt.1D')

# orientation <- 'cw'
# if ( orientation=='ccw' ) {
# theta = linspace( 0-pi/2, 2*pi-pi/2, 11 );
# }
# if ( orientation=='cw' ) {
# theta = linspace( 0-pi/2, 2*pi-pi/2, 11 );
# theta = theta[ c( linspace( length(theta), 1, length(theta) ) ) ];
# }

theta01 <- linspace( 0-pi/2, 2*pi-pi/2-0.571, 12 );
theta02 <- linspace( 0-pi/2+0.571, 2*pi-pi/2, 12 );
theta02 <- theta02[ c( linspace( length(theta02), 1, length(theta02) ) ) ];

#theta <- c(0, theta01, 0, 0, theta02, 0, 0, 0, theta02+pi, 0, 0, theta01+pi, 0, 0 ) 

wrapFunction <- function( idxIn, argIn, limitIn ) {
  if( argIn[idxIn] < limitIn ) 
    { return( argIn[idxIn] ) 
    }else {
    return( argIn[idxIn]-limitIn )
  }
}
#sapply( 1:length(seq(0,3*pi,0.2)), wrapFunction, argIn=seq(0,3*pi,0.2), limitIn=2*pi )

theta <- c(0, theta02+pi, 0, 0, theta01+pi, 0, 0, 0, theta01, 0, 0, theta02, 0, 0 ) 

cartCoords <- pol2cart( cbind( theta, rep(1,length(theta)) ) )

polCoords <- cart2pol( cartCoords )

#theta <- c(0, theta02, 0, 0, theta01, 0, 0, 0, theta01, 0, 0, theta02, 0, 0 )

thetaExp <- round( ( polCoords[,1] + pi*3/2 ), 2 )
thetaExp <- sapply( 1:length(thetaExp), wrapFunction, argIn=thetaExp, limitIn=2*pi )
thetaExp[ theta==0 ] <- 0

x11( width=4, height=4 )
for ( counterImg in 1:length(thetaExp) ) {
  if (abs( thetaExp[counterImg] ) < 0.001 ) {
    polar( thetaExp[counterImg], 1, type='n', col='black'  ); par(new=FALSE); Sys.sleep( 0.2 )
  }else{
    polar( thetaExp[counterImg], 1, type='o', col='red', cex=2, pch=16 ); par(new=FALSE); Sys.sleep( 0.2 )
  }
}

# for ( counterImg in 1:length(thetaExp) ) {
#   if (abs( thetaExp[counterImg] ) < 0.001 ) {
#     plot( cartCoords[counterImg,1], cartCoords[counterImg,3], type='n', col='black', xlim=c(-5,5), ylim=c(-5,5)  ); par(new=FALSE); Sys.sleep( 0.2 )
#   }else{
#     plot( cartCoords[counterImg,2], cartCoords[counterImg,4], type='o', col='red', cex=2, pch=16, xlim=c(-5,5), ylim=c(-5,5) ); par(new=FALSE); Sys.sleep( 0.2 )
#   }
# }


#this part of the code builds a matrix with all the possible prediction tested
if (fineFit==1) {
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  hrfb1Elements <- 1
  hrfb2Elements <- 1
  hrfcElements <- 1
  hrf_a1 <- seq( 6, 7, length.out=hrfDelayOnsetElements )
  hrf_a2 <- seq( 11, 16, length.out=hrfDelayUnderShootElements )
  hrf_b1 <- seq( 0.9, 1, length.out=hrfb1Elements  ) 
  hrf_b2 <- seq( 0.9, 1, length.out=hrfb2Elements  ) 
  hrf_c <- seq( 0.35, 0.5, length.out=hrfcElements  ) 
  muPar <- c( seq( 0.02, 6.26, length.out = 22 ) ) #58
  #muPar <- c( seq( 6.28-(6.28*0.99), 6.28*0.99, length.out = 45 ) ) #58
  kappaPar <- seq( 0.05, 0.95, length.out = 18 )
  nonLin <- 1
}
if (fineFit==0) {
  hrfDelayOnsetElements <- 4
  hrfDelayUnderShootElements <- 3
  hrfb1Elements <- 1
  hrfb2Elements <- 1
  hrfcElements <- 1
  hrf_a1 <- seq( 5, 10, length.out=hrfDelayOnsetElements )
  hrf_a2 <- seq( 11, 16, length.out=hrfDelayUnderShootElements )
  hrf_b1 <- seq( 0.9, 1.2, length.out=hrfb1Elements  ) 
  hrf_b2 <- seq( 0.9, 1.2, length.out=hrfb2Elements  ) 
  hrf_c <- seq( 0.35, 0.6, length.out=hrfcElements  ) 
  muPar <- c( seq( 0+0.02, 2*pi-0.02, length.out = 22 ) ) #58
  #muPar <- c( seq( 6.28-(6.28*0.95), 6.28*0.95, length.out = 8 ) ) #58
  kappaPar <- seq( 0.05, 0.95, length.out = 24 )
  nonLin <- seq(1)
}

if (modelType=='baselineModel') { 
  muPar <- c( seq( 0+0.02, 2*pi-0.02, length.out = 1 ) )
  kappaPar <- seq( 0.05, 0.95, length.out = 1 )
}

predictionGridGlobal <- expand.grid( hrf_a1, hrf_a2, hrf_b1, hrf_b2, hrf_c, muPar, kappaPar, nonLin)
dim( predictionGridGlobal )

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
# if ( orientation=='whole' ) {
#   mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
#   timeStimuliOrig01 <- seq( 6.5, (max(mriTime)+as.numeric( trValue[1] )/2)/2-1.5, length.out = length(thetaExp)/2 )
#   timeStimuliOrig02 <- timeStimuliOrig01 + (max(mriTime)+as.numeric( trValue[1] )/2)/2
#   timeStimuliOrig <- c( timeStimuliOrig01, timeStimuliOrig02 )
#   timeStimuli <- seq( 0, max(timeStimuliOrig)+1.5, 0.1 )
#   thetaExpRes <- rep(0,length(timeStimuli))
#   storeIdx <- rep(0,length(timeStimuliOrig))
#   testIdx <- rep(0,length(timeStimuliOrig))
#   for (i in 1:length(timeStimuliOrig)) {
#     storeIdx[i] <- which( round( timeStimuliOrig[i], 2 )-timeStimuli < 0.00001 )[1]
#     testIdx[i] <- timeStimuli[storeIdx[i]]
#   }
#   thetaExpRes[ storeIdx ] <- thetaExp
#   samplingTime <- as.numeric( nDynamics[1]*trValue[1] )/length(thetaExpRes)
# } else {
  mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
  timeStimuliOrig01 <- seq( 0, (max(mriTime)+as.numeric( trValue[1] )/2), length.out = length(thetaExp) )
  timeStimuliOrig <- c( timeStimuliOrig01 )
  timeStimuli <- seq( 0, max(mriTime)+as.numeric( trValue[1] )/2, 0.1 )
  thetaExpRes <- rep(0,length(timeStimuli))
  storeIdx <- rep(0,length(timeStimuliOrig))
  testIdx <- rep(0,length(timeStimuliOrig))
  for (i in 1:length(timeStimuliOrig)) {
    storeIdx[i] <- which( round( timeStimuliOrig[i], 2 )-timeStimuli < 0.00001 )[1]
    testIdx[i] <- timeStimuli[storeIdx[i]]
  }
  thetaExpRes[ storeIdx ] <- thetaExp
  samplingTime <- as.numeric( nDynamics[1]*trValue[1] )/length(thetaExpRes)
#}
#plot( thetaExpRes ~ timeStimuli, type='h'  )
#plot( ifelse(thetaExpRes>0,1,0) ~ timeStimuli, type='h'  )

#### function to generate the predictions ####
generatePrediction <- function( indexPrediction, inputPredictionGrid, samplingTime, thetaExpRes, scaleData, timeStimuli, mriTime, modelType ) {
  
  library( neuRosim )
  library( pracma )
  library( circular )
  prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
  
  #tapply( thetaExp, list( directionsChange ), mean )
  #modelType <- 'kinematic_baselineModel'
  #prfPar <- c( 5, 10, 0.9, 0.9, 0.35, 4.5, 0.8, 1  ); graphics.off()
  
  hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[1], a2=prfPar[2], b1=prfPar[3], b2=prfPar[4], c=prfPar[5]), verbose=FALSE )
  
  #tuningFunction <- suppressWarnings( dwrappednormal( seq(0,2*pi,0.05), mu=prfPar[6], rho=prfPar[7] ) )
  #angleOverTimeAll01 <- suppressWarnings( dwrappednormal( thetaExp, mu=prfPar[6], rho=prfPar[7] ) )
  #angleOverTimeAll01[ baseline ] <- 0
  #angleOverTimeAll <- suppressWarnings( dwrappednormal( thetaExp, mu=prfPar[6], rho=prfPar[7] ) )
  #angleOverTimeAll <- suppressWarnings( dvonmises( thetaExp, mu=prfPar[6], kappa=prfPar[7] ) )
  #angleOverTimeAll <- suppressWarnings( vonMisesDist( thetaExp, mu=prfPar[6], k=prfPar[7], scaleOutput = 1 ) )
  #spFit <- smooth.spline( seq(-pi,pi,0.05), vonMisesDist( seq(-pi,pi,0.05), mu = prfPar[6], k = prfPar[7], scaleOutput = 0 ) )
  #angleOverTimeAll <- predict( spFit, thetaExpRes )$y
  angleOverTimeAll <- suppressWarnings( dwrappednormal( thetaExpRes, mu=prfPar[6], rho=prfPar[7] ) ) #plot( angleOverTimeAll, type='l' )
  angleOverTimeAll[ thetaExpRes==0 ] <- 0 #plot( angleOverTimeAll, type='l' )
  #angleOverTimeAll <- scaleData( angleOverTimeAll, 1, 0 )
  onOffPredictor <- ifelse(thetaExpRes>0,1,0)
  
  if (modelType=='kinematicModel') {
    predictionLoop <- angleOverTimeAll; 
  }
  if (modelType=='baselineModel') {
    predictionLoop <- onOffPredictor; 
  }
  if (modelType=='kinematic_baselineModel') {
    predictionLoop <- scaleData( onOffPredictor, max( angleOverTimeAll ), min( angleOverTimeAll ) ) + angleOverTimeAll; 
  }
  
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]; #plot(pConvTrim01,type='l')

  #predictionLoop01Interp <- interp1( x=timeStimuli, y=pConvTrim01, xi=mriTime, method=c('nearest') )
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('nearest') )
  
  returnPrediction <- round( scaleData( predictionLoop01Interp, 1, 0 ), 5 ); #plot(returnPrediction,type='l')
  
  # x11( width=6, height=6 )
  #par(mfrow=c(3,2))
  #plot( thetaExpRes~timeStimuli, bty='n', las=1, type='l'   )
  #plot( suppressWarnings( dwrappednormal( seq(0,2*pi,0.05), mu=prfPar[6], rho=prfPar[7] ) )~seq(0,2*pi,0.05), type='l', bty='n', las='1' )
  #plot( onOffPredictor ~ timeStimuli, bty='n', las=1, type='l' )
  #plot( angleOverTimeAll ~ timeStimuli, bty='n', las=1, type='l' )
  #plot( returnPrediction ~ mriTime, bty='n', las=1, type='l'  )
  #plot( predictionLoop ~ timeStimuli, bty='n', las=1, type='l'  )
  #plot( hrf~seq(0,30,samplingTime) )
  
  return( returnPrediction )
  
}

predictionWrapper <- function( passIdx, limitsMat, inputGrid ) {
  
  print( sprintf('iteration: %1.0f of %1.0f', passIdx, dim( limitsMat )[1] ) )
  selectedRows <- limitsMat[passIdx, ]
  
  # generate predictions:
  nCores <- 8
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( tsPredictions <- parSapply(cl, selectedRows[1]:selectedRows[2], generatePrediction,
                                                          inputPredictionGrid=inputGrid,
                                                          samplingTime=samplingTime,
                                                          thetaExpRes=thetaExpRes,
                                                          scaleData=scaleData,
                                                          timeStimuli=timeStimuli,
                                                          mriTime=mriTime,
                                                          modelType=modelType ) )
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

