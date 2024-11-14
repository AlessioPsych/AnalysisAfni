args <- commandArgs(T)
print( args )
mainDir <- getwd()
setwd( mainDir )

# to debug
rm( list=ls() ); gc();
setwd('/analyse/Project0226/KastnerModel/data_KastnerClassic/ASM16_KastnerVest')
args <- c('zzz_CCW_2components_05_params.nii.gz','zzz_CW_mean_blur_sm.nii.gz','ttt_nameOut','CW','05','CCW','05','05')
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
inputParams <- args[1]
inputTs <- args[2]
outName <- args[3]
startRotation <- args[4]
startCycles <- args[5]
targetRotation <- args[6]
targetCycles <- args[7]
targetComponent <- args[8]

print( sprintf('input parameters = %s', inputParams ) )
print( sprintf('input TS = %s', inputTs ) )
print( sprintf('outname = %s', outName ) )
print( sprintf('startRotation = %s', targetRotation ) )
print( sprintf('startCycles = %s', targetCycles ) )
print( sprintf('targetRotation = %s', targetRotation ) )
print( sprintf('targetCycles = %s', targetCycles ) )
print( sprintf('targetComponent = %s', targetComponent ) )

# get TR
instr <- sprintf('3dinfo -tr %s > __tttt.1D', inputTs); system( instr )
trValue <- read.table('__tttt.1D'); system('rm __tttt.1D')

# get number of dynamics
instr <- sprintf('3dinfo -nv %s > __tttt.1D', inputTs); system( instr )
nDynamics <- read.table('__tttt.1D'); system('rm __tttt.1D')

########################
# 5 cycles theta array #
########################
thetaCCW = linspace( 0, 2*pi, 11 );
thetaCW = linspace( 0, 2*pi, 11 );
thetaCW = thetaCW[ c( linspace( length(thetaCW), 1, length(thetaCW) ) ) ];

thetaExpCW <- ( as.numeric( repmat( thetaCW[1:length(thetaCW)-1], 1, 5) )  ) 
thetaExpCCW <- ( as.numeric( repmat( thetaCCW[1:length(thetaCCW)-1], 1, 5) ) )

thetaExpCW05 <- thetaExpCW 
thetaExpCCW05 <- thetaExpCCW 

thetaExpCW05 <- as.circular( thetaExpCW05, type='angles', units='radians', template='none', modulo='2pi', zero=0, rotation='counter' )
thetaExpCCW05 <- as.circular( thetaExpCCW05, type='angles', units='radians', template='none', modulo='2pi', zero=0, rotation='counter' )

# visualize 05 cycles
graphics.off()
x11( width=6, height=4 )
par( mfrow=c( 1,2 ) )
for ( counterImg in 1:length(thetaExpCW05) ) {
  plot( thetaExpCW05[counterImg], type='b', col='red', main='CW, 05 cycles' ); 
  plot( thetaExpCCW05[counterImg], type='b', col='red', main='CCW, 05 cycles' );
  par(new=FALSE); Sys.sleep( 0.1 )
}

#########################
# 10 cycles theta array #
#########################

thetaCCW = linspace( 0, 2*pi, 11 );
thetaCW = linspace( 0, 2*pi, 11 );
thetaCW = thetaCW[ c( linspace( length(thetaCW), 1, length(thetaCW) ) ) ];

thetaCW01 <- thetaCW[1:length(thetaCW)-1]
thetaCW01 <- c( thetaCW01[1:5], thetaCW01[6:10]+pi ) #- pi/10
thetaCCW01 <- thetaCCW[1:length(thetaCW)-1]
thetaCCW01 <- c( thetaCCW01[1:5], thetaCCW01[6:10]+pi ) #+ pi/10

thetaExpCW <- ( as.numeric( repmat( thetaCW01, 1, 5) ) ) 
thetaExpCCW <- ( as.numeric( repmat( thetaCCW01, 1, 5) ) ) 

thetaExpCW10 <- thetaExpCW
thetaExpCCW10 <- thetaExpCCW

thetaExpCW10 <- as.circular( thetaExpCW10, type='angles', units='radians', template='none', modulo='2pi', zero=0, rotation='counter' )
thetaExpCCW10 <- as.circular( thetaExpCCW10, type='angles', units='radians', template='none', modulo='2pi', zero=0, rotation='counter' )

# visualize 10 cycles
graphics.off()
x11( width=6, height=4 )
par( mfrow=c( 1,2 ) )
for ( counterImg in 1:length(thetaExpCW10) ) {
  plot( thetaExpCW10[counterImg], type='b', col='red', main='CW, 10 cycles' ); 
  plot( thetaExpCCW10[counterImg], type='b', col='red', main='CCW, 10 cycles' );
  par(new=FALSE); Sys.sleep( 0.1 )
}
graphics.off()

# MRI time and stimuli time
thetaExp <- thetaExpCW05 #just a variable to identify how many saccades have been performed, does not matter whether it is CW, CCW, 05 or 10
mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
timeStimuliOrig01 <- seq( 6.5, (max(mriTime)+as.numeric( trValue[1] )/2)-1.5, length.out = length(thetaExp) )
timeStimuliOrig <- c( timeStimuliOrig01 )
timeStimuli <- seq( 0, max(mriTime)+as.numeric( trValue[1] )/2, 0.1 )
thetaExpRes <- rep(-1,length(timeStimuli))
storeIdx <- rep(0,length(timeStimuliOrig))
testIdx <- rep(0,length(timeStimuliOrig))
for (i in 1:length(timeStimuliOrig)) {
  storeIdx[i] <- which( round( timeStimuliOrig[i], 2 )-timeStimuli < 0.00001 )[1]
  testIdx[i] <- timeStimuli[storeIdx[i]]
}
thetaExpResCW05 <- thetaExpRes
thetaExpResCCW05 <- thetaExpRes
thetaExpResCW10 <- thetaExpRes
thetaExpResCCW10 <- thetaExpRes
thetaExpResCW05[ storeIdx ] <- thetaExpCW05
thetaExpResCCW05[ storeIdx ] <- thetaExpCCW05
thetaExpResCW10[ storeIdx ] <- thetaExpCW10
thetaExpResCCW10[ storeIdx ] <- thetaExpCCW10
samplingTime <- as.numeric( nDynamics[1]*trValue[1] )/length(thetaExpRes)

# get model parameters from volume into 2D array
paramsFile <- read.AFNI( inputParams )
paramsVolume <- paramsFile$brk
idxVoxels <- which( abs( paramsVolume[,,,1] ) > 0 )
paramsArray <- array( 0, c( length( idxVoxels ), dim( paramsVolume )[4] ) )
for ( k in 1:dim( paramsArray)[2] ) {
  volTemp <- paramsVolume[,,,k]
  paramsArray[,k] <- volTemp[ idxVoxels ]    
}

# get ts from target data volume into 2D array
tsFile <- read.AFNI( inputTs )
tsVolume <- tsFile$brk
tsArray <- array( 0, c( length( idxVoxels ), dim( tsVolume )[4] ) )
for ( k in 1:dim( tsArray)[2] ) {
  volTemp <- tsVolume[,,,k]
  tsArray[,k] <- volTemp[ idxVoxels ]    
}

#### function to generate the predictions ####
generatePrediction <- function( inputPrfPar, samplingTime, mriTime, thetaExpResCW05, thetaExpResCCW05, thetaExpResCW10, thetaExpResCCW10, startRotation, targetRotation, startCycles, targetCycles, targetComponent ) {
  
  # the logic is the following: in the modeling exercise we have already accounted for hrf and orientation (CW or CCW)
  # now we take the parameters from half of the data, taken say from the CW part
  # and fit them to the CCW part. In order to do this we just need to take the fitted parameters (all)
  # and apply them to the target orientation (say CCW), to obtain what we would expect from a motor prf
  # presented with the target orientation (say CCW)
  
  plotFlag <- 0
  
  prfPar <- inputPrfPar # test by plotting
  #prfPar <- c(6, 12, 0.9, 0.9, 0.35, pi, 0.5, 1) # test by plotting
  
  if ( targetRotation=='CW' ) {
    if ( targetCycles=='05' ) { thetaExpRes_retest05 <- thetaExpResCW05 }
    if ( targetCycles=='10' ) { thetaExpRes_retest10 <- thetaExpResCW10 }
  }
  if ( targetRotation=='CCW' ) {
    if ( targetCycles=='05' ) { thetaExpRes_retest05 <- thetaExpResCCW05 }
    if ( targetCycles=='10' ) { thetaExpRes_retest10 <- thetaExpResCCW10 }
  }
  
  hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[1], a2=prfPar[2], b1=prfPar[3], b2=prfPar[4], c=prfPar[5]), verbose=FALSE ); #plot(hrf~seq(0,30,samplingTime))
  
  # motor prf 05
  angleOverTime05 <- suppressWarnings( dwrappednormal( thetaExpRes_retest05, mu=prfPar[6], rho=prfPar[7] ) ) 
  angleOverTime05[ thetaExpRes_retest05 == -1 ] <- 0 
  predictionLoop <- angleOverTime05
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('nearest') )
  predictionLoop01Interp05 <- scaleData( predictionLoop01Interp, 1, 0 )

  # motor prf 10
  angleOverTime10 <- suppressWarnings( dwrappednormal( thetaExpRes_retest10, mu=prfPar[6], rho=prfPar[7] ) ) 
  angleOverTime10[ thetaExpRes_retest10 == -1 ] <- 0 
  predictionLoop <- angleOverTime10
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('nearest') )
  predictionLoop01Interp10 <- scaleData( predictionLoop01Interp, 1, 0 )
  
  if (targetComponent=='both') {
    predictionLoop01InterpOut <- predictionLoop01Interp05*prfPar[8] + predictionLoop01Interp10*(1-prfPar[8])
  }
  if (targetComponent=='05') {
    predictionLoop01InterpOut <- predictionLoop01Interp05
  }
  if (targetComponent=='10') {
    predictionLoop01InterpOut <- predictionLoop01Interp10
  }

  returnPrediction <- round( scaleData( predictionLoop01InterpOut, 1, 0 ), 5 ); 
  
  if (plotFlag==1) {
    x11( width=8, height=5 )
    par(mfrow=c(3,3))
    plot( suppressWarnings( dwrappednormal( seq(0,2*pi,0.05), mu=prfPar[6], rho=prfPar[7] ) )~seq(0,2*pi,0.05), type='l', bty='n', las='1' )
    plot( thetaExpRes_retest05 ~ timeStimuli, bty='n', las=1, type='l'  )
    plot( angleOverTime05 ~ timeStimuli, bty='n', las=1, type='l' )
    plot( angleOverTime10 ~ timeStimuli, bty='n', las=1, type='l' )
    plot( predictionLoop01Interp05 ~ mriTime, bty='n', las=1, type='l' )
    plot( predictionLoop01Interp10 ~ mriTime, bty='n', las=1, type='l' )
    plot( returnPrediction ~ mriTime, bty='n', las=1, type='l'  )
    plot( hrf~seq(0,30,samplingTime) )
  }
  
  return( returnPrediction )
  
}


glmFun <- function( ts_observed, ts_predicted ) { #from here, complete the glmfun to compare predicted and observed
  dMat01 <- cbind( rep(1,length(dMat)), dMat ) #column of ones and column of predictor
  a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope)
  expectedTs <- crossprod( t(dMat01), a ) #expected ts
  residualsSquared <- ( selTsVoxel - expectedTs )^2 
  ssRes <- apply( residualsSquared, 2, sum )
  r2 <- 1-(ssRes/ssTot) #r squares
  
}

# fit the predicted data to the observed
system.time(
for ( k in 1:1000 ) { dim( paramsArray)[1]
  motPrfIn <- paramsArray[k,]
  ts_observed <- tsArray[k,]
  ts_predicted <- generatePrediction( motPrfIn, samplingTime, mriTime, thetaExpResCW05, thetaExpResCCW05, thetaExpResCW10, thetaExpResCCW10, startRotation, targetRotation, startCycles, targetCycles, targetComponent )
  modTemp <- lm( ts_observed~ts_predicted )
}
)

a1 <- generatePrediction( paramsArray[150,seq(1,8)], samplingTime, mriTime, thetaExpResCW05, thetaExpResCCW05, thetaExpResCW10, thetaExpResCCW10, startRotation, targetRotation, startCycles, targetCycles, targetComponent )
plot( a1, type='l' )






















predictionWrapper <- function( passIdx, limitsMat, inputGrid ) {
  
  print( sprintf('iteration: %1.0f of %1.0f', passIdx, dim( limitsMat )[1] ) )
  selectedRows <- limitsMat[passIdx, ]
  
  # generate predictions:
  nCores <- 10
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( tsPredictions <- parSapply(cl, selectedRows[1]:selectedRows[2], generatePrediction,
                                                          inputPredictionGrid=inputGrid,
                                                          samplingTime=samplingTime,
                                                          thetaExpRes=thetaExpRes,
                                                          scaleData=scaleData,
                                                          timeStimuli=timeStimuli,
                                                          mriTime=mriTime,
                                                          thetaExpRes10=thetaExpRes10,
                                                          components=components,
                                                          thetaExpRes_side=thetaExpRes_side ) )
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

