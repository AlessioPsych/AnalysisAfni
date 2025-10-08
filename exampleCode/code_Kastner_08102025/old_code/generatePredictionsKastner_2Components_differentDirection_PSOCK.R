args <- commandArgs(T)
print( args )
mainDir <- getwd()
setwd( mainDir )

# to debug
#rm( list=ls() ); gc();
#setwd('/analyse/Project0226/KastnerModel/data_KastnerClassic/ASM16_KastnerVest')
#args <- c('zzz_CW_mean_blur_sm.nii.gz','zzz_cw_del/','0','cw','05')
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
orientation <- args[4]
components <- args[5]

print( sprintf('input file = %s', inputFile ) )
print( sprintf('output directory = %s', outDirectory ) )
print( sprintf('fineFit = %1.0f', fineFit ) )

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

########################
# 5 cycles theta array #
########################
#thetaCCW = linspace( 0-pi/2, 2*pi-pi/2, 11 );
#thetaCW = linspace( 0-pi/2, 2*pi-pi/2, 11 );
#thetaCW = thetaCW[ c( linspace( length(thetaCW), 1, length(thetaCW) ) ) ];

#thetaExpCW <- ( as.numeric( repmat( thetaCW[1:length(thetaCW)-1], 1, 5) ) + 3.14 ) %% (2*pi)
#thetaExpCCW <- ( as.numeric( repmat( thetaCCW[1:length(thetaCCW)-1], 1, 5) ) + 3.14 ) %% (2*pi)
#thetaExp <- c( thetaExpCW, thetaExpCCW )

thetaCCW = linspace( 0, 2*pi, 11 );
thetaCW = linspace( 0, 2*pi, 11 );
thetaCW = thetaCW[ c( linspace( length(thetaCW), 1, length(thetaCW) ) ) ];

thetaExpCW <- ( as.numeric( repmat( thetaCW[1:length(thetaCW)-1], 1, 5) )  ) 
thetaExpCCW <- ( as.numeric( repmat( thetaCCW[1:length(thetaCCW)-1], 1, 5) ) )
thetaExp <- c( thetaExpCW, thetaExpCCW )

if ( orientation=='cw' ) {
  thetaExp <- thetaExpCW 
}
if ( orientation=='ccw' ) {
  thetaExp <- thetaExpCCW
}

thetaExp <- circular( thetaExp, type='angles', units='radians', 
          template='none', modulo='asis', zero=0, rotation='counter' )

graphics.off()
x11( width=7, height=4 )
par(mfrow=c(1,2))
for ( counterImg in 1:length(thetaExp) ) {
  polar( c( 0, thetaExp[counterImg] ), c(0, 1) , 'o', col='red', lwd=2  ); par(new=FALSE); Sys.sleep( 0.3 )
  plot( thetaExp[counterImg] ); Sys.sleep( 0.3 )
}

#test <- circular( thetaExp[1:3], type='angles', units='radians', 
#                     template='none', modulo='asis', zero=0, rotation='counter' )
#muPar <- circular( 1.5, type='angles', units='radians', 
#                      template='none', modulo='asis', zero=0, rotation='counter' )
#rhoPar <- 0.25
#plot( test )
#dwrappednormal( test, muPar=1.5, rho=rhoPar )

#########################
# 10 cycles theta array #
#########################
# thetaCCW = linspace( 0-pi/2, 2*pi-pi/2, 11 );
# thetaCW = linspace( 0-pi/2, 2*pi-pi/2, 11 );
# thetaCW = thetaCW[ c( linspace( length(thetaCW), 1, length(thetaCW) ) ) ];
# 
# thetaCW01 <- thetaCW[1:length(thetaCW)-1]
# thetaCW01 <- c( thetaCW01[1:5], thetaCW01[6:10]+pi ) - pi/10
# thetaCCW01 <- thetaCCW[1:length(thetaCW)-1]
# thetaCCW01 <- c( thetaCCW01[1:5], thetaCCW01[6:10]+pi ) + pi/10
# 
# theta_sideParameter <- rep(0,length(thetaCCW01))
# theta_sideParameter[1:5] <- 1
# theta_sideParameter[6:10] <- 2
# 
# thetaExpCW <- ( as.numeric( repmat( thetaCW01, 1, 5) ) + 3.14 ) %% (2*pi)
# thetaExpCCW <- ( as.numeric( repmat( thetaCCW01, 1, 5) ) + 3.14 ) %% (2*pi)
# thetaExp_sideParameter <- as.numeric( repmat( theta_sideParameter, 1, 5) )

thetaCCW = linspace( 0, 2*pi, 11 );
thetaCW = linspace( 0, 2*pi, 11 );
thetaCW = thetaCW[ c( linspace( length(thetaCW), 1, length(thetaCW) ) ) ];

thetaCW01 <- thetaCW[1:length(thetaCW)-1]
thetaCW01 <- c( thetaCW01[1:5], thetaCW01[6:10]+pi ) #- pi/10
thetaCCW01 <- thetaCCW[1:length(thetaCW)-1]
thetaCCW01 <- c( thetaCCW01[1:5], thetaCCW01[6:10]+pi ) #+ pi/10

theta_sideParameter <- rep(0,length(thetaCCW01))
theta_sideParameter[1:5] <- 1
theta_sideParameter[6:10] <- 2

thetaExpCW <- ( as.numeric( repmat( thetaCW01, 1, 5) ) ) 
thetaExpCCW <- ( as.numeric( repmat( thetaCCW01, 1, 5) ) ) 
thetaExp_sideParameter <- as.numeric( repmat( theta_sideParameter, 1, 5) )

if ( orientation=='cw' ) {
  thetaExp10 <- thetaExpCW
}
if ( orientation=='ccw' ) {
  thetaExp10 <- thetaExpCCW
}

thetaExp10 <- circular( thetaExp10, type='angles', units='radians', 
                      template='none', modulo='asis', zero=0, rotation='counter' )

graphics.off()
x11( width=7, height=4 )
par(mfrow=c(1,2))
for ( counterImg in 1:length(thetaExp) ) {
  polar( c( 0, thetaExp10[counterImg] ), c(0, 1) , 'o', col='red', lwd=2  ); par(new=FALSE); Sys.sleep( 0.3 )
  plot( thetaExp10[counterImg] ); Sys.sleep( 0.3 ); Sys.sleep( 0.3 )
}


#this part of the code builds a matrix with all the possible prediction tested
if (fineFit==1) {
  hrfDelayOnsetElements <- 2
  hrfDelayUnderShootElements <- 2
  hrfb1Elements <- 3
  hrfb2Elements <- 3
  hrfcElements <- 1
  hrf_a1 <- seq( 6, 7, length.out=hrfDelayOnsetElements )
  hrf_a2 <- seq( 11, 16, length.out=hrfDelayUnderShootElements )
  hrf_b1 <- seq( 0.9, 1.2, length.out=hrfb1Elements  ) 
  hrf_b2 <- seq( 0.9, 1.2, length.out=hrfb2Elements  ) 
  hrf_c <- seq( 0.35, 0.5, length.out=hrfcElements  ) 
  muPar <- c( seq( 0.02, 6.26, length.out = 20 ) ) #58
  #muPar <- c( seq( 6.28-(6.28*0.99), 6.28*0.99, length.out = 45 ) ) #58
  kappaPar <- seq( 0.05, 0.95, length.out = 18 )
  if (components=='both') {
    multPar <- seq( 0.1, 0.9, length.out = 5 )
  } else{
    multPar <- 1
  }
  
  
}
if (fineFit==0) {
  hrfDelayOnsetElements <- 1
  hrfDelayUnderShootElements <- 1
  hrfb1Elements <- 3 #3
  hrfb2Elements <- 3 #3
  hrfcElements <- 1
  hrf_a1 <- seq( 6, 10, length.out=hrfDelayOnsetElements )
  hrf_a2 <- seq( 12, 16, length.out=hrfDelayUnderShootElements )
  hrf_b1 <- seq( 0.9, 1.4, length.out=hrfb1Elements  ) 
  hrf_b2 <- seq( 0.9, 1.4, length.out=hrfb2Elements  ) 
  hrf_c <- seq( 0.35, 0.6, length.out=hrfcElements  ) 
  muPar <- c( seq( 0+0.04, 2*pi-0.04, length.out = 22 ) ) #58
  #muPar <- c( seq( 6.28-(6.28*0.95), 6.28*0.95, length.out = 8 ) ) #58
  kappaPar <- seq( 0.05, 0.95, length.out = 18 )
  if (components=='both') {
    multPar <- seq( 0.05, 0.95, length.out = 5 )
  } else{
    multPar <- 1
  }
}

predictionGridGlobal <- expand.grid( hrf_a1, hrf_a2, hrf_b1, hrf_b2, hrf_c, muPar, kappaPar, multPar)
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

# MRI time and stimuli time, 5 cycles
mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
timeStimuliOrig01 <- seq( 6.5, (max(mriTime)+as.numeric( trValue[1] )/2)-1.5, length.out = length(thetaExp) )
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
thetaExpRes <- circular( thetaExpRes, type='angles', units='radians', 
                        template='none', modulo='asis', zero=0, rotation='counter' )


# MRI time and stimuli time, 10 cycles
mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
timeStimuliOrig01 <- seq( 6.5, (max(mriTime)+as.numeric( trValue[1] )/2)-1.5, length.out = length(thetaExp10) )
timeStimuliOrig <- c( timeStimuliOrig01 )
timeStimuli <- seq( 0, max(mriTime)+as.numeric( trValue[1] )/2, 0.1 )
thetaExpRes10 <- rep(0,length(timeStimuli))
storeIdx <- rep(0,length(timeStimuliOrig))
testIdx <- rep(0,length(timeStimuliOrig))
for (i in 1:length(timeStimuliOrig)) {
  storeIdx[i] <- which( round( timeStimuliOrig[i], 2 )-timeStimuli < 0.00001 )[1]
  testIdx[i] <- timeStimuli[storeIdx[i]]
}
thetaExpRes10[ storeIdx ] <- thetaExp10
thetaExpRes10 <- circular( thetaExpRes10, type='angles', units='radians', 
                         template='none', modulo='asis', zero=0, rotation='counter' )

# MRI time and stimuli time, side_parameter
mriTime <- seq( 0, as.numeric( nDynamics[1]*trValue[1] )-as.numeric( trValue[1] ), as.numeric( trValue[1] ) ) + as.numeric( trValue[1] )/2
timeStimuliOrig01 <- seq( 6.5, (max(mriTime)+as.numeric( trValue[1] )/2)-1.5, length.out = length(thetaExp_sideParameter) )
timeStimuliOrig <- c( timeStimuliOrig01 )
timeStimuli <- seq( 0, max(mriTime)+as.numeric( trValue[1] )/2, 0.1 )
thetaExpRes_side <- rep(0,length(timeStimuli))
storeIdx <- rep(0,length(timeStimuliOrig))
testIdx <- rep(0,length(timeStimuliOrig))
for (i in 1:length(timeStimuliOrig)) {
  storeIdx[i] <- which( round( timeStimuliOrig[i], 2 )-timeStimuli < 0.00001 )[1]
  testIdx[i] <- timeStimuli[storeIdx[i]]
}
thetaExpRes_side[ storeIdx ] <- thetaExp_sideParameter


#### function to generate the predictions ####
generatePrediction <- function( indexPrediction, inputPredictionGrid, samplingTime, thetaExpRes, scaleData, timeStimuli, mriTime, thetaExpRes10, components, thetaExpRes_side ) {
  
  library( neuRosim )
  library( pracma )
  library( circular )
  prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
  
  #prfPar <- c( 6, 12, 1.4, 1.4, 0.35, 1/2*pi, 0.8, 0.1  ); components <- '5'; #graphics.off()
  
  hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[1], a2=prfPar[2], b1=prfPar[3], b2=prfPar[4], c=prfPar[5]), verbose=FALSE ); #plot(hrf~seq(0,30,samplingTime))
  
  muPar <- circular( prfPar[6], type='angles', units='radians', 
                                    template='none', modulo='asis', zero=0, rotation='counter' )
  
  # 5 cycles parameter, motor prf
  angleOverTimeAll05 <- dwrappednormal( thetaExpRes, mu=muPar, rho=prfPar[7] ) #plot( angleOverTimeAll, type='l' )
  angleOverTimeAll05[ thetaExpRes==0 ] <- 0 #plot( angleOverTimeAll, type='l' )
  predictionLoop <- angleOverTimeAll05; 
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]; #plot(pConvTrim01,type='l')
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('nearest') )
  predictionLoop01Interp <- scaleData( predictionLoop01Interp, 1, 0 )
  
  # 10 cycles parameter, motor prf
  angleOverTimeAll10 <- suppressWarnings( dwrappednormal( thetaExpRes10, mu=muPar, rho=prfPar[7] ) ) #plot( angleOverTimeAll, type='l' )
  angleOverTimeAll10[ thetaExpRes==0 ] <- 0 #plot( angleOverTimeAll, type='l' )
  predictionLoop <- angleOverTimeAll10; 
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]; #plot(pConvTrim01,type='l')
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp10 <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('nearest') )
  predictionLoop01Interp10 <- scaleData( predictionLoop01Interp10, 1, 0 )

  # side parameters
  angleOverTimeAll_side01 <- as.numeric( thetaExpRes_side == 1 )
  angleOverTimeAll_side01[ thetaExpRes==0 ] <- 0 #plot( angleOverTimeAll, type='l' )
  predictionLoop <- angleOverTimeAll_side01; 
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]; #plot(pConvTrim01,type='l')
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp_side01 <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('nearest') )
  predictionLoop01Interp_side01 <- scaleData( predictionLoop01Interp_side01, 1, 0 )
  
  angleOverTimeAll_side02 <- as.numeric( thetaExpRes_side == 2 )
  angleOverTimeAll_side02[ thetaExpRes==0 ] <- 0 #plot( angleOverTimeAll, type='l' )
  predictionLoop <- angleOverTimeAll_side02; 
  pConv01 <- conv( predictionLoop, hrf )
  pConvTrim01 <- pConv01[ 1 : length( timeStimuli ) ]; #plot(pConvTrim01,type='l')
  pFit <- smooth.spline( timeStimuli, pConvTrim01, spar=0.5 )
  predictionLoop01Interp_side02 <- interp1( x=pFit$x, y=pFit$y, xi=mriTime, method=c('linear') )
  predictionLoop01Interp_side02 <- scaleData( predictionLoop01Interp_side02, 1, 0 )
  
  if (components=='both') {
    predictionLoop01InterpOut <- predictionLoop01Interp*prfPar[8] + predictionLoop01Interp10*(1-prfPar[8])
  }
  if (components=='5') {
    predictionLoop01InterpOut <- predictionLoop01Interp
  }
  if (components=='10') {
    predictionLoop01InterpOut <- predictionLoop01Interp10
  }
  if (components=='10side') {
    predictionLoop01InterpOut <- predictionLoop01Interp_side01*prfPar[8] + predictionLoop01Interp_side02*(1-prfPar[8])
  }
  
  returnPrediction <- round( scaleData( predictionLoop01InterpOut, 1, 0 ), 5 ); #plot(returnPrediction,type='l')
  
  #x11( width=8, height=5 )
  #par(mfrow=c(3,3))
  #plot( suppressWarnings( dwrappednormal( seq(0,2*pi,0.05), mu=prfPar[6], rho=prfPar[7] ) )~seq(0,2*pi,0.05), type='l', bty='n', las='1' )
  #plot( angleOverTimeAll05 ~ timeStimuli, bty='n', las=1, type='l' )
  #plot( angleOverTimeAll10 ~ timeStimuli, bty='n', las=1, type='l' )
  #plot( angleOverTimeAll_side01 ~ timeStimuli, bty='n', las=1, type='l' )
  #plot( angleOverTimeAll_side02 ~ timeStimuli, bty='n', las=1, type='l' )
  #plot( returnPrediction ~ mriTime, bty='n', las=1, type='l'  )
  #plot( hrf~seq(0,30,samplingTime) )
  
  return( returnPrediction )
  
}

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

