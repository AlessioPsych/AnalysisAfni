# to debug
debugFlag <- 0
if ( debugFlag==1) {
  mainDir <- '/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner/code/models'
  setwd(mainDir)
  args <- c('/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner/expConditions/kastnerFinerCw/ASM16_kastnerFinerCw.nii.gz',
            '/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner/modelsOutput/pRF/ASM16/ASM16_pRF_params.nii.gz',
            'cw_completeStimuliKastnerOrig_borderOnlyAfterSaccade_shifted_final.Rdata',
            'cw',
            '/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner/modelsOutput/kastnerClassic_contrastBased/ASM16',
            'ASM16',
            '/media/alessiofracasso/DATADRIVE1/KastnerData/staging_area_Kastner/modelsOutput/kastnerClassic/ASM16/ASM16_model_wAverage_05.nii.gz')
}
if ( debugFlag==0 ) {
  args <- commandArgs(T)
  print( args )
}

print( 'current folder...' )
print( getwd( ) )

inputEyeEPI <- args[1] 
sprintf('input eye mov EPI file:')
sprintf('%s',inputEyeEPI)

inputpRFModel <- args[2]
sprintf('input pRF model file:')
sprintf('%s',inputpRFModel)

inputStimuliFile <- args[3]
sprintf('input stimuli file:')
sprintf('%s', inputStimuliFile )

inputDirection <- args[4]
sprintf('input direction:')
sprintf('%s', inputDirection )

outputFolder <- args[5]
sprintf('output folder:')
sprintf('%s', outputFolder )

participantName <- args[6]
sprintf('participant name:')
sprintf('%s', participantName )

dspRFName <- args[7]
sprintf('input dspRF name:')
sprintf('%s', dspRFName )

generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( neuRosim )
library( parallel )
library( circular )

# rescale function
rescale_var <- function(x, new_min = 0, new_max = 1) {
  old_min <- min(x, na.rm = TRUE)
  old_max <- max(x, na.rm = TRUE)
  scaled <- (x - old_min) / (old_max - old_min)
  scaled * (new_max - new_min) + new_min
}

setwd( outputFolder )
print( getwd() ) 

# get TR
print('get TR...')
instr <- sprintf('3dinfo -tr %s > __tttt.1D', inputEyeEPI); system( instr )
trValue <- read.table('__tttt.1D'); system('rm __tttt.1D')

# get number of dynamics
print('get number of dynamics...')
instr <- sprintf('3dinfo -nv %s > __tttt.1D', inputEyeEPI); system( instr )
nDynamics <- read.table('__tttt.1D'); system('rm __tttt.1D')

# load stimuli
print('load stimuli...')
load( inputStimuliFile ) 
print('stimuli dimension...')
dim( completeStimuliSequence )

# load Kastner classic data
print('load KastnerClassic data...')
tsFile <- read.AFNI( inputEyeEPI )
tsArray <- array( tsFile$brk, c( prod( dim( tsFile$brk )[1:3] ), dim(tsFile$brk)[4] ) )
print('ts dimension...')
print( dim( tsArray ) )

# load pRF data
print('load pRF data...')
pRFFile <- read.AFNI( inputpRFModel )
pRFArray <- array( pRFFile$brk, c( prod( dim( pRFFile$brk )[1:3] ), dim(pRFFile$brk)[4] ) )
print('pRF dimension...')
print( dim( pRFArray ) )

# load dspRF data
print('load dspRF data...')
dspRFFile <- read.AFNI( dspRFName )
dspRFArray <- array( dspRFFile$brk, c( prod( dim( dspRFFile$brk )[1:3] ), dim(dspRFFile$brk)[4] ) )
print('dspRF dimension...')
print( dim( dspRFArray ) )

varExp <- pRFArray[,12] # pRF varExp
varDspRF <- dspRFArray[,9] #dspRF varExp
selectedIdxAll <- which( varExp > 0.03 & varDspRF > 0.03 ) 
#selectedIdx <- selectedIdxAll

#split between chuncks for parallel fitting (about 1000 voxels each)
nVoxels <- length( selectedIdxAll ) 
totalVoxelsIterations <- ceil( nVoxels / 1000 )
limitsSelVoxels <- round( seq(1, nVoxels, length.out=totalVoxelsIterations) ) 
limitsSelVoxels[1] <- 1
limitsSelVoxels[ length(limitsSelVoxels) ] <- nVoxels
limitsSelVoxelsMatrix <- array( 0, c( (length(limitsSelVoxels)-1) , 2 ) )
limitsSelVoxelsMatrix[,1] <- limitsSelVoxels[1:(length(limitsSelVoxels)-1)]
limitsSelVoxelsMatrix[,2] <- limitsSelVoxels[2:(length(limitsSelVoxels))]
limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] <- limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk
runIndex <- seq( 1:dim(limitsSelVoxelsMatrix)[1] )

for ( nIterations in runIndex ) { # nIterations <- 1
  
  selectedVoxelsPassTemp <- limitsSelVoxelsMatrix[ nIterations, 1 ] : limitsSelVoxelsMatrix[ nIterations, 2 ]
  selectedIdx <- selectedIdxAll[ selectedVoxelsPassTemp ]
  pRFArraySelected <- pRFArray[ selectedIdx, ]
  tsArraySelected <- tsArray[ selectedIdx, ]

  tempMessage <- sprintf('nIteration %d of %d ...', nIterations, dim(limitsSelVoxelsMatrix)[1] )
  print( tempMessage )
  
  modelFunction <- function( passIdx, selectedIdx, pRFArraySelected, tsArraySelected, completeStimuliSequence, rescale_var, trValue ) {
    
    library( pracma )
    library( abind )
    library( neuRosim )
    
    outputR2 <- rep(999, dim( pRFArraySelected )[1] )
    outputCoeffs <- array(999, c( dim( pRFArraySelected )[1], 2 ) )
    outputPredictedTs <- array(999, dim( tsArraySelected ) )
    outputDetrendedTs <- array(999, dim( tsArraySelected ) )
    
    xSpace <- seq( -10, 10, length.out=dim( completeStimuliSequence )[1] )
    ySpace <- seq( -7, 7, length.out=dim( completeStimuliSequence )[1] )
    stimuliSequenceMatrixDimensions <- c( dim(completeStimuliSequence)[1]*dim(completeStimuliSequence)[2], dim(completeStimuliSequence)[3] )
    stimuliSequenceMatrix <- array( completeStimuliSequence, stimuliSequenceMatrixDimensions )
    stimuliTime <- seq( 0.25, dim( completeStimuliSequence )[3]/4, 0.25 )
    mriTime <- seq( trValue[1,1], dim( completeStimuliSequence )[3]/4, trValue[1,1] )
    hrf <- canonicalHRF( seq(0,30,0.25), verbose=FALSE ); #plot(hrf~seq(0,30,samplingTime))
    for ( nVoxel in 1:dim(pRFArraySelected)[1] ) { #nVoxel <- 1
      
      pRFXTemp <- pRFArraySelected[nVoxel,1]
      pRFYTemp <- pRFArraySelected[nVoxel,2]
      pRFSizeTemp <- pRFArraySelected[nVoxel,3]
      a <- dnorm( xSpace, pRFXTemp, pRFSizeTemp )
      b <- dnorm( ySpace, pRFYTemp, pRFSizeTemp )
      imgPrf <- a%*%t(b) # this is a 2d image generated by the matrix multiplication or crossproduct of two vectors (a,b, two gaussians), our pRF, a 2d gaussian
      imgPrf <- rescale_var( imgPrf, 0, 1)
      imgPrf[ imgPrf < 0.01 ] <- 0
      
      # compute time series in Kastner experiment expected from pRF only
      arrayImagePrf <- array( imgPrf )
      expectedTsFull <- as.double( array( arrayImagePrf%*%stimuliSequenceMatrix ) )
      expectedTsFull_upsampled_whole <- conv( expectedTsFull, hrf )
      expectedTsFull_upsampled_cut <- expectedTsFull_upsampled_whole[1:length(stimuliTime)]
      relevantQuantiles <- quantile( expectedTsFull_upsampled_cut, c(0.1,0.9) )
      expectedTsFull_upsampled_cut[ expectedTsFull_upsampled_cut <= relevantQuantiles[1] ] <- relevantQuantiles[1]
      expectedTsFull_upsampled_cut[ expectedTsFull_upsampled_cut >= relevantQuantiles[2] ] <- relevantQuantiles[2]
      expectedTsFull_downsampled <- interp1( x=stimuliTime, y=expectedTsFull_upsampled_cut, xi=mriTime, method=c('linear') )
      expectedTsFull_downsampled_scaled <- rescale_var( expectedTsFull_downsampled, 0, 1 )
      
      # kastner time serie
      tsArrayTemp <- tsArraySelected[nVoxel,]
      
      # detrend
      linearTrend <- seq( -1, 1, length.out=length(mriTime) )
      quadTrend <- rescale_var( linearTrend^2, -1, 1 )
      cubicTrend <- rescale_var( linearTrend^3, -1, 1 )
      dMat <- cbind( rep( 1, length( linearTrend ) ), linearTrend, quadTrend, cubicTrend )      
      a <- solve( crossprod(dMat,dMat), crossprod(dMat,tsArrayTemp) ) #beta coefficients (2: intercept and slope)
      expectedTs <- crossprod( t(dMat), a ) #expected ts
      residualsTs <- tsArrayTemp - expectedTs 
      
      # fit expected ts from pRF to residuals
      dMat <- cbind( rep( 1, length( expectedTsFull_downsampled_scaled ) ), expectedTsFull_downsampled_scaled )      
      a <- solve( crossprod(dMat,dMat), crossprod(dMat,residualsTs) ) #beta coefficients (2: intercept and slope)
      expectedTs_from_pRF <- crossprod( t(dMat), a ) #expected ts
      residualsSquared <- ( residualsTs - expectedTs_from_pRF )^2 
      ssRes <- sum( residualsSquared )
      ssTot <- sum( (residualsTs-mean(residualsTs))^2 )
      r2 <- 1-(ssRes/ssTot) #r squared
      
      outputDetrendedTs[ nVoxel, ] <- residualsTs
      outputPredictedTs[ nVoxel, ] <- expectedTs_from_pRF
      outputR2[ nVoxel ] <- r2
      outputCoeffs[ nVoxel, ] <- t( a )
      
    }
    
    # concatenate: the original indexews (to check), outputR, detrendedT2 and PredictedT2 per individual voxel 
    outMat <- cbind( selectedIdx, outputR2, outputCoeffs, outputDetrendedTs, outputPredictedTs )
    return( outMat ) 
    
  }
  
  # run example:
  #outModel <- modelFunction( 1, 
  #               tsArray, 
  #               pRFArray, 
  #               completeStimuliSequence, 
  #               rescale_var, trValue, 
  #               limitsSelVoxelsMatrix, 
  #               selectedIdxAll )
  
  # run code with or without parallelization
  runParallel <- 1
  if ( runParallel == 0 ){
    storeTimePar <- system.time( outputAll <- lapply( 1, modelFunction, 
                                                      selectedIdx=selectedIdx,
                                                      pRFArraySelected=pRFArraySelected,
                                                      tsArraySelected=tsArraySelected,
                                                      completeStimuliSequence=completeStimuliSequence,
                                                      rescale_var=rescale_var,
                                                      trValue=trValue ) )
    print( storeTimePar )
  }
  if ( runParallel == 1 ){
    print( sprintf('parallel fitting...' ) )
    print( sprintf('total passIdx = %d', dim( limitsSelVoxelsMatrix )[1] ) )
    nCores <- 8
    cl <- makeCluster(nCores, type='PSOCK')
    showConnections()
    storeTimePar <- system.time( outputAll <- parLapply( cl, 1, 
                                                         modelFunction, 
                                                         selectedIdx=selectedIdx,
                                                         pRFArraySelected=pRFArraySelected,
                                                         tsArraySelected=tsArraySelected,
                                                         completeStimuliSequence=completeStimuliSequence,
                                                         rescale_var=rescale_var,
                                                         trValue=trValue ) )
    stopCluster(cl)
    showConnections()
    print( storeTimePar )
  }
  
  #store output
  if ( nIterations==1 ) { outputMat <- outputAll[[1]] }
  if ( nIterations>1 ) { outputMat <- rbind( outputMat, outputAll[[1]] ) }
  
}

#rearrange output
selectedIdx <- outputMat[,1]
outputDetrendedTs <- outputMat[,5:166] 
outputPredictedTs <- outputMat[,167:328]
outputR2 <- outputMat[,2:4] # R2, intercept and slope

# get back original TS to check the everything is aligned correctly
tsArraySelected <- tsArray[ selectedIdx, ]

#save original ts, to make sure everything is in the correct space and place
outputFileName <- sprintf( '%s/%s_%s_originalTs_pRFBasedPrediction.nii.gz', outputFolder, participantName, inputDirection )
sprintf( 'writing file: %s ...', outputFileName )
tsArrayEmpty <- array( 0, dim( tsArray ) )
tsArrayEmpty[ selectedIdx, ] <- tsArraySelected 
storeOriginalTs <- array( tsArrayEmpty, c( dim(tsFile$brk) ) )
write.AFNI( outputFileName, brk = storeOriginalTs,
            origin = tsFile$origin, orient = tsFile$orient,
            defhead = tsFile$NI_head)


#save detrended ts
outputFileName <- sprintf( '%s/%s_%s_detrendedTs_pRFBasedPrediction.nii.gz', outputFolder, participantName, inputDirection )
sprintf( 'writing file: %s ...', outputFileName )
tsArrayEmpty <- array( 0, dim( tsArray ) )
tsArrayEmpty[ selectedIdx, ] <- outputDetrendedTs 
storeOriginalTs <- array( tsArrayEmpty, c( dim(tsFile$brk) ) )
write.AFNI( outputFileName, brk = storeOriginalTs,
            origin = tsFile$origin, orient = tsFile$orient,
            defhead = tsFile$NI_head)

#save predicted ts
outputFileName <- sprintf( '%s/%s_%s_predictedTs_pRFBasedPrediction.nii.gz', outputFolder, participantName, inputDirection )
sprintf( 'writing file: %s ...', outputFileName )
tsArrayEmpty <- array( 0, dim( tsArray ) )
tsArrayEmpty[ selectedIdx, ] <- outputPredictedTs 
storeOriginalTs <- array( tsArrayEmpty, c( dim(tsFile$brk) ) )
write.AFNI( outputFileName, brk = storeOriginalTs,
            origin = tsFile$origin, orient = tsFile$orient,
            defhead = tsFile$NI_head)

#save r2 values
outputFileName <- sprintf( '%s/%s_%s_R2_pRFBasedPrediction.nii.gz', outputFolder, participantName, inputDirection )
sprintf( 'writing file: %s ...', outputFileName )
tsArrayEmpty <- array( 0, c( dim( tsArray )[1], 3 ) )
tsArrayEmpty[ selectedIdx, ] <- outputR2 
storeOriginalTs <- array( tsArrayEmpty, c( dim(tsFile$brk)[1:3], 3 ) )
write.AFNI( outputFileName, brk = storeOriginalTs,
            origin = tsFile$origin, orient = tsFile$orient,
            defhead = tsFile$NI_head)

visualizeFlag <- 0
if (visualizeFlag==1) {
  #### try things quickly when debugging ####
  
  # define space, time and hrf
  xSpace <- seq( -10, 10, length.out=dim( completeStimuliSequence )[1] )
  ySpace <- seq( -7, 7, length.out=dim( completeStimuliSequence )[1] )
  stimuliTime <- seq( 0.25, dim( completeStimuliSequence )[3]/4, 0.25 )
  mriTime <- seq( trValue[1,1], dim( completeStimuliSequence )[3]/4, trValue[1,1] )
  hrf <- canonicalHRF( seq(0,30,0.25), verbose=FALSE ); #plot(hrf~seq(0,30,samplingTime))
  
  # define pRF
  nVoxel <- 4
  
  pRFXTemp <- pRFArraySelected[nVoxel,1]
  pRFYTemp <- pRFArraySelected[nVoxel,2]
  pRFSizeTemp <- pRFArraySelected[nVoxel,3]
  a <- dnorm( xSpace, pRFXTemp, pRFSizeTemp )
  b <- dnorm( ySpace, pRFYTemp, pRFSizeTemp )
  imgPrf <- a%*%t(b) # this is a 2d image generated by the matrix multiplication or crossproduct of two vectors (a,b, two gaussians), our pRF, a 2d gaussian
  imgPrf <- rescale_var( imgPrf, 0, 1)
  imgPrf[ imgPrf < 0.01 ] <- 0
  
  # compute time series in Kastner experiment expected from pRF only
  arrayImagePrf <- array( imgPrf )
  expectedTsFull <- as.double( array( arrayImagePrf%*%stimuliSequenceMatrix ) )
  expectedTsFull_upsampled_whole <- conv( expectedTsFull, hrf )
  expectedTsFull_upsampled_cut <- expectedTsFull_upsampled_whole[1:length(stimuliTime)]
  relevantQuantiles <- quantile( expectedTsFull_upsampled_cut, c(0.1,0.9) )
  expectedTsFull_upsampled_cut[ expectedTsFull_upsampled_cut <= relevantQuantiles[1] ] <- relevantQuantiles[1]
  expectedTsFull_upsampled_cut[ expectedTsFull_upsampled_cut >= relevantQuantiles[2] ] <- relevantQuantiles[2]
  expectedTsFull_downsampled <- interp1( x=stimuliTime, y=expectedTsFull_upsampled_cut, xi=mriTime, method=c('linear') )
  expectedTsFull_downsampled_scaled <- rescale_var( expectedTsFull_downsampled, 0, 1 )
  
  # kastner time serie
  tsArrayTemp <- tsArraySelected[nVoxel,]
  linearTrend <- seq( -1, 1, length.out=length(mriTime) )
  quadTrend <- rescale_var( linearTrend^2, -1, 1 )
  cubicTrend <- rescale_var( linearTrend^3, -1, 1 )
  
  # detrend
  dMat <- cbind( rep( 1, length( linearTrend ) ), linearTrend, quadTrend, cubicTrend )      
  a <- solve( crossprod(dMat,dMat), crossprod(dMat,tsArrayTemp) ) #beta coefficients (2: intercept and slope)
  expectedTs <- crossprod( t(dMat), a ) #expected ts
  residualsTs <- tsArrayTemp - expectedTs 
  
  # fit expected ts from pRF to residuals
  dMat <- cbind( rep( 1, length( expectedTsFull_downsampled_scaled ) ), expectedTsFull_downsampled_scaled )      
  a <- solve( crossprod(dMat,dMat), crossprod(dMat,residualsTs) ) #beta coefficients (2: intercept and slope)
  expectedTs_from_pRF <- crossprod( t(dMat), a ) #expected ts
  residualsSquared <- ( residualsTs - expectedTs_from_pRF )^2 
  ssRes <- sum( residualsSquared )
  ssTot <- sum( (residualsTs-mean(residualsTs))^2 )
  r2 <- 1-(ssRes/ssTot) #r squared
  
  # plotting just in case
  par(mfrow=c(3,4))
  image( xSpace, ySpace, completeStimuliSequence[,,7], zlim=c(0,1), las=1, col=gray( seq( 0,1,0.01 ) ) )
  image( xSpace, ySpace, imgPrf, zlim=c(0,max(array(imgPrf))), las=1, col=gray( seq( 0,1,0.01 ) ) )
  plot( expectedTsFull~stimuliTime, type='l', bty='n', las=1, lwd=2 )
  plot( expectedTsFull_upsampled_cut~stimuliTime, type='l', bty='n', las=1, lwd=2 )
  plot( expectedTsFull_downsampled~mriTime, type='l', bty='n', las=1, lwd=2 )
  plot( expectedTsFull_downsampled_scaled~mriTime, type='l', bty='n', las=1, lwd=2 )
  plot( tsArrayTemp~mriTime, type='l', bty='n', las=1, lwd=2 )
  plot( linearTrend~mriTime, type='l', bty='n', las=1, lwd=2 )
  lines( mriTime, quadTrend, lwd=2  )
  lines( mriTime, cubicTrend, lwd=2  )
  plot( expectedTs~mriTime, type='l', bty='n', las=1, lwd=2 )
  plot( residualsTs~mriTime, type='l', bty='n', las=1, lwd=2 )
  plot( residualsTs~mriTime, type='l', bty='n', las=1, lwd=2 )
  lines( mriTime, expectedTs_from_pRF, lwd=2, col='red'  )
  r2
  
}