args <- commandArgs(T)
print( args )
mainDir <- getwd()
setwd( mainDir )

# to debug
#rm( list=ls() ); gc();
#setwd('/analyse/Project0226/KastnerModel/data_KastnerModel/AFH28_201119')
#args <- c('zzz_KASTMDL_stc_mc_tu_dt_avg_blur_smooth.nii.gz', 'zzz_grey_kastner.nii.gz', 'zzz_kastner_model/', 'zzz_kastner_model', 'kinematicModel')
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
inputRoi <- args[2]
inputDir <- args[3]
outputName <-  args[4]
modelType <- args[5]

print( sprintf('input file = %s', inputFile ) )
print( sprintf('input roi = %s', inputRoi ) )
print( sprintf('inputDir = %s', inputDir ) )
print( sprintf('outputName = %s', outputName ) )
print( sprintf('modelType = %s', modelType ) )

# load data
ts <- read.AFNI( inputFile )
tsArray <- array( ts$brk, c( prod( dim( ts$brk )[1:3] ), dim(ts$brk)[4] ) )

# load roi
roi <- read.AFNI( inputRoi )
roiVol <- roi$brk

# ts matrix
indexVol <- roiVol[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
selIdxVoxel <- which( indexArray == 1 )
tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series

# clean up
rm( ts, roi, tsArray, roiVol, tsTransposedAll ); gc()

# prediction files
setwd( inputDir )
predictionFiles <- dir( pattern = '*.RData' )
setwd(mainDir)

# split number of voxels in chunks
nVoxels <- dim(tsTransposedSel)[2]
totalVoxelsIterations <- ceil( nVoxels / 200 )
limitsSelVoxels <- round( seq(1, nVoxels, length.out=totalVoxelsIterations) ) #split between chuncks for parallel fitting (about 200 voxels each)
limitsSelVoxels[1] <- 1
limitsSelVoxels[ length(limitsSelVoxels) ] <- nVoxels
limitsSelVoxelsMatrix <- array( 0, c( (length(limitsSelVoxels)-1) , 2 ) )
limitsSelVoxelsMatrix[,1] <- limitsSelVoxels[1:(length(limitsSelVoxels)-1)]
limitsSelVoxelsMatrix[,2] <- limitsSelVoxels[2:(length(limitsSelVoxels))]
limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] <- limitsSelVoxelsMatrix[2:dim(limitsSelVoxelsMatrix)[1],1] + 1 #matrix with starting ts and ending ts for each chunk
runIndex <- seq( 1:dim(limitsSelVoxelsMatrix)[1] )

# define function
fittingEyeKinematics <- function( idxFile, predictionFiles, tsTransposedSel ) {
  
  load( sprintf('%s/%s', inputDir, predictionFiles[ idxFile ] ) )
  nTRs <- dim( tsTransposedSel )[1] 
  nParameters <- dim( outMatrix )[1]-nTRs
  
  predictionGrid <- t( outMatrix[ 1:nParameters, ] )
  tsPrediction <- t( outMatrix[ (nParameters+1):dim(outMatrix)[1], ] )
  controlPredictions <- apply( abs( tsPrediction ), 1, sum ) > 0.0001 & apply( is.na(tsPrediction), 1, sum ) == 0
  tsPrediction <- tsPrediction[controlPredictions,]
  predictionGrid <- predictionGrid[controlPredictions,]
  print( sprintf('prediction file: %s...', predictionFiles[ idxFile ] ) )
  print( sprintf('number of predictions to fit: %1.0f...', dim( predictionGrid )[1] ) )
  
  voxelModel <- function( passIdx, tsTransposedSel, limitsSelVoxelsMatrix, tsPrediction, predictionGrid ) { #this fits the model on serveral voxels at a time, see limitsSelVoxelsMatrix
    
    library( abind )
    selTsVoxel <- tsTransposedSel[ , limitsSelVoxelsMatrix[passIdx,1]:limitsSelVoxelsMatrix[passIdx,2] ]
    selTsVoxelMean <- apply( selTsVoxel, 2, mean ) #voxels average
    selTsVoxelMean <- matrix( rep( selTsVoxelMean, dim(selTsVoxel)[1] ), nrow = dim(selTsVoxel)[1], byrow = TRUE )
    ssTot <- apply( (selTsVoxel-selTsVoxelMean)^2, 2, sum) #voxels total sum of squares
    
    runLinMod <- function( nIndex ) { #get best fit for every prediction (nIndex = counter of predictions)
      dMat <- cbind( tsPrediction[nIndex,] )      
      dMat01 <- cbind( rep(1,length(dMat)), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope)
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      
      r2 <- 1-(ssRes/ssTot) #r squares
      #r2 <- diag( cor( expectedTs, selTsVoxel, method=c('spearman') ) )
      
      #if (modelType=='kinematicModel'){ adjR2 <- 1 - ( ( (1-r2)*(dim(selTsVoxel)[1]-1) ) / ( dim(selTsVoxel)[1]-4-1 ) ) }
      #if (modelType=='baselineModel'){ adjR2 <- 1 - ( ( (1-r2)*(dim(selTsVoxel)[1]-1) ) / ( dim(selTsVoxel)[1]-2-1 ) ) }
      logLikelihood <- rep(0,length(r2))
      AIC <- rep(0,length(r2))
      adjR2 <- rep(0,length(r2))
      #varResid <- apply( (selTsVoxel - expectedTs), 2, var )
      #logLikelihood <- -( selTsVoxel*log(2*pi*varResid) ) / 2 + sum( ( -residualsSquared) / (2*varResid) )
      #outModelLoop[15, allTs] <- logLikelihood
      #AIC  <- 2*nParams - 2*logLikelihood
      
      return( rbind( r2, a, logLikelihood, AIC, adjR2, expectedTs ) ) # added 3 more parameters: logLik, AIC, adjR2 
    }  
    
    outVoxelList <- lapply( 1:dim(tsPrediction)[1], runLinMod  ) #apply the function for all predictions, get a list as output
    outVoxel3D <- abind( outVoxelList, along=3 )  #reshape the list in a 3dmatrix with fields: output X voxels tested X predictions
    betaPositiveMatrix <- outVoxel3D[3,,] > 0
    r2Matrix <- outVoxel3D[1,,]	
    
    extractData <- function(nSelectedVoxels) {
      indexBetaZero <- betaPositiveMatrix[nSelectedVoxels,]
      if ( sum(indexBetaZero)>0 ) {
        indexBetaZero <- which( indexBetaZero )
        indexVarExp <- which.max( r2Matrix[nSelectedVoxels,indexBetaZero] )[1]
        return( as.numeric( c( predictionGrid[indexBetaZero[indexVarExp],], outVoxel3D[,nSelectedVoxels,indexBetaZero[indexVarExp]] ) ) ) 
      }
      if ( sum(indexBetaZero)==0 ) {
        return( rep(0, dim(predictionGrid)[2]+6+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 3=r2,beta intercept, beta slope # + changed 3 to 6 (LL, AIC, adjR2)
      }
    }	
    
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }
  
  # parallel fitting
  print( sprintf('parallel fitting: %1.0f of %1.0f ...', idxFile, length( predictionFiles ) ) )
  detectCores()
  nCores <- 12
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( outModel <- parLapply( cl, runIndex, voxelModel, 
                                                     tsTransposedSel=tsTransposedSel, 
                                                     limitsSelVoxelsMatrix=limitsSelVoxelsMatrix, 
                                                     tsPrediction=tsPrediction,
                                                     predictionGrid=predictionGrid ) )
  showConnections()
  stopCluster(cl)
  showConnections()
  print( storeTimePar )
  rm(cl); gc();
  #print( sprintf('fitting: %1.0f of %1.0f ...', idxFile, length( predictionFiles ) ) )
  #storeTimePar <- system.time( outModel <- lapply( runIndex, voxelModel ) )
  #print( storeTimePar )
  
  outMatrix <- array(0, c( dim(outModel[[1]])[1], length(selIdxVoxel) ) )
  for ( nElements in 1:length(outModel) ) {
    tempNElementsStart <- limitsSelVoxelsMatrix[nElements,1]
    tempNElementsEnd <- limitsSelVoxelsMatrix[nElements,2]
    outMatrix[,tempNElementsStart:tempNElementsEnd] <- outModel[[nElements]]
  }
  
  return( outMatrix )
  
}

# loop over files
for ( nFile in 1:length(predictionFiles) ) {
  
  load( sprintf('%s/%s', inputDir, predictionFiles[ nFile ] ) )
  nTRs <- dim( tsTransposedSel )[1] 
  nParameters <- dim( outMatrix )[1]-nTRs
  
  outModel <- fittingEyeKinematics( nFile, predictionFiles=predictionFiles, tsTransposedSel=tsTransposedSel )
  
  if (nFile==1) { outModelLoop <- outModel }
  if (nFile>1) { 
    selectedCols <- outModel[ (nParameters+1) ,] > outModelLoop[ (nParameters+1), ] #where r2 is stored 
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
} 

####################################
# log likelihood + AIC , now adjR2 #
####################################
print( 'computing log-lihelihood, AIC and adjR2...' )

if (modelType=='kinematicModel'){ nParams <- 4 }
if (modelType=='kinematic_baselineModel'){ nParams <- 4 }
if (modelType=='baselineModel'){ nParams <- 2}

expectedWinTs <- outModelLoop[ (nParameters+7):dim(outModelLoop)[1], ]
#resSquared <- (tsTransposedSel-expectedWinTs)^2
#varResid <- apply(resSquared,2,var)
tsArrayTime <- dim(tsTransposedSel)[1]

counterTrack <- 1
keepTrackVisually <- floor( seq(1, dim(tsTransposedSel)[2], length.out = 20) )

for (allTs in 1:dim(tsTransposedSel)[2]){
  #this is here because there were empty time series (?)

  if ( outModelLoop[ nParameters+1 , allTs ] == 0 ) {
    outModelLoop[12, allTs] <- 0 #where log-likelihood is stored
    outModelLoop[13, allTs] <- 0 #where AIC is stored
    outModelLoop[14, allTs] <- 0 #where adjR2 is stored
  }
  if ( outModelLoop[ nParameters+1 , allTs ] != 0) {
    
    mod4AIC <- lm( tsTransposedSel[, allTs] ~ expectedWinTs[, allTs] )
    logLik_current <- as.numeric( logLik( mod4AIC ) )
    AIC_current <- AIC( mod4AIC )
    r2 <- outModelLoop[9, allTs]
    adjR2 <- 1 - ( ( (1-r2)*(tsArrayTime-1) ) / ( tsArrayTime-nParams-1 ) ) 
    #logLikelihood <- -( tsArrayTime*log(2*pi*varResid[allTs]) ) / 2 + sum( ( -resSquared[,allTs]) / (2*varResid[allTs]) )
    #AIC  <- 2*nParams - 2*logLikelihood
    outModelLoop[12, allTs] <- logLik_current
    outModelLoop[13, allTs] <- AIC_current
    outModelLoop[14, allTs] <- adjR2
  }
  
  if ( allTs>keepTrackVisually[counterTrack] ) { 
    print( sprintf('computing log-lik, AIC and adjR2, iteration: %1.0f/%1.0f', counterTrack+1, length( keepTrackVisually ) ) )
    counterTrack <- counterTrack + 1
  }
  
}

#saving output:

# load data
ts <- read.AFNI( inputFile )
tsArray <- array( ts$brk, c( prod( dim( ts$brk )[1:3] ), dim(ts$brk)[4] ) )

# load roi
roi <- read.AFNI( inputRoi )
roiVol <- roi$brk

# ts matrix
indexVol <- roiVol[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
selIdxVoxel <- which( indexArray == 1 )
tsTransposedSel <- tsTransposedAll[,selIdxVoxel] #all selected time series

storeAllPred <- t(outModelLoop[ seq(1,(nParameters+6)), ]) #changed 3 to 6
storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ (nParameters+7):dim(outModelLoop)[1], ]
storeAllPredOut <- array( 0, c((nParameters+6), dim(tsTransposedAll)[2] ) ) 

fileDetrended <- sprintf('%s_DetrendedTs.nii.gz',outputName) 
fileTs <- sprintf('%s_PredixtedTs.nii.gz',outputName) 
fileParams <- sprintf('%s_params.nii.gz',outputName)

#save expected ts
rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( fileTs, brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)

# save parameters
storeAllPredOut[,selIdxVoxel] <- t( storeAllPred )
rParameters <- array( storeAllPredOut, c( dim(storeAllPredOut)[1], dim(ts$brk)[1:3] ) )
rParameters <- aperm( rParameters, c(2,3,4,1) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', inputRoi )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( fileParams, brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')

# save original ts (for naming purposes)
write.AFNI( fileDetrended, brk = ts$brk,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)

# relabel
labels <- c('hrf_a1','hrf_a2','hrf_b1','hrf_b2','hrf_c','mu','kappa','nLin','varExp','intercept','slope','logLikelihood', 'AIC', 'adjR2')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}


# clean up
rm( list=ls() ); gc()
