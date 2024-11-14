debug <- 0

if (debug==0) {
  args <- commandArgs(T)
  print( args )
  mainDir <- getwd()
  setwd( mainDir )
}
if (debug==1) {
  rm( list=ls() ); gc();
  setwd('/analyse/Project0226/tests/nilsTest/data/20230622_SHI27')
  args <- c('ppp_epi_smooth.nii.gz', 'ppp_epi_mask.nii.gz', 'ppp_EPI_predictions/', 'ppp_modelOutput', '5')
  mainDir <- getwd()
  setwd( mainDir )
}

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
nParametersModel <- as.numeric( args[5] )

print( sprintf('input file = %s', inputFile ) )
print( sprintf('input roi = %s', inputRoi ) )
print( sprintf('inputDir = %s', inputDir ) )
print( sprintf('outputName = %s', outputName ) )
print( sprintf('nParametersModel = %d', nParametersModel ) )

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
fittingPrf <- function( idxFile, predictionFiles, tsTransposedSel, nParametersModel ) {
  
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
  
  voxelModel <- function( passIdx, tsTransposedSel, limitsSelVoxelsMatrix, tsPrediction, predictionGrid, nParametersModel ) { #this fits the model on serveral voxels at a time, see limitsSelVoxelsMatrix
    
    # passIdx <- 1
    library( abind )
    
    selTsVoxel <- tsTransposedSel[ , limitsSelVoxelsMatrix[passIdx,1]:limitsSelVoxelsMatrix[passIdx,2] ]
    selTsVoxelMean <- apply( selTsVoxel, 2, mean ) #voxels average
    selTsVoxelMean <- matrix( rep( selTsVoxelMean, dim(selTsVoxel)[1] ), nrow = dim(selTsVoxel)[1], byrow = TRUE )
    ssTot <- apply( (selTsVoxel-selTsVoxelMean)^2, 2, sum) #voxels total sum of squares
    
    runLinMod <- function( nIndex, nParametersModel ) { #get best fit for every prediction (nIndex = counter of predictions)
      dMat <- cbind( tsPrediction[nIndex,] )      
      dMat01 <- cbind( rep(1,length(dMat)), dMat ) #column of ones and column of predictor
      a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope)
      expectedTs <- crossprod( t(dMat01), a ) #expected ts
      residualsSquared <- ( selTsVoxel - expectedTs )^2 
      ssRes <- apply( residualsSquared, 2, sum )
      r2 <- 1-(ssRes/ssTot) #r squares
      
      kParameters <- nParametersModel
      nDataPoints <- length( expectedTs )
      r2Adj <- 1 - (1 - r2) * (nDataPoints - 1) / (nDataPoints - kParameters - 1) # adjusted_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - k - 1)
      AIC <- nDataPoints * log(ssRes / nDataPoints) + 2 * kParameters  #aic_value <- n * log(ssRes / n) + 2 * k, where n is the number of datapoints and k is the number of parameters
      
      return( rbind( r2, a, r2Adj, AIC, expectedTs ) ) # added 2 more parameters: ssRes, ssTot
      
      #return( rbind( r2, a, logLikelihood, AIC, expectedTs ) ) # added 2 more parameters: logLik, AIC
    }  
    
    outVoxelList <- lapply( 1:dim(tsPrediction)[1], runLinMod, nParametersModel=nParametersModel  ) #apply the function for all predictions, get a list as output
    outVoxel3D <- abind( outVoxelList, along=3 )  #reshape the list in a 3dmatrix with fields: output X voxels tested X predictions
    r2Matrix <- outVoxel3D[1,,]	
    betaPositiveMatrix <-  array( TRUE, dim( r2Matrix ) ) #outVoxel3D[3,,] > 0 # set it to positive;
    
    extractData <- function(nSelectedVoxel) { #nSelectedVoxel <- 1
      indexBetaZero <- betaPositiveMatrix[nSelectedVoxel,]
      if ( sum(indexBetaZero)>0 ) {
        indexBetaZero <- which( indexBetaZero )
        indexVarExp <- which.max( r2Matrix[nSelectedVoxel,indexBetaZero] )[1]
        return( as.numeric( c( predictionGrid[indexBetaZero[indexVarExp],], outVoxel3D[,nSelectedVoxel,indexBetaZero[indexVarExp]] ) ) ) 
      }
      if ( sum(indexBetaZero)==0 ) {
        return( rep(0, dim(predictionGrid)[2]+5+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 3=r2,beta intercept, beta slope # + changed 3 to 5 (LL, AIC)
      }
    }	
    
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }

  #test
  #outModelTest <- lapply( 1:3, voxelModel, 
  #   tsTransposedSel=tsTransposedSel, 
  #   limitsSelVoxelsMatrix=limitsSelVoxelsMatrix, 
  #   tsPrediction=tsPrediction, 
  #   predictionGrid=predictionGrid,
  #   nParametersModel=nParametersModel )
  
  print( sprintf('parallel fitting: %1.0f of %1.0f ...', idxFile, length( predictionFiles ) ) )
  detectCores()
  nCores <- 6
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( outModel <- parLapply(cl, runIndex, voxelModel, 
                                                     tsTransposedSel=tsTransposedSel, 
                                                     limitsSelVoxelsMatrix=limitsSelVoxelsMatrix, 
                                                     tsPrediction=tsPrediction,
                                                     predictionGrid=predictionGrid,
                                                     nParametersModel=nParametersModel ) )
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
  
  outModel <- fittingPrf( nFile, predictionFiles=predictionFiles, tsTransposedSel=tsTransposedSel, nParametersModel=nParametersModel )

  naIdx <- which( is.na( outModel[ (nParameters+1) , ] ) )
  outModel[ (nParameters+1) , naIdx ] <- 0
  
  if (nFile==1) { outModelLoop <- outModel }
  if (nFile>1) { 
    selectedCols <- outModel[ (nParameters+1) ,] > outModelLoop[ (nParameters+1), ] #where r2 is stored 
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
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

storeAllPred <- t(outModelLoop[ seq(1,(nParameters+5)), ]) #changed 3 to 5

### test from here to remove surround and remove the phase and eccentricity conversion

storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ (nParameters+6):dim(outModelLoop)[1], ]
storeAllPredOut <- array( 0, c( (nParameters+5), dim(tsTransposedAll)[2] ) ) 

print('save linear step...')
storeAllPredOut[,selIdxVoxel] <- t( storeAllPred )
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
labels <- c('hrf_a1','hrf_a2', 'hrf_b1', 'hrf_b2', 'hrf_c', 'varExp', 'intercept','slope', 'r2Adj', 'AIC' ) 
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}


# clean up
rm( list=ls() ); gc()


