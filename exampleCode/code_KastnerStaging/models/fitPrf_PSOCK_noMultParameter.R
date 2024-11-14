args <- commandArgs(T)
print( args )
mainDir <- getwd()
setwd( mainDir )

# to debug
#rm( list=ls() ); gc();
#setwd('/analyse/Project0226/KastnerModel/data_KastnerClassic_pRF/ASM16_BarsEyes')
#args <- c('ppp_meanTs_bars_topUp_dt_blur_sm.nii.gz', 'ppp_grey_bars.nii.gz', 'ppp_savedPredictions_prfSimple/', 'ppp_prf_nocss','0','0')
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
prfType <-   as.numeric( args[5] )
flagSurround <-   as.numeric( args[6] )

print( sprintf('input file = %s', inputFile ) )
print( sprintf('input roi = %s', inputRoi ) )
print( sprintf('inputDir = %s', inputDir ) )
print( sprintf('outputName = %s', outputName ) )
print( sprintf('prfType = %1.0f', prfType ) )
print( sprintf('flagSurround = %1.0f', flagSurround ) )

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
fittingPrf <- function( idxFile, predictionFiles, tsTransposedSel ) {
  
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
      logLikelihood <- rep(0,length(r2))
      AIC <- rep(0,length(r2))
      return( rbind( r2, a, logLikelihood, AIC, expectedTs ) ) # added 2 more parameters: logLik, AIC
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
        return( rep(0, dim(predictionGrid)[2]+5+dim(selTsVoxel)[1] ) ) #no positive betas, return array of zeros, 3=r2,beta intercept, beta slope # + changed 3 to 5 (LL, AIC)
      }
    }	
    
    outModel <- sapply( 1:dim(r2Matrix)[1], extractData )
    return( outModel )
  }

  #test
  #outModelTest <- lapply( 1:3, voxelModel, tsTransposedSel=tsTransposedSel, limitsSelVoxelsMatrix=limitsSelVoxelsMatrix, tsPrediction=tsPrediction )
  print( sprintf('parallel fitting: %1.0f of %1.0f ...', idxFile, length( predictionFiles ) ) )
  detectCores()
  nCores <- 12
  cl <- makeCluster(nCores, type='PSOCK')
  showConnections()
  storeTimePar <- system.time( outModel <- parLapply(cl, runIndex, voxelModel, 
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
  
  outModel <- fittingPrf( nFile, predictionFiles=predictionFiles, tsTransposedSel=tsTransposedSel )
  
  if (nFile==1) { outModelLoop <- outModel }
  if (nFile>1) { 
    selectedCols <- outModel[ (nParameters+1) ,] > outModelLoop[ (nParameters+1), ] #where r2 is stored 
    if ( sum( selectedCols ) > 0 ) {
      outModelLoop[,selectedCols] <- outModel[,selectedCols] 
    }
  }
  
} 

#log likelihood + AIC
# nParams <- 3 
# if (flagSurround==1){ 
#   nParams <- nParams + 2 }
# if (prfType==1){
#   nParams <- nParams + 3}
# 
# expectedWinTs <- outModelLoop[ (nParameters+6):dim(outModelLoop)[1], ]
# resSquared <- (tsTransposedSel-expectedWinTs)^2
# varResid <- apply(resSquared,2,var)
# 
# for (allTs in 1:length(varResid)){
#   #this is here because there were empty time series (?)
#   
#   if ( outModelLoop[ nParameters+1 , allTs ] == 0 ) {
#     outModelLoop[15, allTs] <- 0
#     outModelLoop[16, allTs] <- 0
#   }
#   
#   if ( outModelLoop[ nParameters+1 , allTs ] != 0) {
#     logLikelihood <- -( tsArrayTime*log(2*pi*varResid[allTs]) ) / 2 + sum( ( -resSquared[,allTs]) / (2*varResid[allTs]) )
#     outModelLoop[15, allTs] <- logLikelihood
#     AIC  <- 2*nParams - 2*logLikelihood
#     outModelLoop[16, allTs] <- AIC
#   }
#   
# }

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
#### extract surround size ####
#### extract surround size ####
print('surround size...')
if (flagSurround==1) { 
  print('get FWHM...')
  progress <- 0.05
  FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
  xFWHM <- seq(-75,75,0.0125)
  for (k in 1:dim(storeAllPred)[1]) {
    
    parameters <- storeAllPred[k,]
    
    if ( sum(parameters) != 0 ) {
      a <- dnorm( xFWHM, 0, parameters[3] )
      b <- dnorm( xFWHM, 0, parameters[4] )
      #r <- scaleData(a,1,0) - scaleData(b,1,0)*parameters[7]
      #parameters <- c(0,0,5,5)
      #a <- dnorm( xFWHM, 0, parameters[3] ); b <- dnorm( xFWHM, 0, parameters[4] );
      if ( parameters[3] != parameters[4] ) {
        r <- a - b
      }else{
        r <- a
      }
      
      maxIdx <- which.max(r)
      rHalf <- r[ maxIdx:length(xFWHM) ]
      xHalf <- xFWHM[ maxIdx:length(xFWHM) ]
      halfMax <- max(rHalf)/2
      xMid <- xHalf[ rHalf<=halfMax  ]  # x point at half max 
      FWHMcenter <- min( xMid )*2 #twice x point at half max
      if ( sum(r<0)>0 ) { #if there is a detectable surround, compute surround size
        rMin <- which.min(rHalf)
        FWHMsurround <- xHalf[rMin[1]]*2
      }
      if ( sum(r<0)==0 ) { #otherwise set it to zero
        FWHMsurround <- 0
      }
      FWHM[k,] <- c(FWHMcenter,FWHMsurround)
      #plot(rHalf~xHalf ); abline(v=min(xMid)); abline(v=xHalf[rMin])
      rm('rMin','xMid','FWHMcenter','FWHMsurround')
    }
    
    if (k>dim(storeAllPred)[1]*progress) {
      cat(paste('*',""))
      progress <- progress + 0.05
    }
  }
}

if (flagSurround==0) { 
  FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
  FWHM[,1] <- storeAllPred[,3]
  FWHM[,2] <- 1000
}

storeAllExpectedTs <- array( 0, dim(tsTransposedAll) )
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ (nParameters+6):dim(outModelLoop)[1], ]
storeAllPredOut <- array( 0, c((nParameters+9), dim(tsTransposedAll)[2] ) ) #changed 7 to 9

print('save linear step...')
polCoords <- cart2pol( storeAllPred[,c(2,1)] ) 
storeAllPredOut[,selIdxVoxel] <- t( cbind( storeAllPred, polCoords, FWHM ) )
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
labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','hrf_a1','hrf_a2','multpar','statNonLin','n_cycles', 'angle', 'phase', 'varExp','intercept','slope', 'logLikelihood', 'AIC', 'theta','radius','fwhmCenter','surroundSize') #added 4 more labels
for (k in 1:length(labels)) {t( cbind( storeAllPred, polCoords, FWHM ) )
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}


# clean up
rm( list=ls() ); gc()
