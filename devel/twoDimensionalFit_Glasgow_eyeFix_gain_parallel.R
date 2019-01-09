

args <- commandArgs(T)
print( args )

#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')
#args <- c('bars_mask_res.nii.gz','barsTs_res.nii.gz','outTest_delme.nii.gz','/packages/afni/17.0.13','1')

## prepare the stimuli with the portrait as well

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

mainDir <- getwd()
args <- c('meanTs_eyeMovement_topUp_res_mask.nii', 'meanTs_eyeMovement_topUp_res.nii', 'output_test_delme','/home/alessiof/abin', '1', '4','1','0.166','1')
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
#instr <- '3dresample -dxyz 4 4 4 -orient RAI -rmode Lin -prefix barsTs_res.nii.gz -inset meanTsBars_topUp.nii'
#system( instr )
#instr <- '3dTstat -mean -prefix bars_mean_res.nii.gz barsTs_res.nii.gz'
#system( instr )
#instr <- '3dAutomask -prefix bars_mask_res.nii.gz bars_mean_res.nii.gz'
#system( instr )
#instr <- '3dresample -dxyz 4 4 4 -orient RAI -rmode Lin -prefix eyeTs_res.nii.gz -inset meanTsEye_topUp.nii'
#system( instr )




source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( neuRosim )
library( parallel )

epiFile <- args[1]
inputFile <- args[2]
outSuffix <- args[3]
flagSurround <- as.numeric( args[5] )
fitIntercept <- as.numeric( args[7] )
polortArg <- as.numeric( args[6] )
samplingTime <- as.numeric( args[8] )
stimType <- as.numeric( args[9] )
fineFit <- 0

# load the data
# get orientation
print('detrend and load data...')
instr <- sprintf('3dinfo -orient %s > __tttt.1D', inputFile)
system( instr )
orientValue <- read.table('__tttt.1D')
system('rm __tttt.1D')
# get TR
instr <- sprintf('3dinfo -tr %s > __tttt.1D', inputFile)
system( instr )
trValue <- read.table('__tttt.1D')
system('rm __tttt.1D')

instr <- sprintf( '3dresample -prefix _ttt_orient.nii.gz -inset %s', inputFile )
system( instr)
instr <- sprintf('3dDetrend -polort %1.0f -prefix _ttt_data.nii.gz _ttt_orient.nii.gz', polortArg)
system( instr )
instr <- sprintf('3dresample -orient %s -prefix %s_detrendedData.nii.gz -inset _ttt_data.nii.gz', orientValue[1,1], outSuffix )
system( instr )
tsOriginal <- read.AFNI( '_ttt_orient.nii.gz' )
tsOriginalArray <- array( tsOriginal$brk, c( prod( dim( tsOriginal$brk )[1:3] ), dim(tsOriginal$brk)[4] ) )
ts <- read.AFNI( '_ttt_data.nii.gz' )
tsArray <- array( ts$brk, c( prod( dim( ts$brk )[1:3] ), dim(ts$brk)[4] ) )
system( 'rm _ttt_orient.nii.gz' )
system( 'rm _ttt_data.nii.gz' )

instr <- sprintf( '3dresample -prefix _ttt_mask.nii.gz -inset %s', epiFile )
system( instr)
meanEpi <- read.AFNI( '_ttt_mask.nii.gz' )
system( 'rm _ttt_mask.nii.gz' )


# load stimuli definition
print('get prfStimuli.RData...')
#load( file='prfStimuli.RData' )
#setwd('/analyse/Project0226/scanner_train - Copy/orientationgrating_points_fit_long_glas_large_get_stimuli')
#arrayStim <- scan( 'eyeMovingStim.txt' )
#setwd('/analyse/Project0226/scanner_train - Copy/orientationgrating_points_fit_long_glas_large_fixation_get_stimuli_border')
#arrayStim <- scan( 'eyeFixStim_border.txt' )

setwd('/analyse/Project0226/dataSummary')
if (stimType==1) { 
	arrayStim <- scan( 'eyeMovingStim.txt' )
	setwd(mainDir)
	stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==2) { 
	arrayStim <- scan( 'eyeFixStim.txt' )
	setwd(mainDir)
	stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==3) { 
	arrayStim <- scan( 'eyeFixStim_border.txt' )
	setwd(mainDir)
	stimMat <- aperm( array( arrayStim, c(200,1510,150) ), c( 3, 1, 2 ) ) # eye movement
}
if (stimType==4) { 
	arrayStim <- scan( 'prfStim.txt' )
	setwd(mainDir)
	stimMat <- aperm( array( arrayStim, c(200,1860,150) ), c( 3, 1, 2 ) ) # eye movement
}

#stimMat <- aperm( array( arrayStim, c(128,620,96) ), c( 1, 3, 2 ) )
#stimMat <- aperm( array( arrayStim, c(200,1860,150) ), c( 3, 1, 2 ) ) # prf
#stimMatInterp2 <- array( 0,c(60,60,1510) )
#xStim <- seq(0,200,length.out=200)	
#yStim <- seq(0,150,length.out=150)	
#xStimInt <- seq(0,200,length.out=60)	
#yStimInt <- seq(0,150,length.out=60)	
#for (interpCounter in 1:dim(stimMat)[2] ) {
#	imageLoop <- stimMat[,,interpCounter]
#	stimMatInterp2[,,interpCounter] <- interp2( xStim, yStim, imageLoop, xStimInt, yStimInt, method=c('linear') )
#}
#stimMat <- stimMatInterp2
stimMatFlip <- stimMat[ ,dim(stimMat)[2]:1, ]
#image( stimMatFlip[,,556], axes=FALSE )
#image( stimMat[,,56], axes=FALSE )
x11( width=3, height=3 )
for ( snap in 1:dim(stimMat)[3] ) {
  image( stimMat[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01) 
}
stimSeq <- stimMatFlip
x <- seq(-7,7,length.out = dim(stimSeq)[1] )
y <- seq(-10,10,length.out = dim(stimSeq)[2] )

#### get the standard hrf ####
print('get hrf...')
hrf <- canonicalHRF( seq(0,30,samplingTime) )


multVector <- c(0.7,0.8,1,1.1,1.2)
for (modelFitCounter in 1:length(multVector)) {

#### prepare to generate the predictions ####
addSpace <- abs( min(x) )*0.5
print('build prediction...')
xPosFit <- seq( min(x)-addSpace, max(x)+addSpace, length.out=6 ) * multVector[ modelFitCounter ]
yPosFit <- seq( min(y)-addSpace, max(y)+addSpace, length.out=6 ) * multVector[ modelFitCounter ]
sigmaArrayPositive <- seq( 0.25, 7, length.out=6 ) * multVector[ modelFitCounter ]

if (flagSurround==1) { sigmaArrayNegative <- sigmaArrayPositive }
if (flagSurround==0) { sigmaArrayNegative <- 1000 }
predictionGridTemp <- expand.grid( xPosFit, yPosFit, sigmaArrayPositive, sigmaArrayNegative )
keepPredictionIdx <- predictionGridTemp[ ,3] < predictionGridTemp[ ,4]
predictionGrid <- predictionGridTemp[ keepPredictionIdx, ]

#tsPrediction <- array( 0, c( dim(predictionGrid)[1], dim( ts$brk )[4] ) )
stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
incrementalCounter <- dim(predictionGrid)[1] * 0.05
timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(ts$brk)[4] )

#### function to generate the predictions ####
generatePrediction <- function( indexPrediction, inputPredictionGrid ) {
  prfPar <- as.numeric( inputPredictionGrid[indexPrediction,] )
  
  #prfPar <- c(-2,-4,1,1.1)
  a <- dnorm( x, prfPar[1], prfPar[3] )
  b <- dnorm( y, prfPar[2], prfPar[3] )
  #a <- normFun( x, prfPar[1], prfPar[3] )
  #b <- normFun( y, prfPar[1], prfPar[3] )
  imgCenter <- tcrossprod(a,b)
  #image( imgCenter )
  
  a <- dnorm( x, prfPar[1], prfPar[4] )
  b <- dnorm( y, prfPar[2], prfPar[4] )
  #a <- normFun( x, prfPar[1], prfPar[4] )
  #b <- normFun( y, prfPar[1], prfPar[4] )
  imgSurround <- tcrossprod(a,b)
  
  if (flagSurround==1) { r <- imgCenter - imgSurround }
  if (flagSurround==0) { r <- imgCenter }
  rMat <- array(r)

  #trMat <- t(stimSeqMat)
  #system.time( predictionLoop <- as.numeric( trMat%*%rMat ) )
  predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
  pConv <- conv( predictionLoop, hrf )
  pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
  tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('linear') )
  #tsPrediction[nPrediction,] <- tsPredictionMriInterp / max( tsPredictionMriInterp )
  #return( round( tsPredictionMriInterp / max( tsPredictionMriInterp ), 5 ) ) #### scale predictions to 1 and round them to 5 digits
  return( round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 ) ) #### scale predictions betweeo 0 and 1 and round them to 5 digits
}

#### generate predictions in parallel ####
library(parallel)
detectCores()
nCores <- 8
cl <- makeCluster(nCores, type='FORK')
storeTimePar <- system.time( tsPredictionTransposed <- parSapply(cl, 1:dim(predictionGrid)[1], generatePrediction, inputPredictionGrid=predictionGrid ) )
stopCluster(cl)
#storeTimeSerial <- system.time( sapply(1:100, generatePrediction, inputPredictionGrid=predictionGrid ) )
print( storeTimePar )
#print( storeTimeSerial )

#### clean up ts predictions and prediction grid ####
tsPrediction <- t( tsPredictionTransposed )
controlPredictions <- apply( tsPrediction, 1, sum ) != 0 
tsPrediction <- tsPrediction[controlPredictions,]
predictionGrid <- predictionGrid[controlPredictions,]
print( dim( predictionGrid ) )


#### linear fitting in parallel ####
print('linear fitting...')
indexVol <- meanEpi$brk[,,,1]
indexArray <- array( indexVol, prod( dim( indexVol ) ) )
tsTransposedAll <- t( tsArray )
selIdxVoxel <- which( indexArray == 1 )
tsTransposedSel <- tsTransposedAll[,selIdxVoxel] 
voxelModel <- function(passIdx) {
	selTsVoxel <- tsTransposedSel[ , passIdx ]
	ssTot <- sum((selTsVoxel-mean(selTsVoxel))^2)
	runLinMod <- function(nIndex) { 
      		dMat <- cbind( tsPrediction[nIndex,] )      
      		dMat01 <- cbind( rep(1,length(dMat)), dMat ) 
		a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) )
      		#a <- solve( qr(dMat01), selTsVoxel)
      		expectedTs <- crossprod( t(dMat01), a )
      		residualsSquared <- ( selTsVoxel - expectedTs )^2
      		ssRes <- sum( residualsSquared )
      		r2 <- 1-(ssRes/ssTot)
      		return( c( r2, a, expectedTs ) )
	}
	outVoxel <- sapply( 1:dim(tsPrediction)[1], runLinMod  )
	indexBetaZero <- outVoxel[3,] > 0
	if ( sum(indexBetaZero)>0 ) {
		indexBetaZero <- which( indexBetaZero )
		indexVarExp <- which.max( outVoxel[1,indexBetaZero] )
		return( as.numeric( c( predictionGrid[indexBetaZero[indexVarExp],], outVoxel[,indexBetaZero[indexVarExp]] ) ) ) 
	}
	if ( sum(indexBetaZero)==0 ) {
		return( rep(0, dim(predictionGrid)[2]+3+length(selTsVoxel) ) ) #no positive betas, return array of zeros, 3=r2,beta intercept, beta slope
	}
}
#system.time( aaa <- sapply( 1:length(selIdxVoxel), voxelModel ) )
library(parallel)
detectCores()
nCores <- 4
cl <- makeCluster(nCores, type='FORK')
storeTimePar <- system.time( outModel <- parSapply(cl, 1:length(selIdxVoxel), voxelModel ) )
stopCluster(cl)
print( storeTimePar )


if (modelFitCounter==1) { outModelLoop <- outModel }
if (modelFitCounter>1) { 
selectedCols <- outModel[5,] > outModelLoop[5,]
if ( sum( selectedCols ) > 0 ) {
outModelLoop[,selectedCols] <- outModel[,selectedCols] 
}
}

}


if (fineFit==1) {
}


storeAllPred <- t(outModelLoop[ seq(1,7), ])
#### extract surround size ####
print('surround size...')
if (flagSurround==1) { 
  print('get FWHM...')
  progress <- 0.05
  FWHM <- array(0, c( dim(storeAllPred)[1], 2 ) )
  xFWHM <- seq(-25,25,0.0125)
  for (k in 1:dim(storeAllPred)[1]) {
    
    parameters <- storeAllPred[k,]
    
    if ( sum(parameters!=0) ) {
      a <- dnorm( xFWHM, 0, parameters[3] )
      b <- dnorm( xFWHM, 0, parameters[4] )
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
storeAllExpectedTs[,selIdxVoxel] <- outModelLoop[ 8:dim(outModelLoop)[1], ]	
storeAllPredOut <- array( 0, c(11, dim(tsTransposedAll)[2] ) )

print('save linear step...')
polCoords <- cart2pol( storeAllPred[,c(1,2)] )
storeAllPredOut[,selIdxVoxel] <- t( cbind( storeAllPred, polCoords, FWHM ) )
fileTs <- sprintf('%s_PredixtedTs.nii.gz',outSuffix) 
fileParams <- sprintf('%s_params.nii.gz',outSuffix)

rStoreFit <- array( storeAllExpectedTs, c( dim(ts$brk)[4], dim(ts$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( '__tt_expectedTs.nii.gz', brk = rStoreFit,
            origin = ts$origin, orient = ts$orient,
            defhead = ts$NI_head)
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_expectedTs.nii.gz', orientValue[1,1], fileTs )
system( instr)
system( 'rm __tt_expectedTs.nii.gz' )

rParameters <- array( storeAllPredOut, c( dim(storeAllPredOut)[1], dim(ts$brk)[1:3] ) )
rParameters <- aperm( rParameters, c(2,3,4,1) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', epiFile )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( '__tt_parameters.nii.gz', brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_parameters.nii.gz', orientValue[1,1], fileParams )
system( instr)
system( 'rm __tt_parameters.nii.gz' )

labels <- c('x-Pos','y-Pos','sigmaPos','sigmaNeg','varExp','intercept','slope','theta','radius','fwhmCenter','surroundSize')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}

