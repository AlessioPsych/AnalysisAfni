args <- commandArgs(T)
print( args )

#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#args <- c('maxVarEye.nii.gz','meanTs_bars_res.nii','outTest_delme.nii.gz','/packages/afni/17.0.13','1')

## prepare the stimuli with the portrait as well

#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_HNA10_FEF_19102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti')
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')


#args <- c('meanTs_bars_res_mask.nii', 'meanTs_bars_res.nii', 'output_test_delme_bars','/home/alessiof/abin', '1', '1','1','0.166','3')
#args <- c('maxVarEye.nii.gz','outputTs','5','0.166','eye_fine_oblique_eyeFixBorder_noSurr_params.nii.gz','eye_fine_oblique_eyeFixBorder_noSurr_detrendedData.nii.gz',3)

mainDir <- getwd()
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

maskFile <- args[1]
outSuffix <- args[2]
fineFit <- as.numeric( args[3] )
samplingTime <- as.numeric( as.numeric(args[4]) )
paramsFile <- args[5]
newTsTest <- args[6]
stimType <- as.numeric( args[7] )

# get orientation
instr <- sprintf('3dinfo -orient %s > __tttt.1D', maskFile)
system( instr )
orientValue <- read.table('__tttt.1D')
system('rm __tttt.1D')
# get TR
instr <- sprintf('3dinfo -tr %s > __tttt.1D', maskFile)
system( instr )
trValue <- read.table('__tttt.1D')
system('rm __tttt.1D')


# load the mask
# get orientation
print('load mask...')
maskFileData <- read.AFNI( maskFile )
maskVolume <- maskFileData$brk
maskMatrix <- array( maskVolume, c( prod( dim( maskVolume )[1:3] ), dim(maskVolume)[4] ) )
maskIdx <- which(maskMatrix==1)

# already fitted model
print('load parameters...')
paramsPredicted <- read.AFNI( paramsFile )
paramsArrayPredicted <- array( paramsPredicted$brk, c( prod( dim( paramsPredicted$brk )[1:3] ), dim(paramsPredicted$brk)[4] ) )

# new ts to test
print('load ts to test predictions...')
tsToTest <- read.AFNI( newTsTest )
tsToTestArray <- array( tsToTest$brk, c( prod( dim( tsToTest$brk )[1:3] ), dim(tsToTest$brk)[4] ) )

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
  arrayStim <- scan( 'eyeFixStim_border.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) )
}
if (stimType==4) { 
  arrayStim <- scan( 'prfStim.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1860,135) ), c( 3, 1, 2 ) ) #
}
if (stimType==5) { 
  arrayStim <- scan( 'eyeFixStim_border_disappear.txt' )
  setwd(mainDir)
  stimMat <- aperm( array( arrayStim, c(240,1510,135) ), c( 3, 1, 2 ) ) #
}

print('show new stimulus...')
stimMatFlip <- aperm( stimMat[ dim(stimMat)[1]:1,, ], c(2,1,3) )
x11( width=3, height=3 )
for ( snap in 1:dim(stimMat)[3] ) {
  image( stimMatFlip[,,snap], axes=FALSE ); par(new=TRUE); Sys.sleep(0.01) 
}
stimSeq <- stimMatFlip
x <- seq(-10,10,length.out = dim(stimSeq)[1] )
y <- seq(-5.5,5.5,length.out = dim(stimSeq)[2] )
outMesh <- meshgrid(y, x)
outMesh$X <- scaleData( outMesh$X, 1, 0 )
outMesh$Y <- scaleData( outMesh$Y, 1, 0 )
rm( stimMat )
rm( stimMatFlip )

# define new stimuli and timing
print( sprintf( 'generate and test predictions per voxels'  ) )
stimSeqMat <- array( stimSeq, c( length(x)*length(y), dim(stimSeq)[3] ) )
timeStimuli <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(stimSeq)[3] )
mriTime <- seq( 0, dim(stimSeq)[3]*samplingTime, length.out = dim(tsToTest$brk)[4] )

generateAndTest <- function(voxelIdx) { #voxelIdx <- 10
  
  singleVoxelIdx <- maskIdx[ voxelIdx ]
  prfPar <- paramsArrayPredicted[singleVoxelIdx,]
  
  hrf <- canonicalHRF( seq(0,30,samplingTime), param=list(a1=prfPar[5], a2=prfPar[6], b1=0.9, b2=0.9, c=0.35), verbose=FALSE )
  
  a <- dnorm( x, prfPar[1], prfPar[3] )
  b <- dnorm( y, prfPar[2], prfPar[3] )
  imgCenter <- scaleData( tcrossprod(a,b), 1, 0 )
  
  a <- dnorm( x, prfPar[1], prfPar[4] )
  b <- dnorm( y, prfPar[2], prfPar[4] )
  imgSurround <- scaleData( tcrossprod(a,b), 1, 0)*prfPar[7]
  
  if (fineFit==-1 | fineFit==-2) { #no gain
    r <- imgCenter - imgSurround
  }
  if (fineFit==0 | fineFit==1 ) { #plane gain
    imgGain <- outMesh$X*prfPar[8] + outMesh$Y*prfPar[9]
    r <- (imgCenter - imgSurround) + imgGain
  }
  if (fineFit==2 | fineFit==3) { #gaussian gain
    a <- dnorm( x, prfPar[8], prfPar[10] )
    b <- dnorm( y, prfPar[9], prfPar[10] )
    imgGain <- scaleData( tcrossprod(a,b), 1, 0 )
    r <- (imgCenter - imgSurround) + imgGain
  }
  if (fineFit==4 | fineFit==5) {
    meshOut = meshgrid(y, x);
    X1 <- meshOut$X
    X2 <- meshOut$Y
    sigma1 <- prfPar[8];
    sigma2 <- prfPar[3];
    Theta <- prfPar[10];
    mu01 <- prfPar[2];
    mu02 <- prfPar[1];
    
    a <- ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
    b <- -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
    c <- ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));
    
    A <- 1;
    rPos <- A*exp(-(a*(X1 - mu01)^2 + 2*b*(X1 - mu01)*(X2 - mu02) + c*(X2 - mu02)^2))
    
    sigma1 <- prfPar[9];
    sigma2 <- prfPar[4];
    Theta <- prfPar[10];
    mu01 <- prfPar[2];
    mu02 <- prfPar[1];
    
    a <- ((cosd(Theta)^2) / (2*sigma1^2)) + ((sind(Theta)^2) / (2*sigma2^2));
    b <- -((sind(2*Theta)) / (4*sigma1^2)) + ((sind(2*Theta)) / (4*sigma2^2));
    c <- ((sind(Theta)^2) / (2*sigma1^2)) + ((cosd(Theta)^2) / (2*sigma2^2));
    
    A <- 1;
    rNeg <- A*exp(-(a*(X1 - mu01)^2 + 2*b*(X1 - mu01)*(X2 - mu02) + c*(X2 - mu02)^2))
    r <- rPos-rNeg*prfPar[7]
  }
  
  r[ r<quantile(r,0.1) ] <- 0
  r <- scaleData(r,1,0)
  rMat <- array(r)
  
  predictionLoop <- as.numeric( crossprod( stimSeqMat,rMat ) ) #### this is the slow step, this is the reason for the parallel computing ####
  pConv <- conv( predictionLoop, hrf )
  pConvTrim <- pConv[ 1 : dim(stimSeq)[3] ]
  tsPredictionMriInterp <- interp1( x=timeStimuli, y=pConvTrim, xi=mriTime, method=c('nearest') )
  returnPrediction <- round( scaleData( tsPredictionMriInterp, 1, 0 ), 5 ) * prfPar[13] + prfPar[12]
  #par(mfrow=c(2,2))
  #plot( r[,70], type='l' )
  #image( r, col=gray.colors(500) )
  #plot( returnPrediction, type='l' )
  #image( imgGain, col=gray.colors(500)  )
  
  selTsVoxel <- tsToTestArray[singleVoxelIdx,]
  selTsVoxelMean <- mean( selTsVoxel )
  ssTot <- sum((selTsVoxel-selTsVoxelMean)^2 ) #voxels total sum of squares
  dMat <- cbind( returnPrediction )      
  dMat01 <- cbind( rep(1,length(dMat)), dMat ) #column of ones and column of predictor
  a <- solve( crossprod(dMat01,dMat01), crossprod(dMat01,selTsVoxel) ) #beta coefficients (2: intercept and slope)
  expectedTs <- crossprod( t(dMat01), a ) #expected ts
  residualsSquared <- ( selTsVoxel - expectedTs )^2
  ssRes <- sum( residualsSquared )
  r2 <- 1-(ssRes/ssTot) #r squares
  
  #to test r2must be equal to prfPar[11], new r2==oldR2: c(r2, prfPar[11])
  #the coefficients (intercept and beta must be close to 0 and 1,respectively): c(a[1],a[2])
  
  return( c(r2, a, returnPrediction) )

}
#bbb <- sapply( 100:500, generateAndTest )
detectCores()
nCores <- 6
cl <- makeCluster(nCores, type='FORK')
runIndex <- 1:length(maskIdx)
storeTimePar <- system.time( outPredictions <- parSapply(cl, runIndex, generateAndTest ) )
stopCluster(cl)
print( storeTimePar )

storeAllExpectedTs <- array( 0, dim( t(tsToTestArray) ) )
storeAllExpectedTs[,maskIdx[runIndex]] <- outPredictions[ 4:dim(outPredictions)[1], ]	
storeAllPredOut <- array( 0, c(3, dim( t(tsToTestArray) )[2] ) )
storeAllPredOut[,maskIdx[runIndex]] <- outPredictions[ 1:3, ]

print('save linear step...')
fileTs <- sprintf('%s_PredixtedTs.nii.gz',outSuffix) 
fileParams <- sprintf('%s_params.nii.gz',outSuffix)

rStoreFit <- array( storeAllExpectedTs, c( dim(tsToTest$brk)[4], dim(tsToTest$brk)[1:3] ) )
rStoreFit <- aperm( rStoreFit, c(2,3,4,1) )
write.AFNI( '__tt_expectedTs.nii.gz', brk = rStoreFit,
            origin = tsToTest$origin, orient = tsToTest$orient,
            defhead = tsToTest$NI_head)
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_expectedTs.nii.gz', orientValue[1,1], fileTs )
system( instr)
system( 'rm __tt_expectedTs.nii.gz' )

rParameters <- array( storeAllPredOut, c( dim(storeAllPredOut)[1], dim(tsToTest$brk)[1:3] ) )
rParameters <- aperm( rParameters, c(2,3,4,1) )
instr <- sprintf( '3dcalc -a %s -expr \u0027step(a)\u0027 -prefix __tt.nii.gz -datum float', maskFile )
system( instr )
appDataset <- read.AFNI('__tt.nii.gz')
write.AFNI( '__tt_parameters.nii.gz', brk = rParameters,
            origin = appDataset$origin, orient = appDataset$orient,
            defhead = appDataset$NI_head )
system('rm __tt.nii.gz')
instr <- sprintf( '3dresample -orient %s -prefix %s -inset __tt_parameters.nii.gz', orientValue[1,1], fileParams )
system( instr)
system( 'rm __tt_parameters.nii.gz' )

labels <- c('varExp','intercept','slope')
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], fileParams )
  print( instr )
  system( instr )
}


