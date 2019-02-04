args <- commandArgs(T)
print( args )

#setwd('/analyse/Project0226/GN18NE278_GVW19_FEF_05102018_nifti/ANATOMY')
#args <- c('leftCalcarineSelected.nii.gz', 'profilesOut_profiles.nii.gz', 'output_test_profiles01')

# get input parameters and libraries
print('get input parameters and libraries...')
mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
inputRoi <- args[1]
inputProfiles <- args[2]
outputFilename <- args[3]
#library(pracma)
library(parallel)

# get input data
print('get input data...')
inputRoiFile <- read.AFNI( inputRoi )
inputProfilesFile <- read.AFNI( inputProfiles )
inputRoiVolume <- inputRoiFile$brk[,,,1]
inputProfilesVolume <- inputProfilesFile$brk

# arrange profiles in the roi from volume to matrix form
print('arrange profiles in the roi from volume to matrix form...')
roiIdx <- which( inputRoiVolume==1 )
profilesMatRaw <- array(0, c( length(roiIdx), dim(inputProfilesVolume)[4] ) )
for (levelCounter in 1:dim(inputProfilesVolume)[4]) {
  tempVolume <- inputProfilesVolume[,,,levelCounter]
  profilesMatRaw[,levelCounter] <- tempVolume[ roiIdx ]
}
profilesMat <- profilesMatRaw

#corticalDepthRaw <- seq( 0, 1, length.out = dim(profilesMatRaw)[2] )
#corticalDepth <- seq( 0, 1, length.out = 20 )
#profilesMat <- array(0, c( dim(profilesMatRaw)[1], length( corticalDepth ) ) )
#for (nProfiles in 1:dim(profilesMatRaw)[1]) {
#  profilesMat[nProfiles,] <- interp1( corticalDepthRaw, profilesMatRaw[nProfiles,], corticalDepth, method=c('linear') )
#}

xPos <- seq( -0.25, 1.25, length.out = 12 )
sdPar <- seq( 0.05, 1.25, length.out = 12 )
testMat <- expand.grid( xPos, sdPar )

runModel <- function(testIndex,runIndex) { #run all models on a single profile
  a <- testMat[testIndex,1] #xPosition
  b <- testMat[testIndex,2] #standard deviation
  inData <- profilesMat[runIndex,]
  gaussArray <- dnorm( seq(0,1,length.out = dim(profilesMat)[2]), a, b )
  gaussArray <- scaleData( gaussArray, 1, 0)
  linArray <- seq(0,1,length.out = dim(profilesMat)[2])
  modLLBase <- lm( inData ~ poly(linArray,2) )
  modLL <- lm( inData ~ gaussArray + poly(linArray,2) )
  comp <- anova(modLLBase,modLL)
  targetProfile <- modLL$fitted.values
  return( c( modLL$fitted.values, a, b, comp$F[2], comp$`Pr(>F)`[2], modLL$coefficients, summary(modLL)$r.squared ) )
}

runAllProfiles <- function(nProfile) { #run all profiles
  singleProfileWholeResults <- sapply( 1:dim(testMat)[1], runModel, runIndex=nProfile )
  positiveGaussParLogical <- singleProfileWholeResults[ dim(singleProfileWholeResults)[1]-3, ] > 0 # are there gaussian parameters > zero
  if (sum( positiveGaussParLogical )>0) { #if there are gaussian parameters larger than zero:
    positiveGaussPar <- which( positiveGaussParLogical )
    winningModelIdx <- which.max( singleProfileWholeResults[ dim(singleProfileWholeResults)[1], positiveGaussPar ] )
    singleProfileResults <- singleProfileWholeResults[, positiveGaussPar[ winningModelIdx ] ]
    return( singleProfileResults )
  }
  if (sum( positiveGaussParLogical )==0) { #otherwise
    return( rep(0, dim(profilesMat)[2]+8 ) ) #parameters stored: the fitted profile + gaussianCenter, sd, f comparison, p comparison, intercept, slope gauss, slope lin, slope quad, r2 
  }
}
#system.time( allProfilesResults <- sapply( 1:15, runAllProfiles ) ); dim(allProfilesResults)
print( sprintf( 'parallel fitting...' )  )
no_cores <- 6
cl <- makeCluster(no_cores, type='FORK')
storeTime <- system.time( outParallel <- parSapply( cl, 1:dim(profilesMat)[1], runAllProfiles ) )
print( storeTime )
stopCluster(cl)

### to debug and test
#nProfile <- 99
#singleProfileWholeResults <- sapply( 1:dim(testMat)[1], runModel, runIndex=nProfile )
#positiveGaussPar <- which( singleProfileWholeResults[ dim(singleProfileWholeResults)[1]-3, ] > 0 )
#winningModelIdx <- which.max( singleProfileWholeResults[ dim(singleProfileWholeResults)[1], positiveGaussPar ] )
#singleProfileResults <- singleProfileWholeResults[, positiveGaussPar[ winningModelIdx ] ]
#singleProfileResults[ length(singleProfileResults) ]
#fitIndex <- nProfile
#predX <- seq( 0,1,length.out = dim(profilesMat)[2] )
#plot( profilesMat[fitIndex,]~predX, bty='n', type='l' )
#lines( predX, singleProfileResults[1:dim(profilesMat)[2]], type='l', col='red' )
#abline( lm(profilesMat[fitIndex,]~predX), col='blue')
#singleProfileResults[ c(9:12,length(singleProfileResults)) ]

predictedProfiles <- outParallel[ 1:dim(profilesMat)[2], ]
parametersProfiles <- outParallel[ (dim(profilesMat)[2]+1):dim(outParallel)[1],  ]
# saving predicted profiles
print('saving predicted profiles...')
storePredictedProfiles <- array( 0, dim(inputProfilesVolume) )
for (levelCounter in 1:dim(inputProfilesVolume)[4]) {
  tempVolume <- storePredictedProfiles[,,,levelCounter]
  tempVolume[ roiIdx ] <- predictedProfiles[levelCounter,]
  storePredictedProfiles[,,,levelCounter] <- tempVolume
}
predictedProfilesName <- sprintf( '%s_predictedProfiles.nii.gz' , outputFilename )
write.AFNI( predictedProfilesName, brk = storePredictedProfiles,
            origin = inputProfilesFile$origin, orient = inputProfilesFile$orient,
            defhead = inputProfilesFile$NI_head)

# saving profiles parameters
print('saving profiles parameters...')
storeParameterProfiles <- array( 0, c( dim(inputProfilesVolume)[1:3], dim(parametersProfiles)[1] ) )
for (levelCounter in 1:dim(parametersProfiles)[1]) {
  tempVolume <- storeParameterProfiles[,,,levelCounter]
  tempVolume[ roiIdx ] <- parametersProfiles[levelCounter,]
  storeParameterProfiles[,,,levelCounter] <- tempVolume
}
predictedProfilesName <- sprintf( '%s_parameterProfiles.nii.gz' , outputFilename )
write.AFNI( predictedProfilesName, brk = storeParameterProfiles,
            origin = inputProfilesFile$origin, orient = inputProfilesFile$orient,
            defhead = inputProfilesFile$NI_head)
labels <- c('gaussianCenter', 'sd', 'f_comparison', 'p_comparison', 'intercept', 'slope_gauss', 'slope_lin', 'slope_quad', 'r2' )
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], predictedProfilesName )
  print( instr )
  system( instr )
}
















































# 
# 
# 
# LL <- function(a,b,muLL,sigmaLL,inData) {
#   gaussArray <- dnorm( seq(0,1,length.out = dim(profilesMat)[2]), a, b )
#   gaussArray <- scaleData( gaussArray, 1, 0)
#   linArray <- seq(0,1,length.out = dim(profilesMat)[2])
#   modLL <- lm( inData ~ gaussArray + poly(linArray,2) )
#   targetProfile <- modLL$fitted.values
#   R <- inData - targetProfile
#   R <- suppressWarnings( dnorm(R, muLL, sigmaLL, log=TRUE) )
#   return( -sum(R) )
# }
# plotPredProfile <- function( a, b , inData) {
#   gaussArray <- dnorm( seq(0,1,length.out = dim(profilesMat)[2]), a, b )
#   gaussArray <- scaleData( gaussArray, 1, 0)
#   linArray <- seq(0,1,length.out = dim(profilesMat)[2])
#   modLL <- lm( inData ~ gaussArray + linArray )
#   summaryStat <- summary( modLL )
#   targetProfile <- modLL$fitted.values
#   return( c( targetProfile, coefficients(modLL), summaryStat$r.squared, a, b ) )
# }
# arrayFit <- function( idx, inputProfiles) {
#   inputData <- array( inputProfiles[idx,] )
#   fitLL <- tryCatch( mle( LL, start = list( a = 0.5, b = 0.1 ),
#                           fixed = list(muLL = 0, sigmaLL = 1, inData=inputData),
#                           nobs = length(inputData),
#                           method='L-BFGS-B',
#                           lower=c(-1,0.005),
#                           upper=c(2,1.25),
#                           control=list(maxit=500) ),
#                      error = function(e) 'error' )
#   if (!isS4(fitLL)) { fitParLoop <- rep( 0, c( length(inputData) + 6 ) ) }
#   if (isS4(fitLL)) { fitParLoop <- plotPredProfile( coef(fitLL)[1], coef(fitLL)[2], inData=inputData ) }
#   return( as.numeric( fitParLoop ) )
# }
# 
# # to debug and test
# fitIndex <- 99
# fitSol <- arrayFit( fitIndex, profilesMat )
# predX <- seq( 0,1,length.out = dim(profilesMat)[2] )
# plot( profilesMat[fitIndex,]~predX, bty='n', type='l' )
# lines( predX, plotPredProfile( fitSol[ length(fitSol)-1 ], fitSol[ length(fitSol) ], inData=profilesMat[fitIndex,] )[1:dim(profilesMat)[2]], type='l', col='red' )
# abline( lm(profilesMat[fitIndex,]~predX), col='blue')
# fitSol
# 
# #outParallel <- sapply( 1:100, arrayFit, inputProfiles=profileMatrixRight )
# #sapply all the profiles withthe mle fit and extract the parameters
# 
# #sapply(1:1000, arrayFit, inputProfiles=profileMatrixRight )
# print( sprintf( 'processing n: %1.0f ... fitting ...', nDataset )  )
# #library( parallel )
# no_cores <- 6
# cl <- makeCluster(no_cores, type='FORK')
# storeTime <- system.time( outParallelRight <- parSapply( cl, 1:dim(profileMatrixRight )[1], arrayFit, inputProfiles=profileMatrixRight ) )
# print( storeTime )
# storeTime <- system.time( outParallelLeft <- parSapply( cl, 1:dim(profileMatrixLeft )[1], arrayFit, inputProfiles=profileMatrixLeft ) )
# print( storeTime )
# stopCluster(cl)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# LL <- function(a,b,muLL,sigmaLL,inData) {
#   gaussArray <- dnorm( seq(0,1,length.out = dim(profilesMat)[2]), a, b )
#   gaussArray <- scaleData( gaussArray, 1, 0)
#   linArray <- seq(0,1,length.out = dim(profilesMat)[2])
#   modLL <- lm( inData ~ gaussArray + linArray )
#   targetProfile <- modLL$fitted.values
#   R <- inData - targetProfile
#   R <- suppressWarnings( dnorm(R, muLL, sigmaLL, log=TRUE) )
#   return( -sum(R) )
# }
# plotPredProfile <- function( a, b , inData) {
#   gaussArray <- dnorm( seq(0,1,length.out = dim(profilesMat)[2]), a, b )
#   gaussArray <- scaleData( gaussArray, 1, 0)
#   linArray <- seq(0,1,length.out = dim(profilesMat)[2])
#   modLL <- lm( inData ~ gaussArray + linArray )
#   summaryStat <- summary( modLL )
#   targetProfile <- modLL$fitted.values
#   return( c( targetProfile, coefficients(modLL), summaryStat$r.squared, a, b ) )
# }
# arrayFit <- function( idx, inputProfiles) {
#   inputData <- array( inputProfiles[idx,] )
#   fitLL <- tryCatch( mle( LL, start = list( a = 0.5, b = 0.1 ),
#                           fixed = list(muLL = 0, sigmaLL = 1, inData=inputData),
#                           nobs = length(inputData),
#                           method='L-BFGS-B',
#                           lower=c(-1,0.005),
#                           upper=c(2,1.25),
#                           control=list(maxit=500) ),
#                      error = function(e) 'error' )
#   if (!isS4(fitLL)) { fitParLoop <- rep( 0, c( length(inputData) + 6 ) ) }
#   if (isS4(fitLL)) { fitParLoop <- plotPredProfile( coef(fitLL)[1], coef(fitLL)[2], inData=inputData ) }
#   return( as.numeric( fitParLoop ) )
# }
# 
# # to debug and test
# fitIndex <- 99
# fitSol <- arrayFit( fitIndex, profilesMat )
# predX <- seq( 0,1,length.out = dim(profilesMat)[2] )
# plot( profilesMat[fitIndex,]~predX, bty='n', type='l' )
# lines( predX, plotPredProfile( fitSol[ length(fitSol)-1 ], fitSol[ length(fitSol) ], inData=profilesMat[fitIndex,] )[1:dim(profilesMat)[2]], type='l', col='red' )
# abline( lm(profilesMat[fitIndex,]~predX), col='blue')
# fitSol
# 
# #outParallel <- sapply( 1:100, arrayFit, inputProfiles=profileMatrixRight )
# #sapply all the profiles withthe mle fit and extract the parameters
# 
# #sapply(1:1000, arrayFit, inputProfiles=profileMatrixRight )
# print( sprintf( 'processing n: %1.0f ... fitting ...', nDataset )  )
# #library( parallel )
# no_cores <- 6
# cl <- makeCluster(no_cores, type='FORK')
# storeTime <- system.time( outParallelRight <- parSapply( cl, 1:dim(profileMatrixRight )[1], arrayFit, inputProfiles=profileMatrixRight ) )
# print( storeTime )
# storeTime <- system.time( outParallelLeft <- parSapply( cl, 1:dim(profileMatrixLeft )[1], arrayFit, inputProfiles=profileMatrixLeft ) )
# print( storeTime )
# stopCluster(cl)
# 
# 
# 
