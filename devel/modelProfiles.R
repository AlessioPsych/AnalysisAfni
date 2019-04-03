args <- commandArgs(T)
print( args )

#setwd('/data1/projects/myelin/myelinData/hemiBackup/controldata/testModel')
#args <- c('del_left_erode02.nii.gz', 'LIZ_HC_01_pdCorrectRegularT1_noBlur_strippedcortexMask_calcarine_car_box_seg_anatomy_intensity_orig.nii.gz', 'output_test_profiles01')

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
modelType <- as.numeric( args[4] )
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

xPos <- seq( -0.25, 1.25, length.out = 16 )
sdPar <- seq( 0.05, 1.25, length.out = 14 )
testMat <- expand.grid( xPos, sdPar )

runModel <- function(testIndex,runIndex) { #run all models on a single profile
  a <- testMat[testIndex,1] #xPosition
  b <- testMat[testIndex,2] #standard deviation
  inData <- profilesMat[runIndex,]
  gaussArray <- dnorm( seq(0,1,length.out = dim(profilesMat)[2]), a, b )
  gaussArray <- scaleData( gaussArray, 1, 0)
  linArray <- seq(0,1,length.out = dim(profilesMat)[2])
  #modLLBase <- lm( inData ~ poly(linArray,2) )
  
  if (modelType==1) {
    modLL <- lm( inData ~ gaussArray + poly(linArray,2) )
    sumLL <- summary( modLL )
    tGauss <- sumLL$coefficients[2,3]
    pGauss <- sumLL$coefficients[2,4]
    #comp <- anova(modLLBase,modLL)
    targetProfile <- modLL$fitted.values
    return( c( modLL$fitted.values, a, b, tGauss, pGauss, modLL$coefficients, summary(modLL)$r.squared ) )
  }
  if (modelType==2) {
    modLL <- lm( inData ~ gaussArray + poly(linArray,1) )
    sumLL <- summary( modLL )
    tGauss <- sumLL$coefficients[2,3]
    pGauss <- sumLL$coefficients[2,4]
    #comp <- anova(modLLBase,modLL)
    targetProfile <- modLL$fitted.values
    return( c( modLL$fitted.values, a, b, tGauss, pGauss, modLL$coefficients, 999, summary(modLL)$r.squared ) )
  }
  
  #return( c( modLL$fitted.values, a, b, comp$F[2], comp$`Pr(>F)`[2], modLL$coefficients, summary(modLL)$r.squared ) )
  
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
    return( rep(0, dim(profilesMat)[2]+9 ) ) #parameters stored: the fitted profile + gaussianCenter, sd, f comparison, p comparison, intercept, slope gauss, slope lin, slope quad, r2 
  }
}
#system.time( allProfilesResults <- sapply( 1:15, runAllProfiles ) ); dim(allProfilesResults)
print( sprintf( 'parallel fitting...' )  )
no_cores <- 6
cl <- makeCluster(no_cores, type='FORK')
#storeTime <- system.time( outParallel <- parSapply( cl, 1:500, runAllProfiles ) )
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
labels <- c('gaussianCenter', 'sd', 't_gauss', 'p_gauss', 'intercept', 'slope_gauss', 'slope_lin', 'slope_quad', 'r2' )
for (k in 1:length(labels)) {
  instr <- sprintf('3drefit -sublabel %1.0f %s %s', round(k-1,0), labels[k], predictedProfilesName )
  print( instr )
  system( instr )
}

