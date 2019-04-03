args <- commandArgs(T)
print( args )

#setwd('/data1/projects/myelin/myelinData/hemiBackup/glaucomadata/test/V6710')
#args <- c('rightCalcarine_roi_test_axis.dset_clust.nii.gz', 'rightMaps01','100')

# get input parameters and libraries
print('get input parameters and libraries...')
mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library(pracma)
inputRoi <- args[1]
outputFilename <- args[2]
scalingFactor <- as.numeric( args[3] )

# get input data
print('get input data...')
inputRoiFile <- read.AFNI( inputRoi )
eccVolume <- inputRoiFile$brk[,,,1]
polVolume <- inputRoiFile$brk[,,,2]
eccValues <- eccVolume[ abs(eccVolume)>0.0001 ]
polValues <- polVolume[ abs(eccVolume)>0.0001 ]
eccValuesScaled <- ( eccValues + abs( min(eccValues) ) )
polValuesScaled <- scaleData( polValues, pi, -pi )

i <- 1i
a <- 0.75
b <- 90
w <- log( (  ( eccValuesScaled  *  exp(i*(polValuesScaled))  ) + a) / (( eccValuesScaled  *  exp(i*polValuesScaled)  ) + b ));
#w <- log( (  ( polValuesScaled  *  exp(i*(eccValuesScaled))  ) + a) ) ;

logAngle <- angle(w);
logRad <- scaleData( abs(w), 100, 0 );

eccMap <- array(0, dim( inputRoiFile$brk[,,,1] ) )
polMap <- array(0, dim( inputRoiFile$brk[,,,1] ) )
eccMap[ abs(eccVolume)>0.0001 ] <- logRad
polMap[ abs(eccVolume)>0.0001 ] <- logAngle

# saving maps
print('saving maps...')
storeMaps <- array( 0, dim( inputRoiFile$brk ) )
storeMaps[,,,1] <- eccMap
storeMaps[,,,2] <- polMap
mapsName <- sprintf( '%s.nii.gz' , outputFilename )
write.AFNI( mapsName, brk = storeMaps,
            origin = inputRoiFile$origin, orient = inputRoiFile$orient,
            defhead = inputRoiFile$NI_head)

