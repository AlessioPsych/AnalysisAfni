args <- commandArgs(T)
print( args )
mainDir <- getwd()

#rm(list=ls())
#mainDir <- '/home/fracasso/data/ArjanProject/Arjan_fMRI_scripts/results'
#setwd( mainDir )
#args <- c('roiBoundary00.1D.roi','0.1','boundary00','mapOnSurface_smoothed_8/','15','/packages/afni/16.1.13','/home/fracasso/analysisAfni')

#mainDir <- '/analyse/Project0226/stripesData/P1/analysis_figures' 
#setwd( mainDir )
#args <- c( 'V2DorsalLeft_filled_interp_surfaces_folder/roiOut_0_filled.1D.roi', '0', 'boundary00', 'phaseSpecCohMax_interp_surfaces_folder/boundary00_sm_phaseSpecCohMax_surf.1D.dset','7','6','v2DorsalOut','1','/home/alessiof/abin','/home/alessiof/Programs/AnalysisAfni')

#rm(list=ls())
#mainDir <- '/analyse/Project0226/stripesData/P2/analysis_figures'
#setwd( mainDir )
#args <- c('test02_interp_surfaces_folder/boundary03_roiCleanCluster.1D.roi',#'V2DorsalLeft_filled_interp_surfaces_folder/roiOut_1_filled.1D.roi', test02_interp_surfaces_folder/boundary01_roiCleanCluster.1D.roi, test02.1D.roi
#'49200', 
#'boundary03',
#'anatomy_interp_surfaces_folder/boundary03_sm_anatomy_surf.1D.dset', # phaseSpecCohMax anatomy
#'7', 
#'1', 
#'test02',
#'1',
#'/home/alessiof/abin',
#'/home/alessiof/Programs/AnalysisAfni')

source( sprintf('%s/AFNIio.R', args[9]) )
source( sprintf('%s/generalPurpose/load1DData.R', args[10]) )

# copy roi in current directory
instr <- sprintf( 'cp %s ttt_tempRoi.1D.roi', args[1] )
system( instr )

# which column index in the surface dataset are we looking for?
indexSurfaceColumn <- as.numeric( args[5] )

# read data in
setwd(mainDir)
roiData <- read.table( 'ttt_tempRoi.1D.roi' )
##################################
################################## for some reason the roi indexes need to be sorted for the code to work, not clear why, but it works
##################################
roiData[,1] <- sort( roiData[,1] )
##################################
##################################
##################################
mapData <- read.table( args[4] )

# get surface intensity
roiLines <- which( is.element( mapData[,1], roiData[,1] ) )
#roiIntensity <- mapData[ roiData[,1] ,indexSurfaceColumn]
roiIntensity <- mapData[ roiLines ,indexSurfaceColumn]

# select first portion of the roi based on the criteria
selNodesIdxNegative <- which( roiIntensity < as.numeric( args[2] ) )
selNodesIdxNegative <- selNodesIdxNegative[ seq(1,length(selNodesIdxNegative),by=as.numeric(args[6])) ] #downsample
selNodesNeg <- roiData[selNodesIdxNegative,1]
write.table( cbind( selNodesNeg, rep(1, length(selNodesNeg) ) ), file = 'ttt_boundNeg2.1D.roi', row.names = FALSE, col.names = FALSE  ) 

# select second portion of the roi based on the criteria
selNodesIdxPositive <- which( roiIntensity > as.numeric( args[2] ) )
selNodesIdxPositive <- selNodesIdxPositive[ seq(1,length(selNodesIdxPositive),by=as.numeric(args[6])) ] #downsample
selNodesPos <- roiData[selNodesIdxPositive,1]
write.table( cbind( selNodesPos, rep(1, length(selNodesPos) ) ), file = 'ttt_boundPos2.1D.roi', row.names = FALSE, col.names = FALSE  ) 

# prepare file to compute distances betwee roi subparts
nodePairs01 <- rep( selNodesPos, rep( length(selNodesNeg), length(selNodesPos) ) )
nodePairs02 <- rep( selNodesNeg, length(selNodesPos) )
write.table( cbind( nodePairs01, nodePairs02 ), file='ttt_datasetStart.1D', col.names=FALSE, row.names=FALSE  )

print('compute distance from starting nodes to positive roi nodes, it might take a while (see downsampling factor, the higher, the faster)...')

instr <- sprintf( 'SurfDist -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -input ttt_datasetStart.1D > ttt_exampleStart.1D', args[3], args[3] )
system.time( system( instr ) )

print('compute smallest distance from each node and the start line...')

fromStart <- read.table('ttt_exampleStart.1D')

uniqueRoiPoints <- unique( fromStart[,2] )
uniqueRoiPointsStart <- uniqueRoiPoints
mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
for (k in 1:length( uniqueRoiPoints) ) {
  indexSel <- which( fromStart[,2]==uniqueRoiPoints[k] )
  selDF <- fromStart[indexSel,]
  indexMin <- which( selDF[,3]==min( selDF[,3] ) )[1]
  mapOutStart[k,] <- as.numeric( selDF[indexMin,c(2,3)] )
}

#mat01 <- cbind( selNodesNeg, rep(1, length(selNodesNeg) ), mapOutStart )
mat01 <- cbind( mapOutStart )
write.table( mat01, file = 'ttt_boundDist01.dset', row.names = FALSE, col.names = FALSE  ) 

uniqueRoiPoints <- unique( fromStart[,1] )
uniqueRoiPointsEnd <- uniqueRoiPoints
mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
for (k in 1:length( uniqueRoiPoints) ) {
  indexSel <- which( fromStart[,1]==uniqueRoiPoints[k] )
  selDF <- fromStart[indexSel,]
  indexMin <- which( selDF[,3]==min( selDF[,3] ) )[1]
  mapOutStart[k,] <- as.numeric( selDF[indexMin,c(1,3)] )
}

#mat02 <- cbind( selNodesPos, rep(1, length(selNodesPos) ), mapOutStart )
mat02 <- cbind( mapOutStart )
write.table( mat02, file = 'ttt_boundDist02.dset', row.names = FALSE, col.names = FALSE  ) 

mat01Corr <- mat01
mat01Corr[,2] <- mat01Corr[,2]*-1
matOut <- rbind(mat01Corr, mat02)

filenameOutDownsampled <- sprintf('%s_%s_downsampled.dset', args[7], args[3])

write.table( matOut, file = filenameOutDownsampled, row.names = FALSE, col.names = FALSE  ) 

##################################################
################################################## uncomment to load the distance map just computed, extract the data from the map and plot it:
##################################################
#distanceIn <- read.table( filenameOutDownsampled )
#mapIntensity <- mapData[ distanceIn[,1]+1, 7 ]
#plot( mapIntensity ~ distanceIn[,2] )
##################################################
##################################################
##################################################

#dataTemp <- cbind( c( mapData[ mat02[,1]+1 ,indexSurfaceColumn], mapData[ mat01[,1]+1 ,indexSurfaceColumn] ), c( mat01[,4], mat02[,4] ) )
#dataTemp <- cbind( c( mapData[ mat02[,1]+1 ,indexSurfaceColumn], mapData[ mat01[,1]+1 ,indexSurfaceColumn] ), c( mat01[,4], mat02[,4] ) )

#roiLines <- which( is.element( mapData[,1], matOut[,1] ) )
#roiIntensity <- mapData[ roiLines ,indexSurfaceColumn]

#plot( dataTemp[,1] ~ dataTemp[,2] )
#x <- rnorm( 100,1,1)
#xPos <- x[ x>0 ]
#median(x); median(xPos)
# problem 1 - possible problem with assignment between intensity and distance
# problem 2 - possible problem with downsampling
# save and load data to test...

## add smoothing over the surface

emptyMap <- array(0, c( dim(mapData)[1], 2) )
emptyMap[ matOut[,1]+1, 1 ] <- 1
emptyMap[ matOut[,1]+1, 2 ] <- matOut[,2]

print('writing distances computed on the downsampled nodes...')

write.table( matOut[,1], file='ttt_outMap_indices.1D', row.names=FALSE, col.names=FALSE )

print( sprintf( 'interpolate over the surface (kernel diameter = %smm )...', args[8] ) )

system( 'mkdir surfInterp' )
instr <- sprintf( 'ROIgrow -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -roi_labels PER_NODE -roi_nodes ttt_outMap_indices.1D -lim %f -prefix roiMap001Fill', args[3], args[3], as.numeric(args[8])/2 )
system( instr )
system('mv roiMap001Fill* surfInterp/')
setwd('surfInterp/')
files <- dir(pattern='roiMap001Fill*')

for ( k in 1:length( files ) ) {
  fileContent <- read.table( files[k], as.is=TRUE )
  distValues <- emptyMap[ fileContent[,1], 2 ]
  meanValue <- median( distValues[ abs( distValues ) > 0 ] ) # mean distance value != 0 inside the kernel around the selected node
  meanValues <- rep( meanValue, length( distValues ) ) # rep mean distance value != 0 inside the kernel around the selected node for the number of nodes where distValue != 0
  if (k==1) {
    storeDataMean <- meanValues
    storeDataNodes <- fileContent[,1]
  }
  if (k>1) {
    storeDataMean <- c( storeDataMean, meanValues )
    storeDataNodes <- c( storeDataNodes, fileContent[,1] )
  }
}

# for all the unique nodes in storeDataNodes, extract the mean distance and report the average
uniqueNodes <- unique( storeDataNodes ) 
storeNodes <- rep(0,length(uniqueNodes))
for (k in 1:length( uniqueNodes ) ) {
  indexNodes <- which( uniqueNodes[k]==storeDataNodes )
  meanDist <- mean( storeDataMean[ indexNodes ], na.rm=TRUE )
  storeNodes[k] <- meanDist
} 

# clean up in case of NAs
if ( sum( is.na(storeNodes) )>0 ) {
  indxNan <- which( is.na(storeNodes) )
  storeNodes[indxNan] <- 0
}

setwd(mainDir)

outMap <- cbind( uniqueNodes, storeNodes )
filenameOutInterpolated <- sprintf('%s_%s_interpolated.dset', args[7], args[3])

print( sprintf( 'write interpolated distance map (filename = %s)...', filenameOutInterpolated) )

write.table( outMap, file=filenameOutInterpolated, row.names=FALSE, col.names=FALSE )

system('rm -R surfInterp/')
system('rm ttt_*')


