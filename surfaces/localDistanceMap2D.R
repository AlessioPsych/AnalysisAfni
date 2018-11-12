#args <- commandArgs(T)
#print( args )

# to debug
setwd('/home/fracasso/data/ArjanProject/Arjan_fMRI_scripts/results')
BOUNDARYNAME='boundary03'
DISTANCEROI='test_2_2D.1D.roi'
PARENTVOLUME='targetVolume_results_zp.nii.gz'
DEPTH ='target_depth.nii.gz'
args <- c( BOUNDARYNAME, DISTANCEROI, PARENTVOLUME, DEPTH, '/packages/afni/16.1.13','/home/fracasso/analysisAfni/surfaces')

source( sprintf('%s/AFNIio.R', args[5] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[6] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[6] ) )
library(RANN)


# ################## #  
# parse distance roi #
# ################## #

print('parsing roi file...')

startLineN_nodesList <- scan( file=args[2], skip=2, nlines=1, n=4, what='character'  )
startLineN_nodes <- as.numeric( startLineN_nodesList[4] )
startLine_nodesList <- scan( file=args[2], skip=9, nlines=startLineN_nodes, what=list('numeric','numeric') )
startLine_nodes <- as.numeric( startLine_nodesList[[1]] )

skip_line_2 <- startLineN_nodes + 9 + 2 + 2
startLineN_nodesList <- scan( file=args[2], skip=skip_line_2, nlines=1, n=4, what='character'  )
startLineN_nodes02 <- as.numeric( startLineN_nodesList[4] )
startLine_nodesList02 <- scan( file=args[2], skip=skip_line_2+7, nlines=startLineN_nodes02, what=list('numeric','numeric') )
startLine_nodes02 <- as.numeric( startLine_nodesList02[[1]] )


roi_nodes_skip <- startLineN_nodes + startLineN_nodes02 + (9 + 2)*2 + 2
roi_nodesList <- scan( file=args[2], skip=roi_nodes_skip, nlines=1, n=4, what='character'  )
roiN_nodes <- as.numeric( roi_nodesList[4] )
roi_nodes_skip_to_nodes <- startLineN_nodes + startLineN_nodes02 + (9 + 2)*2 + 9
roi_nodesList <- scan( file=args[2], skip=roi_nodes_skip_to_nodes, nlines=roiN_nodes, what=list('numeric','numeric') )
roi_nodes <- as.numeric( roi_nodesList[[1]] )


for ( indStartNodes in c(1,2) ) { 
  
  if (indStartNodes==1) { selStartNodes <- startLine_nodes }
  if (indStartNodes==2) { selStartNodes <- startLine_nodes02 }
  
  
  #selStartNodes <- startLine_nodes02
  
  print('downsampling roi nodes...')
  
  sel_roi_nodes <- roi_nodes[ seq( 1, length(roi_nodes), by=3 ) ] # donwsampling roi nodes, otherwise there are too many nodes -> too long
  rep_roi_nodes_start <- rep( sel_roi_nodes, length(selStartNodes)  )
  
  startNodes <- rep( selStartNodes, rep( length(sel_roi_nodes), length(selStartNodes) ) )
  
  write.table( cbind( startNodes, rep_roi_nodes_start ), file='datasetStart.1D', col.names=FALSE, row.names=FALSE  )
  
  print('compute distance from starting nodes to all the roi nodes selected after downsampling...')
  
  instr <- sprintf( '%s/SurfDist -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -input datasetStart.1D > exampleStart.1D', args[5], args[1], args[1], selStartNodes )
  system.time( system( instr ) )
  
  print('compute smallest distance from each node (downsampled) and the start line...')
  
  fromStart <- read.table('exampleStart.1D')
  
  uniqueRoiPoints <- unique( fromStart[,2] )
  mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
  for (k in 1:length( uniqueRoiPoints) ) {
    indexSel <- which( fromStart[,2]==uniqueRoiPoints[k] )
    selDF <- fromStart[indexSel,]
    indexMin <- which( selDF[,3]==min( selDF[,3] ) )[1]
    mapOutStart[k,] <- as.numeric( selDF[indexMin,c(2,3)] )
  }
  
  print('loading surface coordinates...')
  
  # load coords 
  coords <- read.table( sprintf( 'surfaces_folder/%s_sm.1D.coord', args[1] ) )
  coordsIndices <- seq( 0, dim(coords)[1]-1 )
  emptyMap <- array(0, c( length(coordsIndices), 2) )
  emptyMap[ mapOutStart[,1]+1, 1 ] <- 1
  emptyMap[ mapOutStart[,1]+1, 2 ] <- mapOutStart[,2]
  outMap <- cbind( coordsIndices, emptyMap )
  
  print('writing distances computed on the downsampled nodes...')
  
  write.table( outMap, file='outMap.1D', row.names=FALSE, col.names=FALSE )
  write.table( mapOutStart[,1], file='outMap_indices.1D', row.names=FALSE, col.names=FALSE )
  
  print('interpolate over the surface (kernel diameter = 3mm )...')
  
  system( 'mkdir surfInterp' )
  instr <- sprintf( '%s/ROIgrow -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -roi_labels PER_NODE -roi_nodes outMap_indices.1D -lim 1.5 -prefix roiMap001Fill', args[5], args[1], args[1] )
  system( instr )
  system('mv roiMap001Fill* surfInterp/')
  setwd('surfInterp/')
  files <- dir(pattern='roiMap001Fill*')
  
  for ( k in 1:length( files ) ) {
    fileContent <- read.table( files[k], as.is=TRUE )
    distValues <- emptyMap[ fileContent[,1], 2 ]
    meanValue <- median( distValues[distValues>0] ) # mean distance value != 0 inside the kernel around the selected node
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
    meanDist <- median( storeDataMean[ indexNodes ], na.rm=TRUE )
    storeNodes[k] <- meanDist
  } 
  
  # clean up in case of NAs
  
  if ( sum( is.na( storeNodes ) )>0 ) {
    indxNan <- which( is.na( storeNodes ) )
    storeNodes[indxNan] <- 0
  }
  
  outMap <- cbind( uniqueNodes, storeNodes )
  outName <- strsplit( as.character( args[2] ), split='.1D')
  
  if (indStartNodes==1) { outNameFile <- sprintf('%s_distance01.1D', outName[[1]][1]) }
  if (indStartNodes==2) { outNameFile <- sprintf('%s_distance02.1D', outName[[1]][1]) }
  
  #outNameFile <- sprintf('%s_distance.1D', outName[[1]][1])
  
  print( sprintf( 'write interpolated distance map (filename = %s)...', outNameFile) )
  
  write.table( outMap, file=sprintf('../%s', outNameFile), row.names=FALSE, col.names=FALSE )
  
  print('cleaning up ')
  setwd('..')
  system('rm -R surfInterp/')
  system('rm outMap.1D')
  system('rm outMap_indices.1D')
  system('rm exampleStart.1D')
  
}




outNameFileVolume <- sprintf('%s_distance.nii.gz', outName[[1]][1])
print( sprintf( 'write interpolated distance map on a volume, only around the selected surface (filename = %s)...', outNameFileVolume) )
commandLine <- sprintf( '%s/3dSurf2Vol -spec surfaces_folder/spec.surfaces.smoothed -surf_A surfaces_folder/%s_sm.1D.coord -sv %s -grid_parent %s -map_func ave -prefix %s -sdata_1D %s', args[5], args[1], args[3], args[3], outNameFileVolume, outNameFile )
system( commandLine )


print( sprintf( 'load parent, depth and the newly created distance volume (%s)...', outNameFileVolume) )
parentVolume <- read.AFNI( args[3] )
depthVolume <- read.AFNI( args[4] )
distanceVolume <- read.AFNI( outNameFileVolume )

parentData <- parentVolume$brk[,,,1]
depthData <- depthVolume$brk[,,,1]
distanceData <- distanceVolume$brk[,,,1]
emptyVolume <- distanceData

print( sprintf( 'project the newly created distance volume across depth...', outNameFileVolume) )

stepLimit <- 0.1
limit1 <- seq(0,1-stepLimit,by=stepLimit)
limit2 <- limit1 + stepLimit
for (lim in 1:length(limit1) ) {
  ind1 <- which( distanceData!=0 )
  ind2 <- which( depthData>limit1[lim] & depthData<=limit2[lim] )
  coordsWM <- matrix( t( coordinateFromLinearIndex( ind1, dim( emptyVolume ) ) ), ncol=3 )
  coordsGM <- matrix( t( coordinateFromLinearIndex( ind2, dim( emptyVolume ) ) ), ncol=3 )
  nnOut <- nn2( coordsGM, coordsWM,  k=2 )
  coords <- coordsGM[ array(nnOut$nn.idx), ]
  indexVol <- linearIndexFromCoordinate( t( coords ), dim( emptyVolume ) )
  emptyVolume[ indexVol ] <- distanceData[ distanceData!=0 ]
}

outNameFileVolumeDepth <- sprintf('%s_distance_depth.nii.gz', outName[[1]][1])
print( sprintf( 'write distance volume across depth volume (file = %s)...', outNameFileVolumeDepth) )
write.AFNI(outNameFileVolumeDepth,
           brk=emptyVolume,
           label=NULL,
           view='+orig',
           orient=parentVolume$orient,
           origin=parentVolume$origin,
           delta=parentVolume$delta,
           defhead=parentVolume$NI_head )

