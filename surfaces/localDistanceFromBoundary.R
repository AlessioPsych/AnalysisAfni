args <- commandArgs(T)
print( args )
mainDir <- getwd()

#rm(list=ls())
#mainDir <- '/home/fracasso/data/ArjanProject/Arjan_fMRI_scripts/results'
#setwd( mainDir )
#args <- c('roiBoundary00.1D.roi','0.1','boundary00','mapOnSurface_smoothed_8/','15','/packages/afni/16.1.13','/home/fracasso/analysisAfni')

source( sprintf('%s/AFNIio.R', args[6]) )
source( sprintf('%s/load1DData.R', args[7]) )
source( sprintf('%s/generateProfiles.R', args[7]) )

instr <- sprintf( 'cp %s ttt_tempRoi.1D.roi', args[1] )
system( instr )

# which column index does the boundary corresponds to for the ROI and surfaces dataset?
boundarySplit <- strsplit( args[3], '' )[[1]]
columnIndexSplit <- as.numeric( boundarySplit[ c( length( boundarySplit )-1, length( boundarySplit ) ) ] )
columnIndex <- ifelse( columnIndexSplit[1]==0, columnIndexSplit[2], as.numeric( paste( as.character( columnIndexSplit ), collapse = "" ) ) ) + 1

# which column index in the surface dataset are we looking for?
indexSurface <- as.numeric( args[5] )

setwd(mainDir)
roiData <- read.table( 'ttt_tempRoi.1D.roi' )
#instr <- sprintf( 'mappingBetweenSurfaces_ROI.sh ttt_tempRoi.1D.roi %s', as.character( columnIndex-1 ) )  
#system( instr )

setwd(mainDir)
setwd( args[4]  )
filesMap <- dir( pattern = '*.1D.dset' )
fileNumber <- seq(1,length(filesMap),1)
setwd(mainDir)
mapValListStat <- load1DData( args[4], fileNumber, 5, '*.1D.dset' )

#setwd(mainDir)
#fileNumber <- seq(1,length(filesMap),1)
#mapValListStatOriginal <- load1DData( 'mapOnSurface/', fileNumber, 5, '*.1D.dset' ) ##########################################

#setwd(mainDir)
#firstPart <- fileNumber[1:columnIndex]
#secondPart <- seq( (firstPart[ length(firstPart) ]+2), (length(fileNumber))+1, 1 )
#fileNumberRoi <- c( firstPart, secondPart )
#fileNumberRoi <- c(1,2,3,4,6,7,8)
#surfValListRoi <- load1DData( 'ttt_tempRoi_folder/', fileNumberRoi, 0, '^roiSurface.*.1D.roi' )

#setwd(mainDir)
#source( sprintf('%s/generateProfiles.R', args[6]) )
#pRoi_01 <- generateProfiles( surfValListRoi, mapValListStat, indexSurface )
#hist( pRoi_01$intensity[,4] )

###### ###### #######
#selNodesIdx <- which( abs( pRoi_01$intensity[,columnIndex] ) < as.numeric( args[2] ) )
#selNodes <- surfValListRoi[[columnIndex]][selNodesIdx,1]
#write.table( cbind( selNodes, rep(1, length(selNodes) ) ), file = 'ttt_bound.1D.roi', row.names = FALSE, col.names = FALSE  ) 
setwd(mainDir)
roiIntensity <- mapValListStat[[columnIndex]][ roiData[,1]+1  ,indexSurface]
selNodesIdx <- which( abs( roiIntensity ) <= as.numeric( args[2] ) )
selNodes <- roiData[selNodesIdx,1]
write.table( cbind( selNodes, rep(1, length(selNodes) ) ), file = 'ttt_bound.1D.roi', row.names = FALSE, col.names = FALSE  ) 



#selNodesIdx <- mapValListStat[[columnIndex]][,indexSurface]

#dim( pRoi_01$intensity )

#selNodes <- surfValListRoi[[1]][,1]
#write.table( cbind(selNodes, outSetMiddle02Filt, corrDist ), file = 'dist00.dset', row.names = FALSE, col.names = FALSE  ) 
# selNodesIdx <- which( abs( pRoi_01$intensity[,5] ) < as.numeric( args[2] ) )
# selNodes <- surfValListRoi[[5]][selNodesIdx,1]
# write.table( cbind( selNodes, rep(1, length(selNodes) ) ), file = 'ttt_bound_04.1D.roi', row.names = FALSE, col.names = FALSE  ) 
# selNodesIdx <- which( abs( pRoi_01$intensity[,3] ) < as.numeric( args[2] ) )
# selNodes <- surfValListRoi[[3]][selNodesIdx,1]
# write.table( cbind( selNodes, rep(1, length(selNodes) ) ), file = 'ttt_bound_02.1D.roi', row.names = FALSE, col.names = FALSE  ) 
# selNodesIdx <- which( abs( pRoi_01$intensity[,2] ) < as.numeric( args[2] ) )
# selNodes <- surfValListRoi[[2]][selNodesIdx,1]
# write.table( cbind( selNodes, rep(1, length(selNodes) ) ), file = 'ttt_bound_01.1D.roi', row.names = FALSE, col.names = FALSE  ) 


# negative distances from border

#selNodesIdx <- which( pRoi_01$intensity[,4] < 0.1 )
#selNodes <- surfValListRoi[[4]][selNodesIdx,1]

## first step

setwd(mainDir)
roiIntensity <- mapValListStat[[columnIndex]][ roiData[,1]+1  ,indexSurface]
selNodesIdx <- which( roiIntensity < as.numeric( args[2] ) )

#selNodesIdx <- which( pRoi_01$intensity[,columnIndex] < as.numeric( args[2] ) )
selStartNodes <- selNodes
sel_roi_nodes <- roiData[selNodesIdx,1]
rep_roi_nodes_start <- rep( sel_roi_nodes, length(selStartNodes)  )
startNodes <- rep( selStartNodes, rep( length(sel_roi_nodes), length(selStartNodes) ) )
write.table( cbind( startNodes, rep_roi_nodes_start ), file='ttt_datasetStart.1D', col.names=FALSE, row.names=FALSE  )

print('compute distance from starting nodes to positive roi nodes...')

instr <- sprintf( 'SurfDist -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -input ttt_datasetStart.1D > ttt_exampleStart.1D', args[3], args[3] )
system.time( system( instr ) )

print('compute smallest distance from each node and the start line...')

fromStart <- read.table('ttt_exampleStart.1D')

uniqueRoiPoints <- unique( fromStart[,2] )
mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
for (k in 1:length( uniqueRoiPoints) ) {
  indexSel <- which( fromStart[,2]==uniqueRoiPoints[k] )
  selDF <- fromStart[indexSel,]
  indexMin <- which( selDF[,3]==min( selDF[,3] ) )[1]
  mapOutStart[k,] <- as.numeric( selDF[indexMin,c(2,3)] )
}

mat01 <- cbind( sel_roi_nodes, rep(1, length(sel_roi_nodes) ), mapOutStart )
write.table( mat01, file = 'ttt_boundDist01.dset', row.names = FALSE, col.names = FALSE  ) 
#write.table( mat01, file = 'ttt_delMe.dset', row.names = FALSE, col.names = FALSE  ) 


#mat01 <- cbind( sel_roi_nodes, rep(1, length(sel_roi_nodes) ), mapOutStart )
#write.table( mat01, file = 'delMe001.dset', row.names = FALSE, col.names = FALSE  ) 

#mat01 <- cbind( sel_roi_nodes, rep(1, length(sel_roi_nodes) ), round( rnorm( length(sel_roi_nodes) ), 2 ) )
#write.table( mat01, file = 'delMe00.dset', row.names = FALSE, col.names = FALSE  ) 

# second step

setwd(mainDir)
roiIntensity <- mapValListStat[[columnIndex]][ roiData[,1]+1  ,indexSurface]
selNodesIdx <- which( roiIntensity > as.numeric( args[2] ) )

#selNodesIdx <- which( pRoi_01$intensity[,columnIndex] > 0.1 )
selStartNodes <- selNodes
sel_roi_nodes <- roiData[selNodesIdx,1]
rep_roi_nodes_start <- rep( sel_roi_nodes, length(selStartNodes)  )
startNodes <- rep( selStartNodes, rep( length(sel_roi_nodes), length(selStartNodes) ) )
write.table( cbind( startNodes, rep_roi_nodes_start ), file='ttt_datasetStart.1D', col.names=FALSE, row.names=FALSE  )

print('compute distance from starting nodes to negative roi nodes...')

instr <- sprintf( 'SurfDist -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -input ttt_datasetStart.1D > ttt_exampleStart.1D', args[3], args[3] )
system.time( system( instr ) )

print('compute smallest distance from each node and the start line...')

fromStart <- read.table('ttt_exampleStart.1D')

uniqueRoiPoints <- unique( fromStart[,2] )
mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
for (k in 1:length( uniqueRoiPoints) ) {
  indexSel <- which( fromStart[,2]==uniqueRoiPoints[k] )
  selDF <- fromStart[indexSel,]
  indexMin <- which( selDF[,3]==min( selDF[,3] ) )[1]
  mapOutStart[k,] <- as.numeric( selDF[indexMin,c(2,3)] )
}

mat02 <- cbind( sel_roi_nodes, rep(1, length(sel_roi_nodes) ), mapOutStart )
write.table( mat02, file = 'ttt_boundDist02.dset', row.names = FALSE, col.names = FALSE  ) 


### third step
setwd(mainDir)
firstSet <- read.table( 'ttt_boundDist01.dset' ) 
selFirstSet <- firstSet[ firstSet[,4] < 1, 1 ]


roiIntensity <- mapValListStat[[columnIndex]][ roiData[,1]+1  ,indexSurface]
selNodesIdx <- which( abs(roiIntensity) < as.numeric( args[2] ) )

#selNodesIdx <- which( abs(pRoi_01$intensity[,columnIndex]) < 0.1 )
selStartNodes <- selFirstSet
sel_roi_nodes <- roiData[selNodesIdx,1]
rep_roi_nodes_start <- rep( sel_roi_nodes, length(selStartNodes)  )
startNodes <- rep( selStartNodes, rep( length(sel_roi_nodes), length(selStartNodes) ) )
write.table( cbind( startNodes, rep_roi_nodes_start ), file='ttt_datasetStart.1D', col.names=FALSE, row.names=FALSE  )

print('compute distance from starting nodes in the boundary - step 01...')

instr <- sprintf( 'SurfDist -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -input ttt_datasetStart.1D > ttt_exampleStart.1D', args[3], args[3] )
system.time( system( instr ) )

print('compute smallest distance from each node and the start line...')

fromStart <- read.table('ttt_exampleStart.1D')

uniqueRoiPoints <- unique( fromStart[,2] )
mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
for (k in 1:length( uniqueRoiPoints) ) {
  indexSel <- which( fromStart[,2]==uniqueRoiPoints[k] )
  selDF <- fromStart[indexSel,]
  selDFSortedIdx <- sort( selDF[,3], index.return = TRUE  )
  selDFSorted <- selDF[ selDFSortedIdx$ix, ]
  indexMin <- which( selDFSorted[,3]>0 )[1]
  mapOutStart[k,] <- as.numeric( selDFSorted[indexMin,c(2,3)] )
}

mat02 <- cbind( sel_roi_nodes, rep(1, length(sel_roi_nodes) ), mapOutStart )
write.table( mat02, file = 'ttt_boundDist03.dset', row.names = FALSE, col.names = FALSE  ) 


#fourth step

firstSet <- read.table( 'ttt_boundDist02.dset' ) 
selFirstSet <- firstSet[ firstSet[,4] < 1, 1 ]

roiIntensity <- mapValListStat[[columnIndex]][ roiData[,1]+1  ,indexSurface]
selNodesIdx <- which( abs(roiIntensity) < as.numeric( args[2] ) )

#selNodesIdx <- which( abs(pRoi_01$intensity[,columnIndex]) < 0.1 )
selStartNodes <- selFirstSet
sel_roi_nodes <- roiData[selNodesIdx,1]
rep_roi_nodes_start <- rep( sel_roi_nodes, length(selStartNodes)  )
startNodes <- rep( selStartNodes, rep( length(sel_roi_nodes), length(selStartNodes) ) )
write.table( cbind( startNodes, rep_roi_nodes_start ), file='ttt_datasetStart.1D', col.names=FALSE, row.names=FALSE  )

print('compute distance from starting nodes in the boundary - step 02...')

instr <- sprintf( 'SurfDist -i_1D surfaces_folder/%s_sm.1D.coord surfaces_folder/%s_or.1D.topo -input ttt_datasetStart.1D > ttt_exampleStart.1D', args[3], args[3] )
system.time( system( instr ) )

print('compute smallest distance from each node (downsampled) and the start line...')

fromStart <- read.table('ttt_exampleStart.1D')

uniqueRoiPoints <- unique( fromStart[,2] )
mapOutStart <- array( 0, c(length( uniqueRoiPoints),2) )
for (k in 1:length( uniqueRoiPoints) ) {
  indexSel <- which( fromStart[,2]==uniqueRoiPoints[k] )
  selDF <- fromStart[indexSel,]
  selDFSortedIdx <- sort( selDF[,3], index.return = TRUE  )
  selDFSorted <- selDF[ selDFSortedIdx$ix, ]
  indexMin <- which( selDFSorted[,3]>0 )[1]
  mapOutStart[k,] <- as.numeric( selDFSorted[indexMin,c(2,3)] )
}

mat02 <- cbind( sel_roi_nodes, rep(1, length(sel_roi_nodes) ), mapOutStart )
write.table( mat02, file = 'ttt_boundDist04.dset', row.names = FALSE, col.names = FALSE  ) 

# fifth step:

# smallest distance in the nodes in the middle
middleSet01 <- read.table( 'ttt_boundDist03.dset' ) 
middleSet02 <- read.table( 'ttt_boundDist04.dset' ) 
keepIndex <- middleSet01[,4] < middleSet02[,4]
middleSet01[,4] <- middleSet01[,4]*-1
keepIndexT <- which( keepIndex==TRUE )
keepIndexF <- which( keepIndex==FALSE )
outSetMiddle00 <- rbind( middleSet01[keepIndexT,], middleSet02[keepIndexF,] )
write.table( outSetMiddle00, file = 'ttt_boundMidd.dset', row.names = FALSE, col.names = FALSE  ) 

middleSet01 <- read.table( 'ttt_boundDist01.dset' ) 
middleSet02 <- read.table( 'ttt_boundDist02.dset' )
middleSet02[,4] <- middleSet02[,4]*-1
outSetMiddle01 <- rbind( middleSet01, middleSet02 )
write.table( outSetMiddle01, file = 'ttt_boundWhole.dset', row.names = FALSE, col.names = FALSE  ) 

outSetMiddle02 <- rbind( outSetMiddle00, outSetMiddle01 )

removeIdx <- which( outSetMiddle02[,4] == 0 )
outSetMiddle02Filt <- outSetMiddle02[ -removeIdx, ]
fNameOut <- strsplit( args[1], '[.]' )[[1]][1]
fNameOutDistance <- sprintf( '%s_boundaryDistance.dset', fNameOut )
#write.table( outSetMiddle02Filt, file = fNameOutDistance, row.names = FALSE, col.names = FALSE  ) 

system('rm -R ttt_*')

#storeDiff <- diff( sort( outSetMiddle02Filt[,4] ) )
#storeDiffFilt <- storeDiff[ storeDiff>0 & storeDiff>quantile(storeDiff,0.025) & storeDiff < quantile(storeDiff,0.975) ]
#hist( storeDiffFilt )
#mean( storeDiffFilt )

cBin <- .bincode( outSetMiddle02Filt[,4] , breaks=quantile( outSetMiddle02Filt[,4], seq(0,1,0.025) ), include.lowest=TRUE )
cMedian <- tapply( outSetMiddle02Filt[,4], list(cBin), median )
dcMedian <- diff( cMedian )
xVar <- cMedian[ seq( 1, length(dcMedian), 1 ) ]
ssp <- smooth.spline( xVar, dcMedian, spar=1 )
#plot( dcMedian~xVar )
#lines( ssp )
storeConst <- min( ssp$y )

corrDist <- round( ifelse( outSetMiddle02Filt[,4] > 0, outSetMiddle02Filt[,4]-storeConst/2, outSetMiddle02Filt[,4]+storeConst/2 ), 2 )

cBin <- .bincode( corrDist , breaks=quantile( corrDist, seq(0,1,0.025) ), include.lowest=TRUE )
cMedian <- tapply( corrDist, list(cBin), median )
dcMedian <- diff( cMedian )
xVar <- cMedian[ seq( 1, length(dcMedian), 1 ) ]
ssp <- smooth.spline( xVar, dcMedian, spar=1 )
#plot( dcMedian~xVar )
#lines( ssp )


write.table( cbind( outSetMiddle02Filt, corrDist ), file = fNameOutDistance, row.names = FALSE, col.names = FALSE  ) 






# rearrange a bit the code :
# and make a new function to:
# generate roi across depth,
# select nodes for each surface from the volumetric roi;
# fill up eventual empty nodes (smooth a bit);
# for each defined roi compute the metric over the surface;



# surfaceName <- 'boundary00'
# instr <- c( '3dVol2Surf', '-spec surfaces_folder/spec.surfaces.smoothed',
#             sprintf('-surf_A surfaces_folder/%s_sm.1D.coord', surfaceName),
#             '-sv leftSmall.1D.roi_clust+orig',
#             '-grid_parent leftSmall.1D.roi_clust+orig',
#             '-map_func mask',
#             '-out_1D delMeSurf04.1D.dset')
# instrConcatenated <- paste( instr, collapse = " "  )
# system( instrConcatenated )
 

