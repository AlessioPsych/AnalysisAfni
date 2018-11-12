args <- commandArgs(T)
print( args )

#to debug
#setwd('/home/fracasso/data/ArjanProject/2017-01-27_HeadcoilSurf_R3_fMRI_AF/ParRec/results')
#BOUNDARYNAME='boundary00'
#DISTANCEROI='testRoi.1D.roi'
#NODEFACTR=7
#KERNELRADIUS=0
#args <- c( BOUNDARYNAME, DISTANCEROI, NODEFACTR, KERNELRADIUS, '/packages/afni/17.0.13','/home/fracasso/analysisAfni/surfaces' )
# fix kernel option


mainDir <- getwd()

source( sprintf('%s/AFNIio.R', args[5] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[6] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[6] ) )
library(RANN)

nodeDownsampling <- as.numeric( args[3] )
kernelRadius <- args[4]
if ( nodeDownsampling == 1 ) { kernelRadius <- 0  }

# load ROI nodes
print('loading roi nodes...')
roi_file <- args[2]
roi_data <- read.table( roi_file, as.is=TRUE )
roi_nodes <- roi_data[,1]

# downsample ROI
print('downsampling roi nodes...')
sel_roi_nodes <- roi_nodes[ seq( 1, length(roi_nodes), by=nodeDownsampling ) ] # donwsampling roi nodes, otherwise there are too many nodes -> too long

# node to node distance matrix
print('node to node distance matrix (1st step, prepare the data)...')
selStartNodes <- sel_roi_nodes
rep_roi_nodes_start <- rep( sel_roi_nodes, length(selStartNodes)  )
startNodes <- rep( selStartNodes, rep( length(sel_roi_nodes), length(selStartNodes) ) )
write.table( cbind( startNodes, rep_roi_nodes_start ), file='datasetStart.1D', col.names=FALSE, row.names=FALSE  )

# node to node distance matrix
print('node to node distance matrix (2nd step, might take longer)...')
instr <- sprintf( '%s/SurfDist -i_1D %s%s_sm.1D.coord %s%s_or.1D.topo -input datasetStart.1D > exampleStart.1D', args[5], args[7], args[1], args[7], args[1] )
system.time( system( instr ) )

print('read node to node distance matrix')
nodeToNodeDistance <- read.table( 'exampleStart.1D', as.is=TRUE )
n <- length(selStartNodes)
D <- array( nodeToNodeDistance[,3], c( length(selStartNodes), length(selStartNodes) ) )
D <- ( D + t(D) ) / 2;
J <- diag(n) - array(1,dim(D))/n
W <- -J%*%(D^2)%*%J
eigenMat <- eigen( W )
S <- eigenMat$values
sortOut <- sort(S,decreasing = TRUE, index.return=TRUE)
U <- eigenMat$vectors[ , sortOut$ix ]
vertexF <- t( U[,1:2] ) * t( array( rep( sqrt( S[1:2] ), c(n,n) ), c(n,2) ) )

outTemp <- round( cbind( selStartNodes, t( vertexF ) ), 3 )
outFileName <- strsplit( args[2], '[.]')[[1]][1]
outFileName <- sprintf( '%s_axisTemp.dset', outFileName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

setwd(mainDir)

if ( kernelRadius != 0 ) {
  # load coords 
  coords <- read.table( sprintf( '%s%s_sm.1D.coord', args[7], args[1] ) )
  coordsIndices <- seq( 0, dim(coords)[1]-1 )
  emptyMap <- array(0, c( length(coordsIndices), 3) )
  emptyMap[ outTemp[,1]+1, 1 ] <- 1
  emptyMap[ outTemp[,1]+1, 2 ] <- outTemp[,2]
  emptyMap[ outTemp[,1]+1, 3 ] <- outTemp[,3]
  outMap <- cbind( coordsIndices, emptyMap )
  
  print( sprintf( 'interpolate over the surface (kernel radius = %smm )...', kernelRadius ) )
  
  write.table( selStartNodes, row.names=FALSE, col.names = FALSE, 'outMap_indices.1D' )
  system( 'mkdir surfInterp' )
  instr <- sprintf( 'ROIgrow -i_1D %s%s_sm.1D.coord %s%s_or.1D.topo -roi_labels PER_NODE -roi_nodes outMap_indices.1D -lim %s -prefix roiMap001Fill', args[7], args[1], args[7], args[1], kernelRadius )
  system( instr )
  system('mv roiMap001Fill* surfInterp/')
  setwd('surfInterp/')
  files <- dir(pattern='roiMap001Fill*')
  
  for ( k in 1:length( files ) ) {
    fileContent <- read.table( files[k], as.is=TRUE )
    
    # Axis 1
    distValues <- emptyMap[ fileContent[,1]+1, 2 ]
    meanValue01 <- median( distValues[distValues!=0] ) # mean distance value != 0 inside the kernel around the selected node
    meanValues01 <- rep( meanValue01, length( distValues ) ) # rep mean distance value != 0 inside the kernel around the selected node for the number of nodes where distValue != 0
    
    #Axis 2
    distValues <- emptyMap[ fileContent[,1]+1, 3 ]
    meanValue02 <- median( distValues[distValues!=0] ) # mean distance value != 0 inside the kernel around the selected node
    meanValues02 <- rep( meanValue02, length( distValues ) ) # rep mean distance value != 0 inside the kernel around the selected node for the number of nodes where distValue != 0
    
    if (k==1) {
      storeDataMean01 <- meanValues01
      storeDataMean02 <- meanValues02
      storeDataNodes <- fileContent[,1]
    }
    if (k>1) {
      storeDataMean01 <- c( storeDataMean01, meanValues01 )
      storeDataMean02 <- c( storeDataMean02, meanValues02 )
      storeDataNodes <- c( storeDataNodes, fileContent[,1] )
    }
  }
  
  # for all the unique nodes in storeDataNodes, extract the mean distance and report the average
  uniqueNodes <- unique( storeDataNodes ) 
  storeNodes01 <- rep(0,length(uniqueNodes))
  storeNodes02 <- rep(0,length(uniqueNodes))
  for (k in 1:length( uniqueNodes ) ) {
    indexNodes <- which( uniqueNodes[k]==storeDataNodes )
    meanDist01 <- median( storeDataMean01[ indexNodes ], na.rm=TRUE )
    meanDist02 <- median( storeDataMean02[ indexNodes ], na.rm=TRUE )
    storeNodes01[k] <- meanDist01
    storeNodes02[k] <- meanDist02
  } 
  
  # clean up in case of NAs
  if ( sum( is.na(storeNodes01) | is.na(storeNodes02) )>0 ) {
    indxNan <- which( is.na(storeNodes01) | is.na(storeNodes02) )
    storeNodes01[indxNan] <- 0
    storeNodes02[indxNan] <- 0
  }
  
  nodeIdx <- function( idx, startingNodes, finalNodes ) {
    return( which( startingNodes[idx]==finalNodes ) )
  }
  whichUniqueNodes <- sapply( seq( 1, length( roi_nodes) ), nodeIdx, startingNodes=roi_nodes, finalNodes=uniqueNodes )
  
  setwd( mainDir )
  out <- round( cbind( uniqueNodes[whichUniqueNodes], storeNodes01[whichUniqueNodes], storeNodes02[whichUniqueNodes] ), 2 )
  outFileName <- strsplit( args[2], '[.]')[[1]][1]
  outFileName <- sprintf( '%s_axis.dset', outFileName )
  write.table( out, row.name=FALSE, col.names=FALSE , outFileName )
  
  setwd( mainDir )
  print('cleaning up ')
  system('rm -R surfInterp/')
  system('rm outMap_indices.1D')
  system('rm exampleStart.1D')
  system('rm datasetStart.1D')
  
}

if ( kernelRadius == 0 ) {
  outFileName <- strsplit( args[2], '[.]')[[1]][1]
  outFileNameIn <- sprintf( '%s_axisTemp.dset', outFileName )
  outFileNameOut <- sprintf( '%s_axis.dset', outFileName )
  instr <- sprintf( 'mv %s %s', outFileNameIn, outFileNameOut )
  system( instr )
  setwd( mainDir )
  print('cleaning up ')
  system('rm exampleStart.1D')
  system('rm datasetStart.1D')
  
}


#idxCenter <- which( selStartNodes == centerNode )
#vertexFC <- vertexF - t( array( rep( vertexF[ , idxCenter ], c(n,n) ), c(n,2) ) )
#vertexFCR <- rbind( vertexFC[1,]*cos(rot) + vertexFC[2,]*sin(rot) , vertexFC[1,]*sin(rot) + vertexFC[2,]*cos(rot) )

#nIter <- 100
#vertexS <- vertexF
#stress <- rep(0,nIter)
#for ( k in 1:nIter ) {
#   D1 <- array( rep( apply( vertexS^2, 2, sum ), n ), c(n,n) )
#   D1 <- sqrt( D1 + t(D1) - 2*t(vertexS)%*%vertexS )
#   B <- -D/max(D1,1e-10)
#   B <- B - diag( apply( B, 1, sum ) )
#   vertexS <- t( (B%*%t(vertexS)) ) / n
#   stress[k] <- sqrt( sum( abs( array(D1)-array(D) )^2 ) / n^2 )
#}
 
plot( stress )


# read node to node distance
#print('read node to node distance matrix')
#nodeToNodeDistance <- read.table( 'exampleStart.1D', as.is=TRUE )
#nodeToNodeDistanceMat <- array( nodeToNodeDistance[,3], c( length(selStartNodes), length(selStartNodes) ) )
#nodeToNodeDistanceMatOver <- 1 / (nodeToNodeDistanceMat)
#diag( nodeToNodeDistanceMatOver ) <- 0
#eigenMat <- eigen( nodeToNodeDistanceMatOver )
#eigenMat01 <- eigen( nodeToNodeDistanceMat )
#out <- round( cbind( selStartNodes, eigenMat$vectors[ ,c(2,3,4,5,6) ] ), 3 )
#write.table( out, row.name=FALSE, col.names=FALSE , 'axis04.dset' )
#out <- round( cbind( selStartNodes, eigenMat01$vectors[ ,c(2,3,4,5,6) ] ), 3 )
#write.table( out, row.name=FALSE, col.names=FALSE , 'axis05.dset' )



# nodeToNodeDistanceMat <- array( nodeToNodeDistance[,3], c( length(selStartNodes), length(selStartNodes) ) )
# W <- 1 / (nodeToNodeDistanceMat)
# diag( W ) <- 0
# D <- array( 0, dim( W ) )
# diag(D) <- apply( W, 2, sum )
# G <- solve( D )
# L <- G%*%(D-W)
# eigenMat <- eigen( L )
# 
# outVect <- round( eigenMat$vectors[,c(2,3)]*1e4, 2 )
# 
# #G <- array( 0, dim( W ) )
# #diag( G ) <- 1/diag(D)
# #Gp <- G^2
# #L <- G*(D-W)
# 
# #eigenMat <- eigen( L )
#   
# out <-  cbind( selStartNodes, outVect  )
# #out <- round( cbind( selStartNodes, eigenMat$vectors[ ,c(1,2) ] ), 3 )
# write.table( out, row.name=FALSE, col.names=FALSE , 'axis02.dset' )
# 
# eigenMat$values
# 
# 
# vals <- eigenMat$vectors[ ,c(2,3) ]
# hist(vals[,1])
# hist(vals[,2])

#eigenMat$vectors[ ,c(2,3) ] * eigenMat$values[c(2,3)]

