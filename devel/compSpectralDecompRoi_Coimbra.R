args <- commandArgs(T)
print( args )

#to debug
# setwd('/media/zohartal/Data/Tools/Analysis/suma_MNI152_2009_princetonAtlas')
# DISTANCEROI='testRoi.1D.roi'
# NODEFACTR='5'
# KERNELRADIUS='1.75'
# SURFACENAME='std.141.lh.smoothwm.gii'
# ROTATION='0.2'
# OUTNAME='5_0.2_axis_lh'
# args <- c( DISTANCEROI, NODEFACTR, KERNELRADIUS, 
#            Sys.getenv(x='AFNI_INSTALLDIR'),
#            Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES'),
#            SURFACENAME, ROTATION, OUTNAME )

mainDir <- getwd()

set.seed(123)
source( sprintf('%s/AFNIio.R', args[4] ) )
source( sprintf('%s/coordinateFromLinearIndex.r', args[5] ) )
source( sprintf('%s/linearIndexFromCoordinate.r', args[5] ) )
library(RANN)

nodeDownsampling <- as.numeric( args[2] )
kernelRadius <- as.numeric( args[3] )
if ( nodeDownsampling == 1 ) { kernelRadius <- 0  }

# load ROI nodes
print('loading roi nodes...')
roi_file <- args[1]
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
instr <- sprintf( 'SurfDist -i_gii %s -input datasetStart.1D > exampleStart.1D', args[6] )
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
#outFileName <- strsplit( args[2], '[.]')[[1]][1]
outName <- args[8]
outFileName <- sprintf( '%s_axisTemp.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

setwd(mainDir)

if ( kernelRadius != 0 ) {
  # load coords 
  
  #coords <- read.table( sprintf( '%s%s_sm.1D.coord', args[7], args[1] ) )
  #coordsIndices <- seq( 0, dim(coords)[1]-1 )
  emptyMap <- array(0, c( 198812, 3) )
  emptyMap[ outTemp[,1]+1, 1 ] <- 1
  emptyMap[ outTemp[,1]+1, 2 ] <- outTemp[,2]
  emptyMap[ outTemp[,1]+1, 3 ] <- outTemp[,3]
  #outMap <- cbind( coordsIndices, emptyMap )
  
  print( sprintf( 'interpolate over the surface (kernel radius = %smm )...', kernelRadius ) )
  
  write.table( selStartNodes, row.names=FALSE, col.names = FALSE, 'outMap_indices.1D' )
  system( 'mkdir surfInterp' )
  instr <- sprintf( 'ROIgrow -i_gii %s -roi_labels PER_NODE -roi_nodes outMap_indices.1D -lim %s -prefix roiMap001Fill', args[6], kernelRadius )
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
  whichUniqueNodes <- unlist( sapply( seq( 1, length( roi_nodes) ), nodeIdx, startingNodes=roi_nodes, finalNodes=uniqueNodes ) )
  
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

# rotated maps
rot <- as.numeric( args[7] )

# rotated map, original

idxCenter <- which.min( abs( vertexF[1,] )  )
vertexFC <- vertexF - t( array( rep( vertexF[ , idxCenter ], c(n,n) ), c(n,2) ) )
vertexFCR <- rbind( vertexFC[1,]*cos(rot) + vertexFC[2,]*sin(rot) , vertexFC[1,]*sin(rot) + vertexFC[2,]*cos(rot) )

outTemp <- round( cbind( selStartNodes, t( vertexFC ) ), 3 )
#outFileName <- paste( strsplit( args[2], '[.]')[[1]][1], '_', as.character(rot), sep='' )
outName <- args[8]
outFileName <- sprintf( '%s_axisTemp_cent.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

outTemp <- round( cbind( selStartNodes, t( vertexFCR ) ), 3 )
outName <- args[8]
outFileName <- sprintf( '%s_axisTemp_cent_rot.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

# rotated map, interpolated
vertexF_in <- t( out[,c(2,3)] )
n_in <- dim(vertexF_in)[2]
idxCenter <- which.min( abs( vertexF_in[1,] )  )
vertexFC <- vertexF_in - t( array( rep( vertexF_in[ , idxCenter ], c(n_in,n_in) ), c(n_in,2) ) )
vertexFCR <- rbind( vertexFC[1,]*cos(rot) + vertexFC[2,]*sin(rot) , vertexFC[1,]*sin(rot) + vertexFC[2,]*cos(rot) )

outTemp <- round( cbind( out[,1], t( vertexFC ) ), 3 )
#outFileName <- paste( strsplit( args[2], '[.]')[[1]][1], '_', as.character(rot), sep='' )
outName <- args[8]
outFileName <- sprintf( '%s_axis_cent.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

outTemp <- round( cbind( out[,1], t( vertexFCR ) ), 3 )
#outFileName <- paste( strsplit( args[2], '[.]')[[1]][1], '_', as.character(rot), sep='' )
outName <- args[8]
outFileName <- sprintf( '%s_axis_cent_rot.dset', outName )
write.table( outTemp, row.name=FALSE, col.names=FALSE , outFileName )

# to test on terminal
# 'compSpectralDecompRoi_Coimbra.sh testRoi.1D.roi 5 1.75 std.141.lh.smoothwm.gii 0.4 5_0.2_axis_lh'

