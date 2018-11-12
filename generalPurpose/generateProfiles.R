generateProfiles <- function( surfValList, mapValList, extractIndex ) {  
  ## generate profiles from intensity files and roi files
  
  nRoiFiles <- length( surfValList )
  intensityValuesMat <- array(0, c( dim(surfValList[[1]])[1], nRoiFiles ) )
  
  for (k in 1:nRoiFiles) {
    
    print( sprintf('reading data from boundary %d', k ) )
    
    map <- as.data.frame( mapValList[[k]] )
    
    surfVal <- as.data.frame( surfValList[[k]] )
    
    unconnectedPoints <- which( surfVal[,1]==-1 | surfVal[,1]==1 | surfVal[,1]==0 )
    
    surfValClean <- ifelse( (surfVal[,1]==-1 | surfVal[,1]==1 | surfVal[,1]==0), 1, surfVal[,1] ) + 1
    
    if ( length( unconnectedPoints )>0 && exists('storeUnconnectedPoints') ) {
      print( unconnectedPoints )
      storeUnconnectedPoints <- c( storeUnconnectedPoints, unconnectedPoints )
    }    
    if ( length( unconnectedPoints )>0 && !exists('storeUnconnectedPoints') ) {
      print( unconnectedPoints )
      storeUnconnectedPoints <- unconnectedPoints
    }
    
    intValue <- map[ surfValClean, extractIndex ]  
    
    intensityValuesMat[,k] <- intValue
    
  }
  
  
  if (exists('storeUnconnectedPoints')) { 
    print( storeUnconnectedPoints ) 
    out <- list( intensity=intensityValuesMat,
                 unconnectedPoints=storeUnconnectedPoints )  
  }
  if (!exists('storeUnconnectedPoints')) { 
    out <- list( intensity=intensityValuesMat,
                 unconnectedPoints=0 )
  }
  
  return( out )
  
}
