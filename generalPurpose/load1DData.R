load1DData <- function(directory, fileNumber, skipLines, filePattern) {  
  setwd( directory )
  fileSurfacesMapping <- dir( pattern=filePattern)
  mapValList <- list(0)  
  print( sprintf('%d files to load', length(fileNumber) ) )
  
  for ( k in 1:length(fileNumber) ) {
    fileIndex <- fileNumber[k]
    print( sprintf('reading file n. %d: %s ...', k, fileSurfacesMapping[ fileIndex ] ) )
    mapValList[[k]] <- read.table( fileSurfacesMapping[ fileIndex ], skip=skipLines, as.is=TRUE )        
  }
  
  print( 'done' )
  return( mapValList )  
  
}