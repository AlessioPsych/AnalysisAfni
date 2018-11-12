args <- commandArgs(T)
print( args )

#stopifnot(FALSE)  
#args <- c('surface_node_normals', 'B0_angle', '0.2', '0.2', '0.2', '/home/alessiofracasso/Dropbox/analysisAfni/')

## checkdir, create if it does not exists, if it exists then halts execution
mainDir <- getwd()
newDir <- args[2]
flagDir <- dir.create( file.path(mainDir, newDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', newDir )
  warning( msg )
  stopifnot(flagDir)  
}

source( sprintf( '%s/load1DData.R', args[6] ) )

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

directory <- sprintf( args[1] ) #input directory

normals <- dir( directory, pattern='*.do' )
B0vector <- as.numeric( c( args[3], args[4], args[5] ) )

normalsData <- load1DData( directory, seq( 1, length(normals) ), 1, '*.do')

setwd(mainDir)

for (k in 1:length(normalsData) ) {
  data01 <- as.matrix( normalsData[[k]] )
  data01Mat <- data01[,2:4]
  out <- apply( data01Mat, 1, angle, y=B0vector  ) * 180 / pi
  if (k<10) { system( sprintf( 'echo surface n. 0%s, file: %s', k, normals[k] ) ) }
  if (k>=10) { system( sprintf( 'echo surface n. %s, file: %s', k, normals[k] ) ) }
  
  splitName <- strsplit( normals[k], '[.]' )
  outName <- sprintf( '%s/%s_theta.1D', newDir, splitName[[1]][1] )
  
  emptyMap <- array( 0, c( dim( data01Mat )[1], 2 ) )
  emptyMap[,1] <- seq( 1, dim( data01Mat )[1], by=1 )
  emptyMap[,2] <- round( out, 2 )
  write.table( emptyMap, file=outName, row.names=FALSE, col.names=FALSE )
  
}


