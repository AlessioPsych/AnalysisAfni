coordinateFromLinearIndex <- function(indices, dims) {  
  coords <- array( 0, c( length(dims), length(indices) ) )
  dArray <- seq( length(dims), 1)
  for ( n in 1:length(dims) ) {
    d <- dArray[n]
    coords[d,] <- floor( (indices-1) / prod( dims[1:d-1] ) ) + 1;
    indices <- indices - (coords[d,]-1) * prod( dims[1:d-1] );
  }
  return( coords )
}