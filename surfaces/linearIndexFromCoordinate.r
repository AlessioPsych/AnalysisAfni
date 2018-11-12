linearIndexFromCoordinate <- function( coords, dims ) {
  indices <- coords[1,]
  for ( d in 2:length(dims) ) {
    indices <- indices + (coords[d,]-1) * prod(dims[1:d-1]);
  }
  return(indices)
}