compute_normal <- function( vertex, face ) {
  
face <- t( face )  
vertex <- t( vertex )
eps <- 1E-12

crossp <- function( x, y ) {
  z <- x;
  z[1,] = x[2,]*y[3,] - x[3,]*y[2,];
  z[2,] = x[3,]*y[1,] - x[1,]*y[3,];
  z[3,] = x[1,]*y[2,] - x[2,]*y[1,];
  return( z )
}
  
nface <- dim(face)[2]
nvert <- dim(vertex)[2]
normal <- array( 0, c(3, nvert) )

normalf <- crossp( vertex[,face[2,]]-vertex[,face[1,]], vertex[,face[3,]]-vertex[,face[1,]] );
d <- sqrt( apply( normalf^2, 2, sum )  )
d[d<eps] <- 1
normalf <- normalf / rbind( d, d, d )

normal <- array( 0, c(3, nvert) )
for (i in 1:nface) {
  f <- face[,i];
  for (j in 1:3) {
    normal[,f[j]] <- normal[,f[j]] + normalf[,i]
  }
}

# normalize
d <- sqrt( apply( normal^2, 2, sum )  )
d[d<eps] <- 1
normal <- normal / rbind( d, d, d )

# enforce that the normal are outward
mVertex <- apply( vertex, 2, mean )
v <- vertex - rbind( mVertex, mVertex, mVertex )
s <- apply( v*normal, 1, sum )
if ( sum(s>0)<sum(s<0) ) {
  #flip
  normal = -normal;
  normalf = -normalf;
}

return( normal )

}