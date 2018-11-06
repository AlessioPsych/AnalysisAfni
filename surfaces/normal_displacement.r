normal_displacement <- function(vertex, face, displacement) {
  niter <- 5
  displacement <- displacement/niter;
  for (l in 1:niter) {
    normals <- compute_normal(vertex,face)
    vertex <- vertex + displacement*normals
  }
  return(vertex)
}