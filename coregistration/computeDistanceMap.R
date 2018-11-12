computeDistanceMap <- function( maskFile ) {
  
  system( sprintf( '3dcopy %s ttt_dataFile+orig', maskFile ) )
  #system( sprintf('3drefit -anat ttt_dataFile+orig') )
  
  #system( '3dcalc -a ttt_dataFile+orig -expr \u0027step(a)\u0027 -datum short -prefix ttt_dataFile01+orig'  )
  
  system( '3dcalc -a ttt_dataFile+orig -expr \u0027step(a)\u0027 -datum float -prefix ttt_dataFile01+orig'  )
  
  #instr <- sprintf( 'IsoSurface -input ttt_dataFile+orig -isoval 1 -o_vec ttt_boundary_distance' )
  #system( instr )
  #instr <- 'quickspec -spec ttt_spec.surfaces -tn 1D ttt_boundary_distance.1D.coord ttt_boundary_distance.1D.topo'
  #system( instr )
  #instr <- sprintf('SurfSmooth -spec ttt_spec.surfaces -surf_A ttt_boundary_distance.1D.coord -met LM -Niter 50 -output ttt_boundary_distance_sm.1D.coord')
  #system( instr )
  
  #surfCoords <- read.table( 'ttt_boundary_distance.1D.coord' )
  #matString <- as.numeric( as.character( innerFile$NI_head$IJK_TO_DICOM$dat ) )
  #convMat <- matrix( matString, ncol=4, nrow = 3, byrow=TRUE )
  #convMat <- rbind( convMat, c(0,0,0,1) )
  #invConvMat <- solve( convMat )
  #idMat <- matrix( c(1,0,0,0,1,0,0,0,1), ncol=3, byrow=TRUE )
  #idMat %*% convMat
  #coordMat <- surfCoords
  #coordMat <- cbind( coordMat, rep(1,dim(coordMat)[1]) )
  #tCoordMat <- matrix( drop( t(coordMat) ), nrow=4 )
  #voxCoord <- invConvMat %*% tCoordMat
  #voxCoord03 <- t( voxCoord[c(1:3),] )
  
  innerFile <- read.AFNI( 'ttt_dataFile01+orig' )
  innerVolume <- innerFile$brk[,,,1]
  
  inside <- which( innerVolume==1 )
  outside <- which( innerVolume==0 ) 
  
  coordsInside <- coordinateFromLinearIndex( inside, dim(innerVolume) )
  coordsOutside <- coordinateFromLinearIndex( outside, dim(innerVolume) )
  
  iPart <- t( coordsInside )
  oPart <- t( coordsOutside ) 
  #nnOut <- nn2( surfCoords, oPart, k=1 )
  #nnOutIn <- nn2( surfCoords, iPart, k=1 )
  #nnOut <- nn2( voxCoord03, oPart, k=1 )
  #nnOutIn <- nn2( voxCoord03, iPart, k=1 )
  
  nnOut <- nn2( iPart, oPart, k=1 )
  nnOutIn <- nn2( oPart, iPart, k=1 )
  
  emptyVol <- array( 0, dim(innerVolume) )
  distOut <- as.numeric( round( nnOut$nn.dists, 2 ) )
  distIn <- as.numeric( round( nnOutIn$nn.dists, 2 )*-1 )
  emptyVol[outside] <- distOut
  emptyVol[inside] <- distIn
  innerFile$brk[,,,1] <- emptyVol
  volFileName <- sprintf( 'ttt_distanceMap_overall+orig' )
  write.AFNI(volFileName,
             brk=emptyVol,
             label=NULL,
             view='+orig',
             orient=innerFile$orient,
             origin=innerFile$origin,
             defhead=innerFile$NI_head )
  
  
  
  system('rm ttt_dataFile*')
  
  
}