args <- commandArgs(T)
print( args )

#rm(list=ls())
#setwd('/home/fracasso/data/ArjanProject/2017-01-27_HeadcoilSurf_R3_fMRI_AF/ParRec/results')
#args <- c('boundariesThr.nii.gz','46','1-3-5','800')

inputBoundaries <- args[1]
smoothPar <- args[2]
inflateSurfaces <- strsplit( args[3], '-' )[[1]]
inflateIterations <- args[4]
remeshFraction <- args[5]
instr <- sprintf( '3dinfo -nv %s > _ttt_npoints.1D', inputBoundaries )
system( instr )
nSurfaces <- as.numeric( read.table( '_ttt_npoints.1D', as.is=TRUE ) )

## generate all isosurfaces (unsmoothed)
for ( k in 1:nSurfaces ) {
  instr <- sprintf( 'IsoSurface -input %s[%s] -remesh %s -isoval 1 -o_vec boundary0%s_or', inputBoundaries, k-1, remeshFraction, k-1 )
  print( instr )
  system( instr )
}

## prepare spec file (unsmoothed)
coordFiles <- dir(pattern='*.coord')
topoFiles <- dir(pattern='*.topo')
instr <- 'quickspec -spec spec.surfaces'
for ( k in 0:(nSurfaces-1) ) {
  part <- sprintf(' -tn 1D %s %s ', coordFiles[k+1], topoFiles[k+1] )
  instr <- paste( c(instr, part), collapse = '' )
  print(instr)
} 
system( instr )

## smooth surfaces
for (k in 0:(nSurfaces-1)) {
  if ( k < 10 ) {
    instr <- sprintf('SurfSmooth -spec spec.surfaces -surf_A %s -met LM -Niter %s -output boundary0%s_sm.1D.coord', coordFiles[k+1], smoothPar, k)
  }
  if ( k >= 10 ) {
    instr <- sprintf('SurfSmooth -spec spec.surfaces -surf_A %s -met LM -Niter %s -output boundary%s_sm.1D.coord', coordFiles[k+1], smoothPar, k)
  }
  print( instr )
  system( instr )
}
  
## inflate surfaces
for (k in 1:length(inflateSurfaces) ) {
  boundaryFile <- coordFiles[ as.numeric( inflateSurfaces[k] ) + 1 ]
  inflatedName <- strsplit( boundaryFile, '_' )[[1]][1]
  instr <- c( sprintf('SurfSmooth -spec spec.surfaces -surf_A %s -met NN_geom ', coordFiles[ as.numeric( inflateSurfaces[k] ) + 1 ] ), 
              sprintf('-output \u0027inflated_%s.1D.coord\u0027 -Niter %s -match_vol 0.01', inflatedName, inflateIterations ) )
  
#  instr <- c( sprintf('SurfSmooth -spec spec.surfaces -surf_A %s -met HEAT_05 ', coordFiles[ as.numeric( inflateSurfaces[k] ) + 1 ] ), 
#              sprintf('-output \u0027inflated_%s.1D.coord\u0027 -Niter %s', inflatedName, inflateIterations ) )
  
  instr <- paste( instr, collapse='' )
  print( instr )
  system( instr )
}

## prepare spec file (smoothed)
coordFiles <- dir(pattern='*_sm.1D.coord')
topoFiles <- dir(pattern='*.topo')
inflatedSurface <- dir(pattern='inflated_*')
coordFilesInflated <- coordFiles[ as.numeric( inflateSurfaces ) + 1 ]
topoFilesInflated <- topoFiles[ as.numeric( inflateSurfaces ) + 1 ]
counterInflatedSurfaces <- 1
instr <- 'quickspec -spec spec.surfaces.smoothed'
for ( k in 0:(nSurfaces+length(inflateSurfaces)-1) ) {
  if (k < nSurfaces) {
    part <- sprintf(' -tn 1D %s %s ', coordFiles[k+1], topoFiles[k+1] )
    instr <- paste( c(instr, part), collapse = '' )
    print(instr)
  }
  if (k >= nSurfaces) {
    nameSplit <- strsplit( topoFilesInflated[counterInflatedSurfaces], 'boundary' )[[1]][2]
    newName <- sprintf('inflated_%s',nameSplit)
    instrCopy <- sprintf( 'cp %s %s', topoFilesInflated[counterInflatedSurfaces], newName )
    system( instrCopy )
    part <- sprintf(' -tsnad 1D S_inflated_%s %s %s N %s', counterInflatedSurfaces, inflatedSurface[counterInflatedSurfaces], newName, coordFilesInflated[counterInflatedSurfaces])
    instr <- paste( c(instr, part), collapse = '' )
    print(instr)
    counterInflatedSurfaces <- counterInflatedSurfaces + 1
  }
}
system( instr )
system('rm *_or.1D.coord spec.surfaces')

print('INFLATED SURFACES:')
for (k in 1:length(inflateSurfaces)) {
  print( inflatedSurface[k] )
}

## clean up
system('mkdir surfaces_folder')
system('mv *.1D.topo surfaces_folder/')
system('mv *.1D.coord surfaces_folder/')
system('mv spec.surfaces.smoothed surfaces_folder/')
system('rm _ttt*')

