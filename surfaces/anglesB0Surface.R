#mainDirSet <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V5029leftright_Ben_ph/results_CBS'
#setwd(mainDirSet)

#args <- commandArgs(T)
#print( args )

newDir <- 'B0Angles_multiple'
## checkdir, create if it does not exists, if it exists then halts execution
mainDir <- getwd()
flagDir <- dir.create( file.path(mainDir, newDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s_folder to proceed', newDir )
  warning( msg )
  stopifnot(flagDir)  
}

possibleOrientations <- expand.grid( c(1,2,-2,-1), c(1,2,-2,-1), c(1,2,-2,-1) )

for ( k in 1:dim(possibleOrientations)[1] ) {
#for ( k in 1:2 ) {
  
  instr1 <- sprintf('B0_normals.sh surface_node_normals B0_angle %d %d %d', possibleOrientations[k,1], possibleOrientations[k,2], possibleOrientations[k,3])
  instr2 <- sprintf('smoothSurfaceMap.sh B0_angle/ 5')
  instr3 <- sprintf('mv B0_angle_smoothed_5/ B0_angle_smoothed_5_%d%d%d/', possibleOrientations[k,1], possibleOrientations[k,2], possibleOrientations[k,3])
  instr4 <- sprintf('mv B0_angle_smoothed_5_%d%d%d/ %s/', possibleOrientations[k,1], possibleOrientations[k,2], possibleOrientations[k,3], newDir)
  instr5 <- sprintf('rm -R B0_angle/')
  
  system( sprintf('echo %s', instr1) )
  system( instr1 )
  system( sprintf('echo %s', instr2) )
  system( instr2 )
  system( sprintf('echo %s', instr3) )
  system( instr3 )
  system( sprintf('echo %s', instr4) )
  system( instr4 )
  system( sprintf('echo %s', instr5) )
  system( instr5 )
  
}


