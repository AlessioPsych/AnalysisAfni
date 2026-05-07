args <- commandArgs(T)
print( args )

# to debug
#setwd('/analyse/Project0226/HNA10/ANATOMY')
#afniInstallDir <- Sys.getenv(x='AFNI_INSTALLDIR')
#args <- c( 'boundary00', 'testRoi.1D.roi', 'surfaces_folder_right/', '2.5', 'roiFilled.1D.roi' )

boundary <- args[1]
roiFile <- args[2]
surfacesDir <- args[3]
diameter <- args[4]
outputName <- args[5]

# ################## #  
# parse distance roi #
# ################## #

print('parsing roi file...')
# parsing the initial roi file
roi_dimen_char <- grep("ni_dimen", readLines( roiFile ), value = TRUE)
for (k in 1:length( roi_dimen_char ) ) {
  roi_dimen_loop <- as.numeric( gsub( "[^0-9.]", "",  roi_dimen_char[k] ) )
  if (k == 1) {
    nodesList <- list()
    nodesList[[k]] <- scan( file=roiFile, skip=9, nlines=roi_dimen_loop, what=list('numeric','numeric') )[[1]]
    nextLine <- 9+roi_dimen_loop+11
  }
  if (k>1) {
    nodesList[[k]] <- scan( file=roiFile, skip=nextLine, nlines=roi_dimen_loop, what=list('numeric','numeric') )[[1]]
    nextLine <- nextLine+roi_dimen_loop+11
  }
}

# growing the initial roi file
for ( k in 1:length( nodesList ) ) {
  nodesLoop <- as.numeric( nodesList[[k]] )
  print( sprintf( 'iteration %1.0f out of %1.0f, n. nodes: %1.0f', k, length( nodesList ), length( nodesLoop ) )  )
  write.table( nodesLoop, file = '_ttt_tempRoi.1D', row.names = FALSE, col.names = FALSE  )
  instr <- sprintf( 'ROIgrow -i_1D %s%s_sm.1D.coord %s%s_or.1D.topo -roi_nodes _ttt_tempRoi.1D -lim %s -prefix _ttt_tempRoiFill', surfacesDir, boundary, surfacesDir, boundary, diameter )
  system( instr )
  if (k == 1) {
    roiGrowTemp <- read.table( file='_ttt_tempRoiFill.1D', as.is=TRUE )
    roiGrowTemp <- roiGrowTemp[,1]
  }
  if (k > 1) {
    roiGrowTempLoop <- read.table( file='_ttt_tempRoiFill.1D', as.is=TRUE )
    roiGrowTemp <- unique( c( roiGrowTemp, roiGrowTempLoop[,1] )  )
  }
  
  system('rm _ttt_tempRoi.1D')
  system('rm _ttt_tempRoiFill.1D')
}

# save data as a surface ROI file for SUMA
write.table( cbind( roiGrowTemp, rep(1,length(roiGrowTemp) ) ), file=outputName, row.names = FALSE, col.names = FALSE )  




