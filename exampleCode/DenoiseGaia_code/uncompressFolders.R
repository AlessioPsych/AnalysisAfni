mainDir <- '/media/alessiof/Expansion/Gaia_David_info/Proj_Gaia_David'
runFlag <- 0
foldersToProcess <- c('EPI_timeSlicedCorrected_Safe','EPI_timeSlicedCorrected_Threat','freesurfer')

for ( k in 1:length( foldersToProcess ) ) {
  setwd( mainDir )
  setwd( foldersToProcess[k] )
  print( getwd() )
  subDirs <- list.files(recursive = FALSE, full.names = FALSE )
  for ( i in 1 : length( subDirs ) ) {
    folderName <- strsplit( subDirs, '.tar.gz' )[[1]][1]
    if ( dir.exists( sprintf('%s', folderName ) ) ) {
      instr <- sprintf('rm -R %s', folderName )
      print( instr )
      if (runFlag==1) { system( instr ) }
    }
    instr <- sprintf('tar -zxvf %s.tar.gz', subDirs[i] )
    print( instr )
    if (runFlag==1) { system( instr ) }
  }
}

#to uncompress
#tar -zxvf sub03.tar.gz 