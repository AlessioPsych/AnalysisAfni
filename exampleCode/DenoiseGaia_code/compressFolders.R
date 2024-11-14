mainDir <- '/media/alessiof/Expansion/Gaia_David_info/Proj_Gaia_David'
runFlag <- 1
foldersToProcess <- c('EPI_timeSlicedCorrected_Safe','EPI_timeSlicedCorrected_Threat','freesurfer')

for ( k in 1:length( foldersToProcess ) ) {
  setwd( mainDir )
  setwd( foldersToProcess[k] )
  print( getwd() )
  subDirs <- list.dirs(recursive = FALSE, full.names = FALSE )
  for ( i in 1 : length( subDirs ) ) {
    if ( file.exists( sprintf('%s.tar.gz', subDirs[i] ) ) ) {
      instr <- sprintf('rm %s.tar.gz', subDirs[i] )
      print( instr )
      if (runFlag==1) { system( instr ) }
    }
    instr <- sprintf('tar -zcvf %s.tar.gz %s', subDirs[i], subDirs[i] )
    print( instr )
    if (runFlag==1) { system( instr ) }
  }
}

#to uncompress
#tar -zxvf sub03.tar.gz 