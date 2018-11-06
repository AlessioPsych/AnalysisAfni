args <- commandArgs(T)
print( args )

#setwd('/data2/projects/Exchange/UWAMB/150702a')
#args <- c('motionCorrect.results','niftiFiles')

inputDir <- args[1]
outputDir <- args[2]

mainDir <- getwd()

targetDir <- outputDir
flagDir <- dir.create( file.path(mainDir, targetDir) )
if (flagDir==FALSE) { # directory already exists!
  msg <- sprintf( 'Remove the directory %s to proceed', targetDir )
  warning( msg )
  stopifnot(flagDir)  
}

setwd( inputDir )
filesAll <- list.files()
filesDir <- filesAll[ grep('volreg+orig.BRIK',filesAll,fixed=T) ]
for (k in 1:length(filesDir)) {
  if (k < 10) {
    filenameOut <- sprintf('0%s.nii',k)
  }
  if (k >= 10) {
    filenameOut <- sprintf('%s.nii',k)
  }
  instr <- sprintf('3dAFNItoNIFTI -prefix ../%s/%s %s', targetDir, filenameOut, filesDir[k] )
  print(instr)
  system( instr )
}
setwd( mainDir )
