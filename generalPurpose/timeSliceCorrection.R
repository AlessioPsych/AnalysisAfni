args <- commandArgs(T)
print( args )

#to debug
#setwd('/media/alessiof/Data/tests/SEF_visual_response/HMR28')
#args <- c('EPI/','*.nii','EPI_jsons/','EPI_tsc/')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
library(rjson)

#arrange input
inputEpiDir <- args[1]
epiFormat <- args[2]
inputJsonDir <- args[3]
outputNameDir <- args[4]

print( sprintf('input EPI directory = %s', inputEpiDir ) )
print( sprintf('EPI file format = %s', epiFormat ) )
print( sprintf('input JSON directory = %s', inputJsonDir ) )
print( sprintf('output directory = %s', outputNameDir ) )

# create output directory
instr <- sprintf('mkdir %s', outputNameDir ); system( instr )

# get epi filenames
setwd( inputEpiDir )
filesWorkDir <- dir(pattern=epiFormat);
setwd( mainDir )

# get json filenames
setwd( inputJsonDir )
filesJsons <- dir(pattern='*.json');
setwd( mainDir )

# check that wehave the same number of epi and json files, abort otherwise
flagCheckNFiles <- length( filesJsons ) == length( filesWorkDir )
stopifnot( flagCheckNFiles )

# perform ts correction
for (k in 1:length(filesWorkDir)) {
  
  print( sprintf('processing epi file: %s', filesWorkDir[k] ) )
  
  print( sprintf('loading json file: %s', filesJsons[k] ) )
  jsonDset <- fromJSON(file = sprintf( '%s%s', inputJsonDir, filesJsons[k] ) )
  
  print( sprintf('writing slice timing for file: %s', filesJsons[k] ) )
  sliceTimingData <- jsonDset$SliceTiming
  write.table( x=t(sliceTimingData), file = sprintf( '%s_ttt_slice_timing_file.1D', inputEpiDir ), row.names = FALSE, col.names = FALSE  )
  
  instr <- sprintf('3dTshift -tpattern @%s_ttt_slice_timing_file.1D -prefix %s%s %s%s', inputEpiDir, outputNameDir, filesWorkDir[k], inputEpiDir, filesWorkDir[k]  ); system(instr)
  
  # clean up
  instr <- sprintf( 'rm %s_ttt_slice_timing_file.1D', inputEpiDir ); system( instr )
  
}
