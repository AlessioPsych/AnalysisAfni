args <- commandArgs(T)
print( args )

#to debug
#setwd('/analyse/Project0226/tests/DanielaProject')
#args <- c('EPI/','*.nii','4','EPI_denoised/')

mainDir <- getwd()
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 

#arrange input
inputEpiDir <- args[1]
epiFormat <- args[2]
ncores <- args[3]
outputNameDir <- args[4]

print( sprintf('input EPI directory = %s', inputEpiDir ) )
print( sprintf('EPI file format = %s', epiFormat ) )
print( sprintf('ncores = %s', ncores ) )
print( sprintf('output directory = %s', outputNameDir ) )

# create output directory
instr <- sprintf('mkdir %s', outputNameDir ); system( instr )

# get epi filenames
setwd( inputEpiDir )
filesWorkDir <- dir(pattern=epiFormat);
setwd( mainDir )

# perform data denoising correction
for (k in 1:length(filesWorkDir)) {
  
  print( sprintf('processing epi file: %s', filesWorkDir[k] ) )
  
  instr <- sprintf('dwidenoise -nthreads %s %s/%s%s %s/%s%s ', ncores, mainDir, inputEpiDir, filesWorkDir[k], mainDir, outputNameDir, filesWorkDir[k]  ); system(instr)
  
}
