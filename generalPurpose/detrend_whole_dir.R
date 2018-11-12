#args <- commandArgs(T)
#print( args )

#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V3813Alessio_copy/prfModel')
#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V3814Serge')
#setwd('/home/fracasso/data/backupStorage02/storage2/prfSizeDepth_2/V4162Ben')
#rm(list=ls())
#setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

#args <- c('meanTsEye_res_mask.nii.gz','meanTsEye_res_detrend.nii.gz','output_eye_sampling_all_positive','/home/alessiof/abin','directionFile.1D' ,'/home/alessiof/Programs/dlDevelopment/analysisAfni/generalPurpose')


#args <- commandArgs(T)
#print( args )
setwd('/analyse/Project0226/GN18NE278_KMA25_FEF_28092018_nifti')

rm(list=ls())
args <- c( 'epiMotionCorrectEye_topUp/', '*.BRIK', '4' )

mainDir <- getwd()
workDir <- strsplit( args[1], '[/]' )[[1]]
fileFormat <- args[2]
detrendOrder <- args[3]
setwd( mainDir )
instr <- sprintf('mkdir %s_detrend/', workDir )
system( instr )
setwd( workDir )
filesWorkDir <- dir(pattern=fileFormat);
setwd( mainDir )

for (k in 1:length(filesWorkDir)) {
  instr <- sprintf('bet %s _ttt_mask.nii -m -f %s', filesCpDir[k], args[2] )
  print( instr )
  system( instr )
  if (k < 10) { instr <- sprintf('3dcalc -a _ttt_mask.nii -b %s -expr \u0027step(a)*b\u0027 -prefix 0%1.0f.nii.gz', filesCpDir[k], k) }
  if (k >= 10) { instr <- sprintf('3dcalc -a _ttt_mask.nii -b %s -expr \u0027step(a)*b\u0027 -prefix %1.0f.nii.gz', filesCpDir[k], k) }
  print( instr )
  system( instr );
}
system('rm _ttt*')
system('rm *.nii')
system('gunzip *.gz')
setwd( mainDir )
