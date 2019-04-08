args <- commandArgs(T)
print( args )
#setwd('/analyse/Project0226/KMA25')

#rm(list=ls())
#args <- c( 'bars_topUp/', '*.BRIK', '3' )

mainDir <- getwd()
workDir <- strsplit( args[1], '[/]' )[[1]]
fileFormat <- args[2]
detrendOrder <- args[3]

setwd( mainDir )
instr <- sprintf('mkdir %s_detrend/', workDir )
system( instr )

setwd( mainDir )
setwd( workDir )
filesWorkDir <- dir(pattern=fileFormat);
setwd( mainDir )

for (k in 1:length(filesWorkDir)) {
  instr <- sprintf('3dTstat -prefix %s_detrend/_ttt_meanTs.nii.gz %s/%s', workDir, workDir, filesWorkDir[k] )
  system( instr )
  instr <- sprintf( '3dcalc -a %s/%s -b %s_detrend/_ttt_meanTs.nii.gz -expr \u0027 min(200, a/b*100)*step(a)*step(b) \u0027 -prefix %s_detrend/_ttt_scaled.nii.gz' , workDir, filesWorkDir[k], workDir, workDir, filesWorkDir[k] )
  system( instr )
  instr <- sprintf( '3dDetrend -polort %s -prefix %s_detrend/%s %s_detrend/_ttt_scaled.nii.gz', args[3], workDir, filesWorkDir[k], workDir )
  print( instr )
  system( instr )
  instr <- sprintf('rm %s_detrend/_ttt*', workDir)
  system( instr )
}
