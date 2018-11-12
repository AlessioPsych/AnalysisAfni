args <- commandArgs(T)
print( args )

#args <- ('TOPUP/')
mainDir <- getwd()
cpDir <- strsplit( args[1], '[/]' )[[1]]
setwd( mainDir )
instr <- sprintf('cp -R %s %s_mask/', args[1], cpDir )
system( instr )
setwd( sprintf( '%s_mask', cpDir ) )
filesCpDir <- dir(pattern='*.nii');
for (k in 1:length(filesCpDir)) {
  if (k < 10) { instr <- sprintf('3dAutomask -clfrac %s -apply_prefix 0%1.0f.nii.gz %s', args[2], k, filesCpDir[k] ) }
  if (k >= 10) { instr <- sprintf('3dAutomask -clfrac %s -apply_prefix %1.0f.nii.gz %s', args[2], k, filesCpDir[k] ) }
  #instr <- sprintf('bet %s _ttt_mask.nii -m -f %s', filesCpDir[k], args[2] )
  #print( instr )
  #system( instr )
  #if (k < 10) { instr <- sprintf('3dcalc -a _ttt_mask.nii -b %s -expr \u0027step(a)*b\u0027 -prefix 0%1.0f.nii.gz', filesCpDir[k], k) }
  #if (k >= 10) { instr <- sprintf('3dcalc -a _ttt_mask.nii -b %s -expr \u0027step(a)*b\u0027 -prefix %1.0f.nii.gz', filesCpDir[k], k) }
  print( instr )
  system( instr )
}
#system('rm _ttt*')
system('rm *.nii')
system('gunzip *.gz')
setwd( mainDir )
