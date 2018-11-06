args <- commandArgs(T)
print( args )

if (args[2]!=0) {
  
  instr <- sprintf('3dBrickStat -percentile %s %s %s %s > temp.1D', args[2], args[2], args[2], args[1])
  print( instr )
  system( instr )
  
  percentileData <- read.table( file='temp.1D', as.is=TRUE )
  
  instr <- sprintf('3dcalc -a %s -expr \u0027a*step(a-%1.2f)\u0027 -prefix tempFile_calcMask+orig', args[1], percentileData[2] )
  print( instr )
  system( instr )
  
  instr <- sprintf('3dAutomask -apply_prefix tempFile_calcMask_auto+orig -clfrac %s -dilate %s -erode %s tempFile_calcMask+orig', args[3], args[4], args[5])
  print( instr )
  system( instr )
    
}

if (args[2]==0) {
  instr <- sprintf('3dcopy %s tempFile_calcMask_auto+orig', args[1])
  print( instr )
  system( instr )
}

instr <- sprintf('3dZeropad -R %s -L %s -A %s -P %s -I %s -S %s -prefix tempFile_calcMask_auto_zp+orig tempFile_calcMask_auto+orig', args[6], args[6], args[6], args[6], args[6], args[6] )
print( instr )
system( instr )

if (args[9]==1) {  
  instr <- sprintf('3dAttribute IJK_TO_DICOM_REAL tempFile_calcMask_auto+orig > temp_obl_matrix.1D')
  print( instr )
  system( instr )

  instr <- sprintf('3drefit -atrfloat IJK_TO_DICOM_REAL temp_obl_matrix.1D tempFile_calcMask_auto_zp+orig')
  print( instr )
  system( instr )
}

if (args[8]==0) {
  instr <- sprintf('3dcopy tempFile_calcMask_auto_zp+orig tempFile_calcMask_auto_zp_unif+orig' )
}
if (args[8]==1) {
  instr <- sprintf('3dUnifize -input tempFile_calcMask_auto_zp+orig -prefix tempFile_calcMask_auto_zp_unif+orig' )
}
if (args[8]==2) {
  instr <- sprintf('3dUnifize -input tempFile_calcMask_auto_zp+orig -prefix tempFile_calcMask_auto_zp_unif+orig -T2 -T2' )
}
if (args[8]==3) {
  instr <- sprintf('3dUniformize -anat tempFile_calcMask_auto_zp+orig -prefix tempFile_calcMask_auto_zp_unif+orig' )
}

print( instr )
system( instr )

instr <- sprintf('3dcopy tempFile_calcMask_auto_zp_unif+orig %s', args[7])
print( instr )
system( instr )

#instr <- sprintf('3dAFNItoNIFTI -prefix %s tempFile_calcMask_auto_zp_unif+orig', args[7] )
#print( instr )
#system( instr )
