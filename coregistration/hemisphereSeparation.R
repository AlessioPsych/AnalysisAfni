args <- commandArgs(T)
print( args )

#args <- c('anatomy_res.nii.gz','segmentation47.nii.gz','/usr/share/afni/atlases/TT_icbm452+tlrc','1')

inputAnat <- args[1]
inputSegmentation <- args[2]
atlasPath <- args[3]
flagLin <- args[4]

instr <- sprintf('conformToTlrc.sh %s %s %s', inputAnat, atlasPath, flagLin)
print( instr )
system( instr )

if (flagLin==1) {
  instr <- sprintf('3dcalc -a Qwarp+tlrc -expr \u0027within(x,1,1000)*1\u0027 -prefix ttt_roi.nii.gz')
  print( instr )
  system( instr )
  instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s ttt_roi.nii.gz 1 Qwarp_WARPINV+tlrc', inputAnat)
  print( instr )
  system( instr )
}
if (flagLin!=1) {
  instr <- sprintf('3dcalc -a MPRAGE_at.nii -expr \u0027within(x,1,1000)*1\u0027 -prefix ttt_roi.nii.gz')
  print( instr )
  system( instr )
  instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s ttt_roi.nii.gz 0', inputAnat)
  print( instr )
  system( instr )
}

instr <- sprintf('3dmerge -1blur_fwhm 1.5 -prefix ttt_maskBlur.nii.gz outputVolume.nii.gz')
print( instr )
system( instr )
instr <- sprintf('3dcalc -a ttt_maskBlur.nii.gz -b %s -expr \u0027within(a,0.5,1)*b\u0027 -prefix leftSeg.nii.gz', inputSegmentation)
print( instr )
system( instr )
instr <- sprintf('3dcalc -a ttt_maskBlur.nii.gz -b %s -expr \u0027within(a,0,0.499)*b\u0027 -prefix rightSeg.nii.gz', inputSegmentation)
print( instr )
system( instr )

system('rm ttt_*')
