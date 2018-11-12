args <- commandArgs(T)
print( args )

#args <- c('_ttt_T1MPRAGE_SS_N4.nii.gz', '/home/alessiof/abin', '0')
#args <- c('pdCorrectRegularT1_noBlur_stripped.nii','/packages/afni/17.0.13','1')
#'TT_152_2009c+tlrc'
#setwd('/home/dijkj/Documents/Segmentation_tests/JeDi/RemoveCerebellum')
#setwd('/home/alessiof/Mount/Data/anatData/anatAnalysis_06032018/20171101_VBE10_preproc')

inputAnat <- args[1]
atlasPath <- args[2]
flagLin <- args[3]

instr <- sprintf('conformToTlrc.sh %s %s/TT_152_2009c+tlrc %s', inputAnat, atlasPath, flagLin)
print( instr )
system( instr )

#now a step using the TT_desai_dk_mpm+tlrc roi definitions, 3dCalc. 
#1,74,75,76,81,82,86,92,93,94

if (flagLin==1) {
  # first for the full anatomy
  instr <- sprintf('3dcalc -a %s/TT_desai_dk_mpm+tlrc -expr \u0027within(a,74,75)+within(a,92,93)+within(a,81,82)\u0027 -prefix ttt_roi.nii.gz', atlasPath, atlasPath)
  print( instr )
  system( instr )
  instr <- sprintf('3dmask_tool -input ttt_roi.nii.gz -prefix ttt_roi_dilate.nii.gz -dilate_input 1') #dilate the roi a bit
  print( instr )
  system( instr )
  instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s ttt_roi_dilate.nii.gz 1 Qwarp_WARPINV+tlrc', inputAnat)
  print( instr ) #warp anatomy back to original space
  system( instr )
  instr <- sprintf('3dcopy outputVolume.nii.gz ttt_fullMask.nii.gz')
  print( instr )
  system( instr )
  system('rm outputVolume.nii.gz')
  
  #Only for one hemisphere
  #instr <- sprintf('3dcalc -a ttt_roi.nii.gz -expr \u0027within(x,1,1000)*1\u0027 -prefix ttt_roiLeft.nii.gz') #x is ta space coordinate
  #print( instr ) #get only the left hemisphere roi
  #system( instr )
  #instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s ttt_roiLeft.nii.gz 1 Qwarp_WARPINV+tlrc', inputAnat)
  #print( instr ) #warp roi back to original space
  #system( instr )
  #instr <- sprintf('3dcopy outputVolume.nii.gz leftMask.nii.gz')
  #print( instr )
  #system( instr )
  #system('rm outputVolume.nii.gz')

}
if (flagLin!=1) {
  # first for the full anatomy
  instr <- sprintf('3dcalc -a %s/TT_desai_dk_mpm+tlrc -expr \u0027within(a,74,75)+within(a,92,93)+within(a,81,82)\u0027 -prefix ttt_roi.nii.gz', atlasPath, atlasPath)
  print( instr )
  system( instr )
  instr <- sprintf('3dmask_tool -input ttt_roi.nii.gz -prefix ttt_roi_dilate.nii.gz -dilate_input 1') #dilate the roi a bit
  print( instr )
  system( instr )
  instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s ttt_roi_dilate.nii.gz 0', inputAnat)
  print( instr ) #warp anatomy back to original space
  system( instr )
  instr <- sprintf('3dcopy outputVolume.nii.gz ttt_fullMask.nii.gz')
  print( instr )
  system( instr )
  system('rm outputVolume.nii.gz')
  
  #Only for one hemisphere
  #instr <- sprintf('3dcalc -a ttt_roi.nii.gz -expr \u0027within(x,1,1000)*1\u0027 -prefix ttt_roiLeft.nii.gz') #x is ta space coordinate
  #print( instr ) #get only the left hemisphere roi
  #system( instr )
  #instr <- sprintf('convertRoi_from_desi_atlas.sh MPRAGE_at.Xat.1D %s ttt_roiLeft.nii.gz 0', inputAnat)
  #print( instr ) #warp roi back to original space
  #system( instr )
  #instr <- sprintf('3dcopy outputVolume.nii.gz leftMask.nii.gz')
  #print( instr )
  #system( instr )
  #system('rm outputVolume.nii.gz')
}

instr <- '3dclust 0 20 ttt_fullMask.nii.gz > out.1D'
system( instr )
clustTable <- read.table( 'out.1D', comment.char = "#" )
system( 'rm out.1D' )
clusteringInst <- sprintf('3dclust -prefix fullMask_clust.nii.gz 0 %1.0f ttt_fullMask.nii.gz', clustTable[1,1] - 1 )
system(clusteringInst)

instr <- sprintf('3dcalc -a %s -b fullMask_clust.nii.gz -expr \u0027a*not(b)\u0027 -prefix anatomy_noCB.nii.gz', inputAnat)
print( instr ) #get the original anatomy without cerebellum
system( instr )

system('rm MPRAGE_at.nii')
system('rm MPRAGE_at.nii.Xaff12.1D')
system('rm MPRAGE_at.nii_WarpDrive.log')
system('rm MPRAGE_at.Xat.1D')
system('rm fullMask_clust.nii.gz')
system('rm ttt_fullMask.nii.gz')
system('rm ttt_roi_dilate.nii.gz')
system('rm ttt_roi.nii.gz')


