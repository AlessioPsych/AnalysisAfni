
mainDir <- '/analyse/Project0370/20230620_CCY20/ANATOMY'
setwd( mainDir )
getwd()

instr <- 'skullStrip_mp2rage.sh 42_20230620_CCY20_mp2rage_sag_p3_0.6mm_20230620094607_42_c32_e1.nii 41_20230620_CCY20_mp2rage_sag_p3_0.6mm_20230620094607_41_c32.nii  0.7'; system( instr )

  runFlag <- 1
  nCores <- 8
  fileAnatomy <- 'anatCopy.nii.gz'
  #redefine SUBJECTS_DIR to current directory
  currentDir <- getwd()
  print( currentDir )
  print( '...' )
  Sys.setenv( SUBJECTS_DIR = paste( currentDir ) )
  print( Sys.getenv( 'SUBJECTS_DIR') )
  print( '...' )
  
  #recon-all instruction
  print( '...' )
  instr <- sprintf( 'recon-all -subjid FreeSeg_result -i %s -all -parallel -openmp %1.0f', fileAnatomy, nCores )
  print( instr )
  if ( runFlag==1 ) {
    system( instr )
  }
  print( '...' )

  #hippocampal subfields segmentation
  instr <- 'recon-all -s FreeSeg_result -hippocampal-subfields-T1'; system( instr )
  
  #suma instruction
  print( '...' )
  instr <- sprintf( '@SUMA_Make_Spec_FS -NIFTI -fspath %s/FreeSeg_result -sid FreeSeg_result', currentDir )
  print( instr )
  if ( runFlag==1 ) {
    system( instr )
  }
  print( '...' )

  subjectDir <- Sys.getenv( 'SUBJECTS_DIR')
  freeSurfDir <- sprintf( '%s/FreeSeg_result/mri/', subjectDir )
  sumaDir <- sprintf( '%s/FreeSeg_result/SUMA/', subjectDir )

  instr <- sprintf('mri_convert %s/rh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz %s/rh.hippoSfLabels-T1.v10.FSvoxelSpace.nii', freeSurfDir, sumaDir ); system( instr )

  instr <- sprintf('mri_convert %s/lh.hippoSfLabels-T1.v10.FSvoxelSpace.mgz %s/lh.hippoSfLabels-T1.v10.FSvoxelSpace.nii', freeSurfDir, sumaDir ); system( instr )

  instr <- sprintf('mri_convert %s/rh.hippoSfLabels-T1.v10.mgz %s/rh.hippoSfLabels-T1.v10.nii', freeSurfDir, sumaDir ); system( instr )

  instr <- sprintf('mri_convert %s/lh.hippoSfLabels-T1.v10.mgz %s/lh.hippoSfLabels-T1.v10.nii', freeSurfDir, sumaDir ); system( instr )

  instr <- sprintf('cp %s/rh.hippoSfVolumes-T1.v10.txt %s/rh.hippoSfVolumes-T1.v10.txt', freeSurfDir, sumaDir ); system( instr )

  instr <- sprintf('cp %s/lh.hippoSfVolumes-T1.v10.txt %s/lh.hippoSfVolumes-T1.v10.txt', freeSurfDir, sumaDir ); system( instr )


instr <- sprintf('cp %s/anatCopy.nii.gz %s/anatCopy.nii.gz', mainDir, sumaDir ); system( instr )

setwd(sumaDir)
instr <- sprintf('3dresample -master anatCopy.nii.gz -input rh.hippoSfLabels-T1.v10.FSvoxelSpace.nii -rmode NN -prefix rh.hippoSfLabels-T1.v10.FSvoxelSpace_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input lh.hippoSfLabels-T1.v10.FSvoxelSpace.nii -rmode NN -prefix lh.hippoSfLabels-T1.v10.FSvoxelSpace_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input rh.hippoSfLabels-T1.v10.nii -rmode NN -prefix rh.hippoSfLabels-T1.v10_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input lh.hippoSfLabels-T1.v10.nii -rmode NN -prefix lh.hippoSfLabels-T1.v10_resampled.nii'); system( instr )

setwd(sumaDir)
instr <- sprintf('3dresample -master anatCopy.nii.gz -input rh.hippoSfLabels-T1.v10.FSvoxelSpace.nii -rmode NN -prefix rh.hippoSfLabels-T1.v10.FSvoxelSpace_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input lh.hippoSfLabels-T1.v10.FSvoxelSpace.nii -rmode NN -prefix lh.hippoSfLabels-T1.v10.FSvoxelSpace_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input rh.hippoSfLabels-T1.v10.nii -rmode NN -prefix rh.hippoSfLabels-T1.v10_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input lh.hippoSfLabels-T1.v10.nii -rmode NN -prefix lh.hippoSfLabels-T1.v10_resampled.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input aparc.a2009s+aseg.nii -rmode NN -prefix aparc.a2009s+aseg_resample.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input aparc.a2009s+aseg_REN_all.nii.gz -rmode NN -prefix aparc.a2009s+aseg_REN_all_resample.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input aparc+aseg.nii -rmode NN -prefix aparc+aseg_resample.nii'); system( instr )
instr <- sprintf('3dresample -master anatCopy.nii.gz -input aparc+aseg_REN_all.nii.gz -rmode NN -prefix aparc+aseg_REN_all_resample.nii'); system( instr )




