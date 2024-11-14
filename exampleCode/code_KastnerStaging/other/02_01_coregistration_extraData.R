rm(list=ls())
epiDir <- '/analyse/Project0226/KastnerModel/staging_area_Kastner/epi_data_extra'
anatomyDir <- '/analyse/Project0226/KastnerModel/staging_area_Kastner/anatomies_KastnerClassic_Freesurfer'

setwd( epiDir )
print( sprintf('current folder: %s', getwd() ) )
epiSessions <- dir( epiDir )
print( sprintf('epi sessions:' ) )
epiSessions

setwd( anatomyDir )
print( sprintf('current folder: %s', getwd() ) )
anatomyFolders <- dir( anatomyDir, pattern='*_ANATOMY' )
print( sprintf('anatomy folders:' ) )
anatomyFolders

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

### align anatomyFolders to epiSession 
epiSessionSubj <- rep('aaa',length(epiSessions) )
for ( nSession in 1:length( epiSessions ) ) { # 
  currentEpiSession <- epiSessions[ nSession ]
  currentEpiSessionSplit <- strsplit( currentEpiSession, '_' )[[1]][1]
  epiSessionSubj[ nSession ] <- currentEpiSessionSplit
}

anatSubj <- rep('aaa',length(anatomyFolders) )
for ( nAnatomy in 1:length( anatomyFolders ) ) { # 
  currentAnatomy <- anatomyFolders[ nAnatomy ]
  currentAnatomySplit <- strsplit( currentAnatomy, '_' )[[1]][1]
  anatSubj[ nAnatomy ] <- currentAnatomySplit
}

anatIndex <- rep( 999, length( epiSessionSubj ) )
for( nSession in 1:length( epiSessionSubj ) ) { # nSession <- 1
  epiTemp <- epiSessionSubj[ nSession ]
  anatIndex[ nSession ] <- which( anatSubj == epiTemp )
}

#anatomyPerSession <- rep( anatomyFolders, array( table( epiSessionSubj ) ) )
anatomyPerSession <- anatomyFolders[ anatIndex ]

### check the alignment:
data.frame( anatomyPerSession, epiSessions )

for ( nSession in 1:9) { # length( epiSessions ) nSession <- 34; check session 20
  
  mainDir <- sprintf('%s/%s', epiDir, epiSessions[ nSession ] )
  anatomyDirLoop <- sprintf('%s/%s', anatomyDir, anatomyPerSession[ nSession ] )
  
  setwd( mainDir )
  getwd()
  print( getwd() )
  print( sprintf('session: %d', nSession ) )
  
  #### create coregister folder, if needed ####
  if ( dir.exists('coregister_01/') ) {
    instr <- 'rm -R coregister_01' 
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }  
  if ( !dir.exists('coregister_01/') ) {
    instr <- 'mkdir coregister_01' 
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }         
  
  ### copy anat
  instr <- sprintf( 'cp %s/FreeSeg_result/SUMA/anatCopy.nii.gz %s/coregister_01/anatCopy.nii.gz', anatomyDirLoop, mainDir )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  ### copy coreg
  coregDir01 <- list.dirs( full.names=FALSE, recursive=FALSE )
  coregDirFlag <- rep(999, length(coregDir01))
  for ( nDirs in 1:length( coregDir01) ) {
    coregFlag <- strfind( coregDir01[nDirs], 'coreg_' )
    if (is.null(coregFlag)) { coregDirFlag[ nDirs ] <- 0 } 
    if (!is.null(coregFlag)) { coregDirFlag[ nDirs ] <- 1 } 
  }
  coregDirIdx <- which( coregDirFlag==1 )[1]
  
  ### copy coregMat
  instr <- sprintf( 'cp %s/coregMat.1D %s/coregister_01/coregMat.1D', coregDir01[coregDirIdx], mainDir )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  ### copy invMat
  instr <- sprintf( 'cp %s/invMat.1D %s/coregister_01/invMat.1D', coregDir01[coregDirIdx], mainDir )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  instr <- 'computeAmplitudeAnatomy_nifti_gz.sh data_tsCorr_denoised_distCorr_motCorr/'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  instr <- 'mv amplitudeAnatomy.nii coregister_01/ttt_amplitudeAnatomy.nii'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  setwd('coregister_01')
  
    # clean up
    instr <- sprintf('rm singleShot.nii.gz')
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'cp ttt_amplitudeAnatomy.nii amplitudeAnatomy.nii' #to keep original starting point
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    #cp /analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Belinda/EPI/14_GN22NE438_Belinda_22082023_MB_GRE_EPI_005VOL_Brainstem_wholeBrain_62Slices_RL_2_2_2_20230822125110_14.nii /analyse/Project0226/tests/Daniela_pilotData/DanielaProject22082023_Belinda/coregistration/ttt_original_EPI_ts.nii
    
    epiOrigFiles <- dir('../EPI')
    instr <- sprintf('cp ../EPI/%s ttt_original_EPI_ts.nii', epiOrigFiles[1] )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # gets the oblique matrix and saves it on file as a 3X4 matrix
    instr <- '3dWarp -disp_obl_xform_only -deoblique ttt_original_EPI_ts.nii  > ttt_mat_obl_transform.1D'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # removes the comments from the oblique matrix file
    instr <- sprintf('grep -o \u0027^[^#]*\u0027 ttt_mat_obl_transform.1D > ttt_mat_obl_transform_no_comments.1D' )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # saves the oblique matrix file as a 1D file (as 3dAllineate requires)
    instr <- 'cat_matvec -ONELINE ttt_mat_obl_transform_no_comments.1D -I > ttt_mat_obl_transform_no_comments_ONELINE.1D'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'applyObliqueMatrix.sh ttt_amplitudeAnatomy.nii ttt_original_EPI_ts.nii'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <-' 3dWarp -deoblique -prefix ttt_amplitudeAnatomy_refit_deob.nii.gz ttt_amplitudeAnatomy_refit.nii.gz'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # apply transformation matrix to deoblique matrix, to check, here I am using the deoblique(d) space as master
    instr <- '3dAllineate -prefix ttt_original_EPI_mean_applyDeobMat.nii -1Dmatrix_apply ttt_mat_obl_transform_no_comments_ONELINE.1D -final wsinc5 -input ttt_amplitudeAnatomy.nii -master ttt_amplitudeAnatomy_refit_deob.nii.gz'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # align centers of EPI and anatCopy.nii.gz
    instr <- '@Align_Centers -cm -base anatCopy.nii.gz -dset ttt_original_EPI_mean_applyDeobMat.nii'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # upsample EPI to resolution of ANATOMY
    #3dresample -prefix _ttt_epi_startingPoint_resample.nii -input ttt_original_EPI_mean_applyDeobMat_shft.nii -rmode Linear -master anatCopy.nii.gz
    instr <- '3dresample -prefix _ttt_epi_startingPoint_resample.nii -input ttt_original_EPI_mean_applyDeobMat_shft.nii -rmode Linear -dxyz 1 1 1'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- '3dresample -prefix _ttt_anatCopy.nii.gz -input anatCopy.nii.gz -rmode Linear -master _ttt_epi_startingPoint_resample.nii'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # ---!!! UPDATE THIS LINE !!!---
    # apply nudge manually:
    #	- afni, underlay, ttt_amplitudeAnatomy.nii
    #	- new, _ttt_epi_startingPoint_resample.nii
    # 	- define datamode, plugins, nudge dataset, choose dataset, _ttt_epi_startingPoint_resample
    #	- print command line and copy-paste here
    instr <- sprintf( '3drotate -quintic -clipit %s -prefix nudgedDataset.nii.gz _ttt_epi_startingPoint_resample.nii', nudgeString )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # store affine transformation in 1D format
    instr <- sprintf('cat_matvec -ONELINE \u0027nudgedDataset.nii.gz::ROTATE_MATVEC_000000\u0027 -I > nudgeMat.1D')
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # perform coregistration
    instr <- 'align_epi_anat.py -epi_base 0 -epi nudgedDataset.nii.gz -anat _ttt_anatCopy.nii.gz -epi2anat -anat_has_skull no -epi_strip None -Allineate_opts -onepass -weight_frac 1.0 -maxrot 4 -maxshf 4 -VERB -warp aff'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # store all transformations in single transform matrix
    instr <- 'cat_matvec -ONELINE nudgedDataset_al_mat.aff12.1D nudgeMat.1D ttt_original_EPI_mean_applyDeobMat_shft.1D ttt_mat_obl_transform_no_comments_ONELINE.1D > coregMat.1D' 
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'cat_matvec -ONELINE coregMat.1D -I > invMat.1D'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # apply single transform matrix to initial input
    instr <- '3dAllineate -prefix singleShot.nii.gz -1Dmatrix_apply coregMat.1D -final linear -input ttt_amplitudeAnatomy.nii -master anatCopy.nii.gz'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    # clean up
    instr <- 'rm ttt*'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'rm _ttt*'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'rm nudgedDataset.nii.gz'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'rm nudgeMat.1D'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'rm anatCopy_al_mat.aff12.1D'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
    
    instr <- 'rm nudgedDataset_*'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  
}


