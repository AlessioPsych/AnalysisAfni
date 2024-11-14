rm( list=ls() ); gc();

#need to load module R, afni, and activate conda environment beforehand

#### Move to the desired folder ####
sprintf('Move to the desired folder:')
mainDir <- '/scratch/af4887/Proj_Gaia_David'
dataDirSafe <- '/scratch/af4887/Proj_Gaia_David/afni_processed_Safe'
dataDirThreat <- '/scratch/af4887/Proj_Gaia_David/afni_processed_Threat'
inputFreesurfer <- '/scratch/af4887/Proj_Gaia_David/freesurfer'
setwd( mainDir )
sprintf( 'Current folder: %s', print( getwd() ) )
singleSubjectFolders_Safe <- list.dirs( dataDirSafe, recursive = FALSE )
singleSubjectFolders_Threat <- list.dirs( dataDirThreat, recursive = FALSE )
singleSubjectFolders_Freesurfer <- list.dirs( inputFreesurfer, recursive=FALSE )
runCodeFlag <- 1

for ( nSubj in 1 : 30  ) { # nSubj <- 1 length( singleSubjectFolders )
  
  setwd( mainDir )
  print( getwd() )

  # define input dataset folder
  dsetsFolder <- sprintf('%s', singleSubjectFolders_Safe[ nSubj ] )
  setwd( dsetsFolder )
  print( '###################' )
  print( '###################' )
  print( '###################' )
  print( getwd() )
  print( '###################' )
  print( '###################' )
  print( '###################' )

  currentSubjFolder <- strsplit( singleSubjectFolders_Safe[ nSubj ], '[/]' )[[1]][6]
  currentSubj <- strsplit( currentSubjFolder, '[.]' )[[1]][1]
  
  # stats safe into SUMA
  if ( file.exists( sprintf('%s/SUMA/stats_safe.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats_safe.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  if ( file.exists( sprintf('%s/SUMA/vr_base_min_outlier_safe.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/vr_base_min_outlier_safe.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  
  instr <- paste('afni 3dAllineate',
                 sprintf('-base %s/SUMA/brain.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 sprintf('-input %s/stats.%s+orig', singleSubjectFolders_Safe[ nSubj ], currentSubj ),
                 sprintf('-1Dmatrix_apply %s/Coregistration/vr_base_min_outlier_al_mat.aff12.1D', singleSubjectFolders_Safe[ nSubj ] ),
                 sprintf('-prefix %s/SUMA/stats_safe.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 '-final linear',
                 '-mast_dxyz 2'
                 )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  instr <- paste('afni 3dAllineate',
                 sprintf('-base %s/SUMA/brain.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 sprintf('-input %s/Coregistration/vr_base_min_outlier+orig', singleSubjectFolders_Safe[ nSubj ] ),
                 sprintf('-1Dmatrix_apply %s/Coregistration/vr_base_min_outlier_al_mat.aff12.1D', singleSubjectFolders_Safe[ nSubj ] ),
                 sprintf('-prefix %s/SUMA/vr_base_min_outlier_safe.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 '-final linear',
                 '-mast_dxyz 2'
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # stats threat info SUMA
  if ( file.exists( sprintf('%s/SUMA/stats_threat.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats_threat.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  if ( file.exists( sprintf('%s/SUMA/vr_base_min_outlier_threat.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/vr_base_min_outlier_threat.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  
  instr <- paste('afni 3dAllineate',
                 sprintf('-base %s/SUMA/brain.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 sprintf('-input %s/stats.%s+orig', singleSubjectFolders_Threat[ nSubj ], currentSubj ),
                 sprintf('-1Dmatrix_apply %s/Coregistration/vr_base_min_outlier_al_mat.aff12.1D', singleSubjectFolders_Threat[ nSubj ] ),
                 sprintf('-prefix %s/SUMA/stats_threat.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 '-final linear',
                 '-mast_dxyz 2'
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  instr <- paste('afni 3dAllineate',
                 sprintf('-base %s/SUMA/brain.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 sprintf('-input %s/vr_base_min_outlier+orig', singleSubjectFolders_Threat[ nSubj ] ),
                 sprintf('-1Dmatrix_apply %s/Coregistration/vr_base_min_outlier_al_mat.aff12.1D', singleSubjectFolders_Threat[ nSubj ] ),
                 sprintf('-prefix %s/SUMA/vr_base_min_outlier_threat.nii.gz', singleSubjectFolders_Freesurfer[ nSubj ] ),
                 '-final linear',
                 '-mast_dxyz 2'
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  setwd( singleSubjectFolders_Freesurfer[ nSubj ] )
  setwd('SUMA')
  print( getwd() )
  
  # surfaces safe lh niml
  if ( file.exists( sprintf('%s/SUMA/stats.safe.lh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.safe.lh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.lh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.lh.pial.gii' ),
                 sprintf('-sv stats_safe.nii.gz' ),
                 sprintf('-grid_parent stats_safe.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_niml stats.safe.lh.niml.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # surfaces safe rh niml
  if ( file.exists( sprintf('%s/SUMA/stats.safe.rh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.safe.rh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.rh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.rh.pial.gii' ),
                 sprintf('-sv stats_safe.nii.gz' ),
                 sprintf('-grid_parent stats_safe.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_niml stats.safe.rh.niml.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # surfaces threat lh niml
  if ( file.exists( sprintf('%s/SUMA/stats.threat.lh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.threat.lh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.lh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.lh.pial.gii' ),
                 sprintf('-sv stats_threat.nii.gz' ),
                 sprintf('-grid_parent stats_threat.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_niml stats.threat.lh.niml.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # surfaces threat rh niml
  if ( file.exists( sprintf('%s/SUMA/stats.threat.rh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.threat.rh.niml.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.rh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.rh.pial.gii' ),
                 sprintf('-sv stats_safe.nii.gz' ),
                 sprintf('-grid_parent stats_threat.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_niml stats.threat.rh.niml.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  # surfaces safe lh 1D
  if ( file.exists( sprintf('%s/SUMA/stats.safe.lh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.safe.lh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.lh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.lh.pial.gii' ),
                 sprintf('-sv stats_safe.nii.gz' ),
                 sprintf('-grid_parent stats_safe.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_1D stats.safe.lh.1D.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # surfaces safe rh 1D
  if ( file.exists( sprintf('%s/SUMA/stats.safe.rh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.safe.rh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.rh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.rh.pial.gii' ),
                 sprintf('-sv stats_safe.nii.gz' ),
                 sprintf('-grid_parent stats_safe.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_1D stats.safe.rh.1D.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # surfaces threat lh 1D
  if ( file.exists( sprintf('%s/SUMA/stats.threat.lh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.threat.lh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.lh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.lh.pial.gii' ),
                 sprintf('-sv stats_threat.nii.gz' ),
                 sprintf('-grid_parent stats_threat.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_1D stats.threat.lh.1D.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # surfaces safe rh 1D
  if ( file.exists( sprintf('%s/SUMA/stats.threat.rh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] ) ) ) {
    instr <- sprintf('rm %s/SUMA/stats.threat.rh.1D.dset', singleSubjectFolders_Freesurfer[ nSubj ] )
    if (runCodeFlag==1) { system( instr ) }
  }
  instr <- paste('afni 3dVol2Surf',
                 sprintf('-spec std.60.FreeSeg_results_both.spec' ),
                 sprintf('-surf_A std.60.rh.smoothwm.gii' ),
                 sprintf('-surf_B std.60.rh.pial.gii' ),
                 sprintf('-sv stats_safe.nii.gz' ),
                 sprintf('-grid_parent stats_threat.nii.gz' ),
                 sprintf('-map_func ave'),
                 sprintf('-f_steps 10'),
                 sprintf('-f_index nodes'),
                 sprintf('-out_1D stats.threat.rh.1D.dset')                 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
} 

