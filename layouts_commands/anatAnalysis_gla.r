rm( list=ls() )
library('ANTsR')
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniAtlasDir <- Sys.getenv(x='AFNI_ATLASDIR')
source( sprintf( '%s/AFNIio.R', afniDir ) )

mainDir <- '/home/fracasso/Downloads/anatData/4Fracasso' 
outputDir <- '/home/fracasso/Downloads/anatData/dataOut' 
setwd( mainDir )
instr <- sprintf( 'mkdir %s/store_preproc/', outputDir )
system( instr )
inDir <- list.dirs( mainDir, recursive = FALSE, full.names = FALSE )
for ( nDir in 1:length( inDir ) ) {

  instr <- sprintf( 'mkdir %s/%s_preproc/', outputDir, inDir[ nDir ] )
  setwd( mainDir )
  system( instr )
  setwd( inDir[ nDir ] )
  setwd('ANAT/T13DMPRAGE')
  instr <- sprintf( 'dcm2nii -o %s/%s_preproc/ *.IMA', outputDir, inDir[ nDir ] )
  system( instr )
  
  setwd( mainDir )
  setwd( inDir[ nDir ] )
  setwd('LOCALISER')
  instr <- sprintf( 'dcm2nii -o %s/%s_preproc/ *.IMA', outputDir, inDir[ nDir ] )
  system( instr )
  
  setwd( sprintf( '%s/%s_preproc/', outputDir, inDir[ nDir ] ) )
  fileToProcess <- dir( pattern='co2' )
  instr <- sprintf('3dresample -orient RAI -prefix _ttt_T1_orient.nii.gz -inset %s', fileToProcess)
  system( instr )
  instr <- sprintf('3dUnifize -prefix _ttt_T1_orient_uni.nii.gz _ttt_T1_orient.nii.gz')
  system( instr )
  instr <- sprintf('3dSkullStrip -input _ttt_T1_orient_uni.nii.gz -blur_fwhm 2 -prefix _ttt_T1_orient_uni_ss.nii.gz')
  system( instr )
  instr <- sprintf('3dcalc -a _ttt_T1_orient_uni_ss.nii.gz -expr \u0027step(a)\u0027 -prefix _ttt_T1_orient_uni_ss_mask.nii.gz')
  system( instr )
  instr <- sprintf('3dmask_tool -input _ttt_T1_orient_uni_ss_mask.nii.gz -prefix _ttt_T1_orient_uni_ss_mask_dil.nii.gz -dilate_input 2')
  system( instr )
  instr <- sprintf('3dcalc -a _ttt_T1_orient_uni_ss_mask_dil.nii.gz -b _ttt_T1_orient.nii.gz -expr \u0027step(a)*b\u0027 -prefix _ttt_T1MPRAGE_SS.nii.gz')
  system( instr )
  
  instr <- sprintf('cp _ttt_T1MPRAGE_SS.nii.gz %s/store_preproc/%s_T1_MPRAGE.nii.gz',outputDir, inDir[ nDir ] )  
  system( instr )
  
  #instr <- sprintf('3dresample -dxyz 0.6 0.6 0.6 -inset _ttt_T1MPRAGE_SS.nii.gz -prefix _ttt_T1MPRAGE_SS_res.nii.gz -rmode Lin')
  #system( instr )
    
  #system('segment_white_matter_params.sh _ttt_T1MPRAGE_SS_res.nii.gz 0.000001 0.00001 5 8-8-8-8-9-9 0')
  #system('segment_gray_matter.sh _ttt_T1MPRAGE_SS_res.nii.gz white_matter_mask.nii.gz 0.65-0.65-0.65-0.65')

  #system('segment_white_matter_params.sh _ttt_T1MPRAGE_SS_res.nii.gz 0.0000001 0.000001 5 8-8-8-8-9-9 0 1 1')
  

  #system('segment_white_matter_params.sh _ttt_T1MPRAGE_SS_N4.nii.gz 0.001 0.001 5 3-3-3-3-3-3 0 1 1')
  #system('segment_gray_matter.sh _ttt_T1MPRAGE_SS_N4.nii.gz white_matter_mask_N4.nii.gz 0.5-0.5-0.5-0.5')
  
  
  volumeAnts <- antsImageRead( '_ttt_T1MPRAGE_SS.nii.gz' ); #here I use ANTsR, modify
  volumeAntsN4 <- n4BiasFieldCorrection( volumeAnts, 5 )
  antsImageWrite( volumeAntsN4, '_ttt_T1MPRAGE_SS_N4.nii.gz' ) 

  system('segment_white_matter_params.sh _ttt_T1MPRAGE_SS_N4.nii.gz 0.000001 0.000001 5 7-7-7-7-7-7 0 0 1')
  system('mv white_matter_mask.nii.gz white_matter_mask_N4.nii.gz')
    
  instr <- sprintf('removeCerebellum.sh _ttt_T1MPRAGE_SS_N4.nii.gz %s 1', afniAtlasDir)
  system( instr )
  
  instr <- sprintf('3dcalc -a anatomy_noCB.nii.gz -expr \u0027step(a)\u027 -prefix mask_no_CB.nii.gz')
  system( instr )
  instr <- sprintf('3dSeg -anat anatomy_noCB.nii.gz -mask mask_no_CB.nii.gz -classes \u0027CSF ; GM ; WM\u0027 -bias_classes \u0027GM ; WM\u0027 -bias_fwhm 0.0 -mixfrac UNI -main_N 15 -blur_meth BFT' )
  system( instr )
  
  setwd('Segsy')
  system('cp Classes+orig* ../')
  setwd( sprintf( '%s/%s_preproc/', outputDir, inDir[ nDir ] ) )
  system('3dAFNItoNIFTI -prefix classes.nii.gz Classes+orig.BRIK')
  system('3dcalc -a classes.nii.gz -expr \u0027within(a,2,3)\u0027 -prefix brainmask.nii.gz')
  #system('3dcalc -a classes.nii.gz -expr \u0027within(a,3,3)\u0027 -prefix wmmask.nii.gz')
  
  instr <- sprintf('3dclust 0 10 brainmask.nii.gz > out.1D' )
  system( instr )
  dataIn <- read.table( 'out.1D', as.is=TRUE )
  pTable <- round( dataIn[,1] / sum( dataIn[,1]), 3 )
  pIdx <- pTable<0.05
  cLimit <- dataIn[which( pIdx )[1],1]
  instr <- sprintf('3dclust -prefix brainmask_corr.nii.gz 0 %1.0f brainmask.nii.gz', cLimit+1 )
  system( instr )
  system('3dcalc -a brainmask_corr.nii.gz -b white_matter_mask_N4.nii.gz -expr \u0027and(a,not(b))\u0027 -prefix ribbon.nii.gz')
  instr <- sprintf('nighres_cruise.sh white_matter_mask_N4.nii.gz ribbon.nii.gz 0.90')
  system( instr )  
  
  setwd('cruise_output')
  system('cp cruise_segmentation_cruise_cortex.nii.gz ../')
  setwd( sprintf( '%s/%s_preproc/', outputDir, inDir[ nDir ] ) )
  system('3dcalc -a cruise_segmentation_cruise_cortex.nii.gz -expr \u0027within(a,2,2)\u0027 -prefix wm_cruise.nii.gz')
  system('3dcalc -a cruise_segmentation_cruise_cortex.nii.gz -expr \u0027within(a,1,1)\u0027 -prefix gm_cruise.nii.gz')
  
  
  #instr <- sprintf('nighres_volumetric.sh wm_cruise.nii.gz gm_cruise.nii.gz')
  #system( instr )  

  #instr <- sprintf('nighres_profiles.sh volumetric_output/volumetric_layering_boundaries.nii.gz anatomy_noCB.nii.gz')
  #system( instr )  
  
  #bmk <- antsImageRead( '_ttt_T1_orient_uni_ss_mask_dil.nii.gz' )
  #img <- antsImageRead( 'anatomy_noCB.nii.gz' )
  #segs <- kmeansSegmentation( img, 3, bmk, mrf=0.15 )
  #antsImageWrite( segs$segmentation, '_ttt_T1MPRAGE_SS_segs_N4.nii.gz' ) 

  #instr <- sprintf('3dcalc -a _ttt_T1MPRAGE_SS_segs_N4.nii.gz -expr \u0027within(a,3,3)\u0027 -prefix whiteMatter_ants.nii.gz')
  #system( instr )  
  #instr <- sprintf('3dcalc -a _ttt_T1MPRAGE_SS_segs_N4.nii.gz -expr \u0027within(a,2,2)\u0027 -prefix grayMatter_ants.nii.gz')
  #system( instr )  

  #instr <- sprintf('3dclust 0 10 whiteMatter_ants.nii.gz > out.1D' )
  #system( instr )
  #dataIn <- read.table( 'out.1D', as.is=TRUE )
  #pTable <- round( dataIn[,1] / sum( dataIn[,1]), 3 )
  #pIdx <- pTable<0.05
  #cLimit <- dataIn[which( pIdx )[1],1]
  #instr <- sprintf('3dclust -prefix whiteMatter_ants_clust.nii.gz 0 %1.0f whiteMatter_ants.nii.gz', cLimit+1 )
  #system( instr )

  #instr <- sprintf('3dclust 0 10 grayMatter_ants.nii.gz > out.1D' )
  #system( instr )
  #dataIn <- read.table( 'out.1D', as.is=TRUE )
  #pTable <- round( dataIn[,1] / sum( dataIn[,1]), 3 )
  #pIdx <- pTable<0.05
  #cLimit <- dataIn[which( pIdx )[1],1]
  #instr <- sprintf('3dclust -prefix grayMatter_ants_clust.nii.gz 0 %1.0f grayMatter_ants.nii.gz', cLimit+1 )
  #system( instr )
  
  #instr <- sprintf('nighres_cruise.sh whiteMatter_ants_clust.nii.gz grayMatter_ants_clust.nii.gz')
  #system( instr )  
  
  #change this param (0.5) and test:
  #library(pracma)
  #plot( sigmoid( seq(-5,5,0.1), 0.5, 0 )~seq(-5,5,0.1) )
  
  #segIn <- antsImageRead( '_ttt_T1MPRAGE_SS_N4.nii.gz' )
  
  
  #volumeAnts <- antsImageRead( '_ttt_T1MPRAGE_SS.nii.gz' ); #here I use ANTsR, modify
  #volumeAntsN4 <- n4BiasFieldCorrection( volumeAnts, 5 )
  #antsImageWrite( volumeAntsN4, '_ttt_T1MPRAGE_SS_N4.nii.gz' ) 

  #bmk <- antsImageRead( '_ttt_T1_orient_uni_ss_mask_dil.nii.gz' )
  #img <- antsImageRead( '_ttt_T1MPRAGE_SS_N4.nii.gz' )
  #img <- n3BiasFieldCorrection( img , 4 )
  #img <- n3BiasFieldCorrection( img , 2 )
  #segs <- kmeansSegmentation( img, 3, bmk )
  #priors <- segs$probabilityimages
  #seg <- geoSeg( img, bmk, priors )
  #antsImageWrite( seg, '_ttt_T1MPRAGE_SS_seg.nii.gz' ) 
  
  #antsImageWrite( segs$segmentation, '_ttt_T1MPRAGE_SS_segs.nii.gz' ) 
  
}