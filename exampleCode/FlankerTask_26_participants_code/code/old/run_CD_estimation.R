rm( list=ls() )
mainDir <- '/media/alessiofracasso/DATADRIVE1/Flanker'
setwd( mainDir )
outputDir <- 'CD_estimation'

partDirs <- list.dirs(full.names = FALSE, recursive = FALSE)
partDirs <- partDirs[grepl("^*sub-", partDirs)]
print( 'participant dirs...' )
print(partDirs)

print( 'clean up and create main output folder in derivatives...' )
dirToCheck <- sprintf('%s/derivatives/%s', mainDir, outputDir)
flagDir <- dir.exists( dirToCheck  )
if ( flagDir==TRUE ) { 
  system( sprintf('rm -R %s', dirToCheck ) )
  dir.create( dirToCheck )
}
if ( flagDir==FALSE ) { dir.create( dirToCheck ) }

#nCores <- 6
runFlag <- 1

##### run freesurfer and suma ####
for ( i in  1:1 ) { # i <- 1 length( partDirs )
  setwd( mainDir )
  setwd( 'derivatives' )
  setwd( 'Freesurfer_output' )
  setwd( partDirs[ i ] )
  setwd( 'Freesurfer_result/SUMA' )
  print( getwd() )

  #system( sprintf('cp %s/CD_estimation.sh %s/CD_estimation.sh', mainDir, getwd() ) )
  # combine left and right upsampled ribbons
  system( sprintf('3dcalc -a rh.ribbon.nii -b lh.ribbon.nii -expr \u0027a+b\u0027 -prefix gray_matter_ribbon.nii.gz') )
  
  # get white matter mask from aparc+aseg segmentation
  system( sprintf('3dcalc -a aparc+aseg.nii -expr \u0027within(a,1.5,2.5)+within(a,40.5,41.5)\u0027 -prefix white_matter_mask_Freesurfer_res.nii.gz') )
  
  # upsample white matter mask from Freesurfer to anatCopy.nii.gz resolution
  system( sprintf('3dresample -master Freesurfer_result_SurfVol.nii -rmode NN -prefix white_matter_mask.nii.gz -inset white_matter_mask_Freesurfer_res.nii.gz') )
  
  # combine white matter mask and gray matter ribbon at anatCopy.nii.gz resolution
  system( sprintf('3dcalc -a white_matter_mask.nii.gz -b gray_matter_ribbon.nii.gz -expr \u0027a+b\u0027 -prefix gray_matter_mask_out.nii.gz') )
  
  currentDir <- getwd()
  instr <- sprintf('docker run -v %s:/testNighres --rm nighres', currentDir )
  print( instr )
  if ( runFlag==1 ) { system( instr ) } 
  print( '...' )
  
}
