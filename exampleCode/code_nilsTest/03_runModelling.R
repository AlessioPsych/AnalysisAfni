rm(list=ls()); gc;
graphics.off()

for (nSubj in 1:1) {
  
  if ( nSubj==1 ) {
    mainDir <- '/analyse/Project0226/tests/nilsTest/data/20230622_SHI27'
    anatomyDir <- '/analyse/Project0226/tests/nilsTest/data/20230622_SHI27/ANATOMY/anat_test/SUMA'
    codeFolder <- '/analyse/Project0226/tests/nilsTest/code/modelling'
  }

  setwd( mainDir )
  print( getwd() )    

  # clean up
  system('rm -R ppp*')
  
  # copy EPI 4 coreg as a base:
  instr <- 'cp coregistration/meanEpi4Coreg.nii.gz ppp_meanEpi4Coreg.nii.gz'
  system( instr )
   
  # copy meanTs coreg as a base:
  instr <- 'cp meanTsFolder/meanTs.nii ppp_meanTs.nii'
  system( instr )
  
  # get mask
  instr <- '3dAutomask -prefix ppp_epi_mask.nii.gz ppp_meanEpi4Coreg.nii.gz'; system( instr )
  
  # smooth data in time
  instr <- '3dTsmooth -hamming 3 -prefix ppp_epi_smooth.nii.gz ppp_meanTs.nii'; system( instr )
  
  # generate predictions _ block design
  instr <- sprintf('Rscript %s/generatePredictions_blockDesign_PSOCK.R stim_data/feat_out.txt ppp_EPI_predictions_blockDesign 300 2', codeFolder )
  print( instr )
  system( instr )
  
  # fit model predictions _ block design 
  instr <- sprintf('Rscript %s/fitModel_PSOCK.R ppp_epi_smooth.nii.gz ppp_epi_mask.nii.gz ppp_EPI_predictions_blockDesign/ ppp_modelOutput_blockDesign 5', codeFolder )
  print( instr )
  system( instr )

  # generate predictions _ block design + amplitude modulation over time
  instr <- sprintf('Rscript %s/blockDesign_scaling_generatePredictions_PSOCK.R stim_data/feat_out.txt ppp_EPI_predictions_ampScaling 300 2', codeFolder )
  print( instr )
  system( instr )
  
  # fit model predictions _ block design + amplitude modulation over time
  instr <- sprintf('Rscript %s/blockDesign_scaling_fitModel_PSOCK.R ppp_epi_smooth.nii.gz ppp_epi_mask.nii.gz ppp_EPI_predictions_ampScaling/ ppp_modelOutput_ampScaling 5 2', codeFolder )
  print( instr )
  system( instr )
  
}

