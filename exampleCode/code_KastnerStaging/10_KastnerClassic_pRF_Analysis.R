rm(list=ls())

# toolbox and functions
generalPurposeDir <- Sys.getenv( x='AFNI_TOOLBOXDIRGENERALPURPOSE' )
afniInstallDir <- Sys.getenv( x='AFNI_INSTALLDIR' ) 
source( sprintf('%s/scaleData.R', generalPurposeDir) )
source( sprintf('%s/AFNIio.R', afniInstallDir ) )
library( pracma )
library( abind )
library( circular )

# folders 
mainFolder <- '/analyse/Project0226/KastnerModel/staging_area_Kastner'
prfDir <- sprintf( '%s/expConditions/bars', mainFolder )
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelsDir <- sprintf('%s/modelsOutput', mainFolder )

setwd( prfDir )
print( sprintf('current folder: %s', getwd() ) )
pRFParticipants <- dir( prfDir )
print( sprintf('pRF participants:' ) )
pRFParticipants

setwd( anatomyDir )
print( sprintf('current folder: %s', getwd() ) )
anatomyFolders <- dir( anatomyDir, pattern='*_ANATOMY' )
print( sprintf('anatomy folders:' ) )
anatomyFolders

# look for corresponding anatomy folders
anatomyFolders_id <- rep('a',length(anatomyFolders))
for (nAnat in 1:length(anatomyFolders)) {
  anatomyFolders_id[nAnat] <- strsplit( anatomyFolders, '_' )[[nAnat]][1]
}
anatomyFolders_matching <- rep(0,length(pRFParticipants))
for ( nPart in 1:length( pRFParticipants ) ) {
  prfPartId <- strsplit( pRFParticipants[ nPart ], '_' )[[1]][1]
  anatIdx <- which( anatomyFolders_id == prfPartId )
  print( anatIdx )
  anatomyFolders_matching[nPart] <- anatIdx
}

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
data.frame( pRFParticipants, selectedAnatFolders )

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

#### create outputFolder folder, pRF ####
dirToCheck <- sprintf('%s/pRF', modelsDir )
if ( dir.exists( dirToCheck ) ) {
  instr <- sprintf( 'rm -R %s', dirToCheck ) 
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
}  
if ( !dir.exists( dirToCheck ) ) {
  instr <- sprintf('mkdir %s', dirToCheck)
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
}         

for ( nPart in 1:length(pRFParticipants) ) { # nPart <- 1 1:length(pRFParticipants)
  
  epiPrfTemp <- sprintf('%s/%s', prfDir, pRFParticipants[nPart] )
  anatomyTemp <- sprintf('%s/%s/FreeSeg_result/SUMA', anatomyDir, selectedAnatFolders[nPart] )
  print( epiPrfTemp )
  print( anatomyTemp )
  participantName <- strsplit( selectedAnatFolders[nPart], '_' )[[1]][1]
  print( participantName )
  
  #### create outputFolder folder, kastnerClassic ####
  dirToCheck <- sprintf('%s/pRF/%s', modelsDir, participantName )
  if ( dir.exists( dirToCheck ) ) {
    instr <- sprintf( 'rm -R %s', dirToCheck ) 
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }  
  if ( !dir.exists( dirToCheck ) ) {
    instr <- sprintf('mkdir %s', dirToCheck)
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }         
  
  setwd( dirToCheck )
  print( getwd() )  
  
  # lh ribbon  
  instr <- paste('3dresample',
                 sprintf( '-input %s/lh.ribbon.nii', anatomyTemp ),
                 sprintf( '-master %s', epiPrfTemp ),  
                 sprintf( '-prefix %s/%s/%s/%s.zzz.lh.ribbon.nii', modelsDir, 'pRF', participantName, selectedAnatFolders[nPart] ), 
                 '-rmode NN' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
  
  # rh ribbon  
  instr <- paste('3dresample',
                 sprintf( '-input %s/rh.ribbon.nii', anatomyTemp ),
                 sprintf( '-master %s', epiPrfTemp ),  
                 sprintf( '-prefix %s/%s/%s/%s.zzz.rh.ribbon.nii', modelsDir, 'pRF', participantName, selectedAnatFolders[nPart] ), 
                 '-rmode NN' )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }  
  
  # combine ribbons
  instr <- paste('3dcalc', 
                 sprintf( '-a %s/%s/%s/%s.zzz.lh.ribbon.nii', modelsDir, 'pRF', participantName, selectedAnatFolders[nPart] ), 
                 sprintf( '-b %s/%s/%s/%s.zzz.rh.ribbon.nii', modelsDir, 'pRF', participantName, selectedAnatFolders[nPart] ), 
                 sprintf('-expr \u0027a+b\u0027'),
                 sprintf( '-prefix %s/%s/%s/%s.zzz.grayMatter.nii', modelsDir, 'pRF', participantName, selectedAnatFolders[nPart] ) ) 
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  # copy pRF data
  instr <- paste('cp',
                 sprintf( '%s', epiPrfTemp ),  
                 sprintf( '%s/%s/%s/%s', modelsDir, 'pRF', participantName, pRFParticipants[nPart] )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #blur pRF
  instr <- paste('3dmerge',
                 '-1blur_fwhm 0.1',
                 '-doall',
                 sprintf('-prefix %s_pRF_blur.nii.gz', participantName ),
                 sprintf( '%s', pRFParticipants[nPart]  ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  #mean pRF
  instr <- paste('3dTstat',
                 sprintf( '-mean -prefix %s_pRF_mean.nii', participantName ),
                 sprintf('%s_pRF_blur.nii.gz', participantName ) )
  print( instr )
  if (runCodeFlag==1) { system( instr ) } 
  
  # scale pRF
  instr <- paste('3dcalc',
                 sprintf('-a %s_pRF_blur.nii.gz', participantName ),
                 sprintf( '-b %s_pRF_mean.nii', participantName ),
                 sprintf( '-expr \u0027 min(200, a/b*100)*step(a)*step(b) \u0027' ),
                 sprintf( '-prefix %s_zzz_pRF_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # detrend cw
  instr <- paste('3dDetrend',
                 sprintf('-polort 3'),
                 sprintf( '-prefix %s_zzz_pRF_detrended.nii', participantName ),
                 sprintf( '%s_zzz_pRF_scaled.nii', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #pRF model predictions 
  instr <- paste('Rscript',
                 sprintf('%s/code/models/generatePredictionsPrf_PSOCK_noMultParameter.R', mainFolder ),
                 sprintf('%s_zzz_pRF_detrended.nii', participantName ),
                 sprintf('4 0 0 0 '),
                 sprintf('%s_pRF_predictions/', participantName )
  )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # get mask
  instr <- sprintf('3dmask_tool -input %s_ANATOMY.zzz.grayMatter.nii -prefix %s_mask_auto.nii.gz -dilate_input 1', participantName, participantName )  
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #pRF model fit
  instr <- paste('Rscript',
                 sprintf('%s/code/models/fitPrf_PSOCK_noMultParameter.R', mainFolder), 
                 sprintf('%s_zzz_pRF_detrended.nii', participantName ), 
                 sprintf('%s_mask_auto.nii.gz', participantName ), 
                 sprintf('%s_pRF_predictions/', participantName ),
                 sprintf('%s_pRF', participantName ), 
                 sprintf('0 0') )
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  # clean up
  setwd( dirToCheck )
  print( getwd() )
  system('rm *zzz*')
  system( sprintf('rm %s', pRFParticipants[nPart] ) )
  system( sprintf('rm %s_pRF_blur.nii.gz', participantName ) )

}

