rm(list=ls())
epiDir <- '/analyse/Project0226/KastnerModel/staging_area_Kastner/epi_data'
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

runCodeFlag <- 0

for ( nSession in 1:length( epiSessions ) ) {
  
  if ( nSession==1 ) {
    mainDir <- sprintf('%s/%s', epiSessions[ nSession ] )
    anatomyDir <- sprintf('%s/%s' anatomyDir, anatomyFolders[1] )
  }
    
  setwd( mainDir )
  getwd()
  
  #### mildly mask the data ####
  if ( dir.exists('EPI_mask/') ) {
	instr <- 'rm -R EPI_mask/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }
  }
  instr <- 'automask_whole_dir_afni.sh EPI/ 0.1'; 
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  #### denoise the data ####
  if ( dir.exists('EPI_denoise/') ) {
  	instr <- 'rm -R EPI_denoise/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }  
  }
  instr <- 'denoiseData.sh EPI_mask/ *.nii 8 EPI_denoise/' 
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #### time slice correction ####
  if ( dir.exists('EPI_timeSlice/') ) { 
  	instr <- 'rm -R EPI_timeSlice/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }    
  }
  instr <- 'timeSliceCorrection_noRJSON.sh EPI_denoise/ *.nii EPI_json/ EPI_timeSlice/'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  #### compute top up correction and correction ####
  if ( dir.exists('topUpDir/') ) { 
  	instr <- 'rm -R topUpDir/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }      
  }
  if ( dir.exists('motionCorrectEpi/') ) { 
   	instr <- 'rm -R motionCorrectEpi/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }      
  }
  if ( dir.exists('motionCorrect_topUp_Epi/') ) { 
     	instr <- 'rm -R motionCorrect_topUp_Epi/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }      
  }
  if ( dir.exists('motionCorrectEpiForTopUp/') ) { 
        instr <- 'rm -R motionCorrectEpiForTopUp/'
	print( instr )
        if (runCodeFlag==1) { system( instr ) }      
  }
  instr <- 'motionCorrect.afni.blip.sh EPI_timeSlice/ TOPUP/ 3 1 3-6-8 1-2-3 7 -1'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }

  instr <- 'motionCorrectEPI.with.topUp.sh EPI_timeSlice/ topUpDir/ 3'
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
    
}


  #### get EPI for coregistration and mean TS for analysis ####
  #if ( file.exists('meanEpi4Coreg.nii.gz') ) { system('rm meanEpi4Coreg.nii.gz') }
  #if ( file.exists('meanTs.nii') ) { system('rm meanTs.nii') }
  #instr <- 'computeMeanTs.sh motionCorrect_topUp_Epi/'
  #system( instr )
  #instr <- '3dTstat -prefix meanEpi4Coreg.nii.gz meanTs.nii'; system( instr )
  #system('mv meanEpi4Coreg.nii.gz coregistration/meanEpi4Coreg.nii.gz')
  #system('cp ANATOMY/anatCopy.nii.gz coregistration/anatCopy.nii.gz')
  #system( sprintf( 'cp ANATOMY/anatCopy.nii.gz %s/anatCopy.nii.gz', anatomyDir ) )
  
  #### compute detrended meanTS and move mean TS in a dedicated folder ####
  #if ( dir.exists('motionCorrect_topUp_Epi_detrend/') ) { system('rm -R motionCorrect_topUp_Epi_detrend/') }
  #if ( dir.exists('meanTsFolder/') ) { system('rm -R meanTsFolder/') }
  #system('mkdir meanTsFolder/')
  #system('rm meanTs.nii') # from previous computeMeanTs.sh instruction, without detrending, on the original data
  #system('detrend_whole_dir.sh motionCorrect_topUp_Epi/ *.BRIK 3')
  #instr <- 'computeMeanTs.sh motionCorrect_topUp_Epi_detrend/'; system( instr )
  #system('mv meanTs.nii meanTsFolder/meanTs.nii')


