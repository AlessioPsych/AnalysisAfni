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

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='6')
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

for ( nSession in 29:35) { # length( epiSessions ) nSession <- 17
  
  if (nSession==1) { # AFH28_17112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==2) { # AFH28_19112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==3) { # ASM16_07112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==4) { # ASM16_08112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==5) { # ASM16_13112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 2-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==6) { # AXS03_07112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-7 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==7) { # AXS03_08112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==8) { # DWR26_16082022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 2-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==9) { # DWR26_18082022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-4-9 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==10) { # EVC01_11082022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==11) { # EVC01_12082022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-4-8 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==12) { # EWE17_06112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 4 1 4-7 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==13) { # EWE17_13112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==14) { # EWE17_19112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==15) { # FSH19_20112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-7 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==16) { # FSH19_21112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 4-8 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==17) { # GKA26_18102022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==18) { # GKA26_22112022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6-9 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==19) { # HMR28_11102019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==20) { # HMR28_12112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==21) { # HMR28_15102019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==22) { # JFE27_08092022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==23) { # JFE27_13092022
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-4-9 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==24) { # KAR27_05112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-5-7 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==25) { # KAR27_29102020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==26) { # KMA25_14032019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==27) { # KMA25_18102019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==28) { # KMA25_20112020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==29) { # KMA25_28092018
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 1 1 1-2-6 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==30) { # LSO17_29102020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6-9 1-2-3'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==31) { # LSO17_31102020
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-7 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==32) { # RDN03_27112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-7 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==33) { # RDN03_28112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==34) { # SCA01_20112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==35) { # SCA01_21112019
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-8 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  
  mainDir <- sprintf('%s/%s', epiDir, epiSessions[ nSession ] )

  setwd( mainDir )
  getwd()
  print( getwd() )
  print( sprintf('session: %d', nSession ) )

  # remove file amplitude anatomy.nii, if it exists
  if ( file.exists('amplitudeAnatomy.nii') ) {
    instr <- 'rm amplitudeAnatomy.nii'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }
  
  #### denoise the data ####
  if ( dir.exists('EPI_denoise/') ) {
    instr <- 'rm -R EPI_denoise/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }
  instr <- 'denoiseData.sh EPI/ *.nii 8 EPI_denoise/' 
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
  
  if ( dir.exists('topUpDir/') ) { 
    instr <- 'rm -R topUpDir/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  instr <- topUpStartInstruction
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  topUpFiles <- dir('topUpDir', pattern='*.nii')
  forwardBlipDataset <- topUpFiles[ idxForwardTopUp ]
  reverseBlipDataset <- topUpFiles[ idxReverseTopUp ]
  if ( dir.exists('subjProcessed.results') ) { 
    instr <- 'rm -R subjProcessed.results'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  if ( file.exists('output.proc.subjProcessed') ) { 
    instr <- 'rm -R output.proc.subjProcessed'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  if ( file.exists('proc.subjProcessed') ) { 
    instr <- 'rm -R proc.subjProcessed'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  
  instr <- paste( 'afni_proc.py',
                  sprintf('-subj_id %s', 'subjProcessed' ),
                  sprintf('-dsets %s/*.nii', 'EPI_timeSlice' ),
                  sprintf('-blocks despike volreg'),
                  sprintf('-blip_forward_dset %s/%s', 'topUpDir', forwardBlipDataset ),
                  sprintf('-blip_reverse_dset %s/%s', 'topUpDir', reverseBlipDataset ),
                  sprintf('-volreg_align_to MIN_OUTLIER'),
                  sprintf('-scr_overwrite'),
                  sprintf('-no_epi_review'),
                  sprintf('-execute') 
  )    
  print( instr )
  if (runCodeFlag==1) { system( instr ) }
  
  if (runCodeFlag==1) {    
    #convert clean files to nifti 
    if ( dir.exists('data_tsCorr_denoised_distCorr_motCorr') ) { 
      instr <- 'rm -R data_tsCorr_denoised_distCorr_motCorr/'
      print( instr )
      if (runCodeFlag==1) { system( instr ) } 
    }	   
    instr <- 'mkdir data_tsCorr_denoised_distCorr_motCorr/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
    
    cleanFilesToConvert <- list.files( path='subjProcessed.results/', pattern='*volreg*' )
    cleanFilesToConvert <- cleanFilesToConvert[seq(1,length(cleanFilesToConvert),2)] # selects only BRIK files
    for (nFile in 1:length(cleanFilesToConvert) ) {
      fileNameTemp01 <- strsplit( cleanFilesToConvert[nFile], '[.]' )
      fileNameTemp <- paste( fileNameTemp01[[1]][1], '_',  fileNameTemp01[[1]][2], '_', fileNameTemp01[[1]][3], '.nii.gz', sep='')
      instr <- sprintf('3dcopy %s/%s %s/%s', 'subjProcessed.results', cleanFilesToConvert[nFile], 'data_tsCorr_denoised_distCorr_motCorr', fileNameTemp )
      print( instr )
      if (runCodeFlag==1) { system( instr ) }    
    }
    
    #copy motion correction files 
    if ( dir.exists('data_tsCorr_denoised_distCorr_motCorr_motion_param_files') ) { 
      instr <- 'rm -R data_tsCorr_denoised_distCorr_motCorr_motion_param_files/'
      print( instr )
      if (runCodeFlag==1) { system( instr ) }    
    }
    instr <- 'mkdir data_tsCorr_denoised_distCorr_motCorr_motion_param_files/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
    filesToCopy <- list.files( path='subjProcessed.results/', pattern='dfile.*' )
    for (nFile in 1:length(filesToCopy) ) {
      instr <- sprintf('cp %s/%s %s/%s', 'subjProcessed.results', filesToCopy[nFile], 'data_tsCorr_denoised_distCorr_motCorr_motion_param_files', filesToCopy[nFile] )
      print( instr )
      if (runCodeFlag==1) { system( instr ) }    
    }
  }
  
  # clean up
  if ( dir.exists('EPI_timeSlice') ) { 
    instr <- 'rm -R EPI_timeSlice/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  if ( dir.exists('EPI_denoise') ) { 
    instr <- 'rm -R EPI_denoise/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  if ( dir.exists('topUpDir') ) { 
    instr <- 'rm -R topUpDir/'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  if ( dir.exists('subjProcessed.results') ) { 
    instr <- 'rm -R subjProcessed.results'
    print( instr )
    if (runCodeFlag==1) { system( instr ) }    
  }
  
}


