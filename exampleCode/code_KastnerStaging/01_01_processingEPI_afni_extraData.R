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

Sys.setenv(OMP_NUM_THREADS='6')
Sys.getenv('OMP_NUM_THREADS')

### align anatomyFolders to epiSession 
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

for ( nSession in 1:length( epiSessions ) ) { # length( epiSessions ) nSession <- 17
  
  if (nSession==1) { # ASM16_22102020_eyebar2_eyebar3
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==2) { # AXS03_24112020_eyebar2_eyebar3
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==3) { # CDA05_02122020_eyebar2
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-4 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==4) { # CDA05_26112020_bar_eyebar3
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==5) { # GMA01_03112020_eyebar2
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 2 1 2-5 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==6) { # GMA01_30102020_bar_eyebar3
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==7) { # HMR28_15102020_eyebar2_eyebar3
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3-6 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==8) { # SON25_19112020_bar_eyebar3
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 4 1 4-7 1-2'
    idxForwardTopUp <- 1
    idxReverseTopUp <- 2
  }
  if (nSession==9) { # SON25_20112020_eyebar2
    topUpStartInstruction <- 'motionCorrect.afni.blip.noWarp.sh EPI_timeSlice/ TOPUP/ 3 1 3 1'
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


