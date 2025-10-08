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

for ( nSession in 1:length( epiSessions ) ) { # nSession <- 1
  
  mainDir <- sprintf('%s/%s', epiDir, epiSessions[ nSession ] )
  anatomyDirLoop <- sprintf('%s/%s', anatomyDir, anatomyPerSession[ nSession ] )
  
  setwd( mainDir )
  getwd()
  print( getwd() )
  print( sprintf('session: %d', nSession ) )
  
  #### create epiOnAanat folder ####
  if ( dir.exists('epiOnAnat/') ) {
    instr <- 'rm -R epiOnAnat' 
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }  
  if ( !dir.exists('epiOnAnat/') ) {
    instr <- 'mkdir epiOnAnat' 
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }         
  
  epiFiles <- list.files( path='data_tsCorr_denoised_distCorr_motCorr', full.names = FALSE, recursive = FALSE )
  
  for ( nFile in 1:length( epiFiles ) ) {
    fileIn <- epiFiles[ nFile ]
    fileOut <- sprintf('%s_anat.nii.gz',strsplit( fileIn, '.nii.gz' )[[1]][1] )
    instr <- paste( '3dAllineate', 
      sprintf('-prefix epiOnAnat/%s', fileOut ),
      sprintf('-1Dmatrix_apply coregister_01/coregMat.1D'),
      sprintf('-final linear') ,
      sprintf('-input data_tsCorr_denoised_distCorr_motCorr/%s', fileIn), 
      sprintf('-master %s/FreeSeg_result/SUMA/anatCopy.nii.gz', anatomyDirLoop),
      sprintf('-mast_dxyz 1.5') )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }
  }
  
}


