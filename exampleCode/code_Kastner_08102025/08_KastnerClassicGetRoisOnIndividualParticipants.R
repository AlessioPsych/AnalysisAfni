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
anatomyDir <- sprintf( '%s/anatomies_KastnerClassic_Freesurfer', mainFolder )
modelOutputFolder <- sprintf('%s/modelsOutput/kastnerClassic', mainFolder )
surfaceFolder <- sprintf('%s/anatomies_KastnerClassic_Freesurfer/surfaceAtlases/suma_MNI152_2009', mainFolder)
roiSurfaceFolderMNI <- sprintf('%s/probatlas_v4/ProbAtlas_v4/subj_surf_all', surfaceFolder)

setwd( modelOutputFolder )
print( sprintf('current folder: %s', getwd() ) )
modelParticipants <- dir( modelOutputFolder )
print( sprintf('model participants:' ) )
modelParticipants

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
anatomyFolders_matching <- rep(0,length(modelParticipants))
for ( nPart in 1:length( modelParticipants ) ) {
  partId <- strsplit( modelParticipants[ nPart ], '_' )[[1]][1]
  anatIdx <- which( anatomyFolders_id == partId )
  print( anatIdx )
  anatomyFolders_matching[nPart] <- anatIdx
}

# check the alignment
selectedAnatFolders <- anatomyFolders[ anatomyFolders_matching ]
data.frame( modelParticipants, selectedAnatFolders )

runCodeFlag <- 1

Sys.setenv(OMP_NUM_THREADS='8')
Sys.getenv('OMP_NUM_THREADS')

for ( nPart in 1:length( modelParticipants )  ) { # nPart <- 1 # length( modelParticipants )
  
  # move into individual subject SUMA folder
  setwd(anatomyDir)
  setwd( selectedAnatFolders[ nPart ] )
  setwd('FreeSeg_result/SUMA')
  print( getwd() )
  
  #### create roi folder ####
  dirToCheck <- sprintf('roisWang2015_KastnerClassic' )
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
  
  for ( nHemi in c('lh','rh') ) { # nHemi <- 'lh'
    
    # get roi surface filename
    roiVolumeFilename <- sprintf( '%s/maxprob_surf_%s.1D.dset', roiSurfaceFolderMNI, nHemi ) 
    allVolumeFiles <- list.files( path='kastnerClassicModelsAnatomy', recursive = FALSE, full.names = FALSE )
    volumeIndividualParticipant <- allVolumeFiles[1]
    
    instr <- paste('3dSurf2Vol',
                   sprintf( '-spec std.141.FreeSeg_result_both.spec'),
                   sprintf( '-surf_A std.141.%s.smoothwm.gii', nHemi ),
                   sprintf( '-surf_B std.141.%s.pial.gii', nHemi ),
                   sprintf( '-sdata_1D %s', roiVolumeFilename ),
                   sprintf( '-sv %s/%s', 'kastnerClassicModelsAnatomy', volumeIndividualParticipant ),
                   sprintf( '-grid_parent %s/%s', 'kastnerClassicModelsAnatomy', volumeIndividualParticipant ),
                   sprintf( '-map_func nzmode' ),
                   sprintf( '-f_steps 15' ),
                   sprintf( '-prefix %s/%s.%s.nii.gz', 'roisWang2015_KastnerClassic', 'WangRoi2015', nHemi )
    )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
    
    # copy original anatomy in roi folder (to visual check only)
    instr <- paste('cp',
                   'anatCopy.nii.gz',
                   'roisWang2015_KastnerClassic/anatCopy.nii.gz')
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
    
  }
  
}
    
