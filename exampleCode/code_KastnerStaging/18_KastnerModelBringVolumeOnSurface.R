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
modelOutputFolder <- sprintf('%s/modelsOutput/kastnerBaseline', mainFolder )

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
  
  # get models names:
  setwd( modelOutputFolder )
  setwd( modelParticipants[ nPart ] )
  print('#####')
  print( getwd() )
  print('#####')
  KastnerModel <- list.files( getwd(), patter='params', full.names = FALSE, recursive = FALSE)
  allModels <- c( KastnerModel )
  
  #### create outputFolder folder, KastnerNewModel ####
  setwd( anatomyDir )
  setwd( selectedAnatFolders[ nPart ] )
  setwd('FreeSeg_result/SUMA')
  print('#####')
  print( getwd() )
  print('#####')
  dirToCheck <- sprintf('KastnerNewModelAnatomy' )
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
  
  #### create outputFolder folder, KastnerNewModel ####
  setwd( anatomyDir )
  setwd( selectedAnatFolders[ nPart ] )
  setwd('FreeSeg_result/SUMA')
  print('#####')
  print( getwd() )
  print('#####')
  dirToCheck <- sprintf('KastnerNewModelSurfacesNiml' )
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
  
  #### create outputFolder folder, pRF ####
  setwd( anatomyDir )
  setwd( selectedAnatFolders[ nPart ] )
  setwd('FreeSeg_result/SUMA')
  print('#####')
  print( getwd() )
  print('#####')
  dirToCheck <- sprintf('KastnerNewModelSurfaces1D' )
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
  
  # copy individual participant, volumetric model files
  for ( nModels in 1:length( allModels ) ) {
    instr <- paste('cp',
                   sprintf('%s/%s/%s', modelOutputFolder, modelParticipants[ nPart ], allModels[ nModels ] ),
                   sprintf('%s/%s/FreeSeg_result/SUMA/KastnerNewModelAnatomy/%s', anatomyDir, selectedAnatFolders[ nPart ], allModels[ nModels ] )
    )
    print( instr )
    if (runCodeFlag==1) { system( instr ) }  
  }
  
  for ( nModels in 1:length( allModels ) ) {
    for (outName in c('niml','1D') ) { # outName <- 'niml'
      for (nHemi in c('lh','rh') ) { # nHemi <- 'lh'
        modelOutputName <- strsplit( allModels[ nModels ], '.nii.gz' )[[1]][1]
        if (outName=='niml') {
          modelOutputNameInstruction <- sprintf( '-out_niml %s/%s.%s.%s.dset', 'KastnerNewModelSurfacesNiml', modelOutputName, nHemi, outName ) 
        }
        if (outName=='1D') {
          modelOutputNameInstruction <- sprintf( '-out_1D %s/%s.%s.%s.dset', 'KastnerNewModelSurfaces1D', modelOutputName, nHemi, outName ) 
        }  
        instr <- paste('3dVol2Surf',
                       sprintf( '-spec std.141.FreeSeg_result_both.spec'),
                       sprintf( '-surf_A std.141.%s.smoothwm.gii', nHemi ),
                       sprintf( '-surf_B std.141.%s.pial.gii', nHemi ),
                       sprintf( '-sv %s/%s', 'KastnerNewModelAnatomy', allModels[ nModels ] ),
                       sprintf( '-grid_parent %s/%s', 'KastnerNewModelAnatomy', allModels[ nModels ] ),
                       sprintf( '-map_func nzave' ),
                       sprintf( '-f_steps 10' ),
                       sprintf( '-f_index nodes' ),                                  
                       modelOutputNameInstruction
        )
        print( instr )
        if (runCodeFlag==1) { system( instr ) }  
      }
    }
  }
}
