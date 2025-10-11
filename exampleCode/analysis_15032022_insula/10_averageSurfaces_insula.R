rm( list=ls() ); gc();
graphics.off();

freesurferDirPart01 <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/part1'
freesurferDirPart02 <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/part2'
mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insula'
atlasDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/suma_MNI152_2009_princetonAtlas'
saveDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insulaSurfaceFromVolumetricAnalysis_poly01'
aHEADDataBaseDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840'
filenameIn <- 'poly01'
setwd( mainDir )

set.seed(1);
library( plotly )
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniAtlasDir <- Sys.getenv(x='AFNI_ATLASDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )
source( sprintf( '%s/linearIndexFromCoordinate.r', afniSurfaces ) )


for (freesurferDirCounter in 1:2) {
  if (freesurferDirCounter==1) {
    freesurferDir <- freesurferDirPart01
    setwd( freesurferDir )
    participantsFolders <- dir( pattern = 'sub-*')
    participantsFolders <- participantsFolders[-c(8,50,51)] #to exclude subject 0007 and 0049 and 0050 due to processing problems
  }else{
    freesurferDir <- freesurferDirPart02
    setwd( freesurferDir )
    participantsFolders <- dir( pattern = 'sub-*')
    participantsFolders <- participantsFolders[-c(25)] #to exclude subject 0075
  }
  
  for ( i in 1:length( participantsFolders ) ) { #i <- 1
    print( sprintf('copy surface data, subject: %s', participantsFolders[i] ) )
    setwd( freesurferDir ); getwd()
    setwd( participantsFolders[i] ); getwd()
    setwd( 'ses-1' ); getwd()
    setwd( 'anat' ); getwd()
    setwd( 'AHEAD_test' ); getwd()
    setwd( 'SUMA' ); getwd()
    
    loadFileNameLh <- sprintf('insulaClusters_%s_lh.1D.dset', filenameIn ) 
    loadFileNameRh <- sprintf('insulaClusters_%s_rh.1D.dset', filenameIn ) 
    
    instr <- sprintf('cp %s %s/%s_insulaClusters_lh.1D.dset', loadFileNameLh, saveDir, participantsFolders[i] ); system( instr )
    instr <- sprintf('cp %s %s/%s_insulaClusters_rh.1D.dset', loadFileNameRh, saveDir, participantsFolders[i] ); system( instr )
    
    instr <- 'ConvertDset -o_1D -input std.141.rh.curv.niml.dset -prepend_node_index_1D -prefix std.141.rh.curv_ppp_.1D.dset'; system( instr )
    instr <- 'ConvertDset -o_1D -input std.141.lh.curv.niml.dset -prepend_node_index_1D -prefix std.141.lh.curv_ppp_.1D.dset'; system( instr )
    instr <- sprintf('cp std.141.rh.curv_ppp_.1D.dset %s/%s_std.141.rh.curv_ppp_.1D.dset', saveDir, participantsFolders[i] ); system( instr )
    instr <- sprintf('cp std.141.lh.curv_ppp_.1D.dset %s/%s_std.141.lh.curv_ppp_.1D.dset', saveDir, participantsFolders[i] ); system( instr )
    
    instr <- 'ConvertDset -o_1D -input std.141.rh.thickness.niml.dset -prepend_node_index_1D -prefix std.141.rh.thickness_ppp_.1D.dset'; system( instr )
    instr <- 'ConvertDset -o_1D -input std.141.lh.thickness.niml.dset -prepend_node_index_1D -prefix std.141.lh.thickness_ppp_.1D.dset'; system( instr )
    instr <- sprintf('cp std.141.rh.thickness_ppp_.1D.dset %s/%s_std.141.rh.thickness_ppp_.1D.dset', saveDir, participantsFolders[i] ); system( instr )
    instr <- sprintf('cp std.141.lh.thickness_ppp_.1D.dset %s/%s_std.141.lh.thickness_ppp_.1D.dset', saveDir, participantsFolders[i] ); system( instr )
    
    system('rm *_ppp_*')
  }
}

# left hemisphere
setwd( saveDir )
surfaceFile <- dir( pattern = '*_lh.1D.dset*')
surfaceFile_thickness <- dir( pattern = '*lh.thickness_ppp_*')
surfaceFile_curvature <- dir( pattern = '*lh.curv_ppp_*')
storeSingleSubjectMap <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_01 <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_02 <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_thickness <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_curvature <- array( 9999, c( 198811, length( surfaceFile ) ) )
for ( i in 1:length(surfaceFile) ) { # i <- 1
  print( sprintf('load subject file: %s', surfaceFile[i] ) )
  surfaceData <- read.table( surfaceFile[i], as.is=TRUE, header=TRUE )
  storeSingleSubjectMap[,i] <- surfaceData[,7]
  dataApp <- surfaceData[,7]
  dataApp <- ifelse( dataApp==1, 1, 0 )
  storeSingleSubjectMap_01[,i] <- dataApp
  dataApp <- surfaceData[,7]
  dataApp <- ifelse( dataApp==2, 1, 0 )
  storeSingleSubjectMap_02[,i] <- dataApp
  
  surfaceData_thickness <- read.table( surfaceFile_thickness[i], as.is=TRUE, header=TRUE )
  storeSingleSubjectMap_thickness[,i] <- surfaceData_thickness[,2]
  
  surfaceData_curvature <- read.table( surfaceFile_curvature[i], as.is=TRUE, header=TRUE )
  storeSingleSubjectMap_curvature[,i] <- surfaceData_curvature[,2]
}

averageLH <- round( apply( storeSingleSubjectMap, 1, mean ), 4 ) # average cluster [1,2]
averageLH_01 <- round( apply( storeSingleSubjectMap_01, 1, mean ), 4 ) # proportion cluster 1
averageLH_02 <- round( apply( storeSingleSubjectMap_02, 1, mean ), 4 ) # proportion cluster 2
averageLH_thickness <- round( apply( storeSingleSubjectMap_thickness, 1, mean ), 4 ) # average thickness
averageLH_curvature <- round( apply( storeSingleSubjectMap_curvature, 1, mean ), 4 ) # average curvature

surfaceData[,7] <- averageLH
surfaceData[,8] <- averageLH_01
surfaceData[,9] <- averageLH_02
surfaceData[,10] <- averageLH_thickness 
surfaceData[,11] <- averageLH_curvature
write.table( surfaceData, file=sprintf( '%s/average_insulaClusters_%s_lh_corrected.1D.dset', atlasDir, filenameIn ), row.names=FALSE, col.names=FALSE )

# right hemisphere
setwd( saveDir )
surfaceFile <- dir( pattern = '*_rh.1D.dset*')
surfaceFile_thickness <- dir( pattern = '*rh.thickness_ppp_*')
surfaceFile_curvature <- dir( pattern = '*rh.curv_ppp_*')
storeSingleSubjectMap <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_01 <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_02 <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_thickness <- array( 9999, c( 198811, length( surfaceFile ) ) )
storeSingleSubjectMap_curvature <- array( 9999, c( 198811, length( surfaceFile ) ) )
for ( i in 1:length(surfaceFile) ) { # i <- 1
  print( sprintf('load subject file: %s', surfaceFile[i] ) )
  surfaceData <- read.table( surfaceFile[i], as.is=TRUE, header=TRUE )
  storeSingleSubjectMap[,i] <- surfaceData[,7]
  dataApp <- surfaceData[,7]
  dataApp <- ifelse( dataApp==1, 1, 0 )
  storeSingleSubjectMap_01[,i] <- dataApp
  dataApp <- surfaceData[,7]
  dataApp <- ifelse( dataApp==2, 1, 0 )
  storeSingleSubjectMap_02[,i] <- dataApp
  
  surfaceData_thickness <- read.table( surfaceFile_thickness[i], as.is=TRUE, header=TRUE )
  storeSingleSubjectMap_thickness[,i] <- surfaceData_thickness[,2]
  
  surfaceData_curvature <- read.table( surfaceFile_curvature[i], as.is=TRUE, header=TRUE )
  storeSingleSubjectMap_curvature[,i] <- surfaceData_curvature[,2]
}

averageRH <- round( apply( storeSingleSubjectMap, 1, mean ), 4 ) # average cluster [1,2]
averageRH_01 <- round( apply( storeSingleSubjectMap_01, 1, mean ), 4 ) # proportion cluster 1
averageRH_02 <- round( apply( storeSingleSubjectMap_02, 1, mean ), 4 ) # proportion cluster 2
averageRH_thickness <- round( apply( storeSingleSubjectMap_thickness, 1, mean ), 4 ) # average thickness
averageRH_curvature <- round( apply( storeSingleSubjectMap_curvature, 1, mean ), 4 ) # average curvature

surfaceData[,7] <- averageRH
surfaceData[,8] <- averageRH_01
surfaceData[,9] <- averageRH_02
surfaceData[,10] <- averageRH_thickness 
surfaceData[,11] <- averageRH_curvature
write.table( surfaceData, file=sprintf( '%s/average_insulaClusters_%s_rh_corrected.1D.dset', atlasDir, filenameIn ), row.names=FALSE, col.names=FALSE )


#####  limit by age ####


aHEADDataBase <- read.csv( sprintf('%s/participants.csv',aHEADDataBaseDir), as.is = TRUE )
table( aHEADDataBase$Group )
#selectedParticipants <- c( aHEADDataBase$ScanName[ aHEADDataBase$Group=='18-30' ], aHEADDataBase$ScanName[ aHEADDataBase$Group=='31-40' ] )
#selectedParticipants <- selectedParticipants[-c(23)] #to exclude subject 0007 and 0049 and 0050 due to processing problems

Age_group1 = c( aHEADDataBase$ScanName[ aHEADDataBase$Group=='18-30' ], aHEADDataBase$ScanName[ aHEADDataBase$Group=='31-40' ] )
Age_group1 <- Age_group1[-c(23,34)] #to exclude subject 0049 and 0075 due to processing problems
Age_group2 = aHEADDataBase$ScanName[!(aHEADDataBase$Group %in% c('18-30', '31-40'))]
Age_group2 =  Age_group2[-c(5,26)] #to exclude subject 0007 and 0050 due to processing problems
AgeGroups = list(Age_group1,Age_group2)

for (group in AgeGroups){
  group_name <- ifelse(group %in% Age_group1, "18_40", "41Plus")
  group_name = group_name[1]
  print(group_name)
  selectedParticipants = group
  
  
  # left hemisphere
  setwd( saveDir )
  surfaceFile <- dir( pattern = '*_lh.1D.dset*')
  surfaceFile_thickness <- dir( pattern = '*lh.thickness_ppp_*')
  surfaceFile_curvature <- dir( pattern = '*lh.curv_ppp_*')
  storeSingleSubjectMap <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_01 <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_02 <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_thickness <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_curvature <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  counter <- 1
  for ( i in 1:length(surfaceFile) ) { # i <- 1
    
    particiantLoop <- strsplit( surfaceFile[i], '_' )[[1]][1]
    if (particiantLoop %in% selectedParticipants) {
      print( sprintf('load subject file: %s', surfaceFile[i] ) )
      surfaceData <- read.table( surfaceFile[i], as.is=TRUE, header=TRUE )
      storeSingleSubjectMap[,counter] <- surfaceData[,7]
      dataApp <- surfaceData[,7]
      dataApp <- ifelse( dataApp==1, 1, 0 )
      storeSingleSubjectMap_01[,counter] <- dataApp
      dataApp <- surfaceData[,7]
      dataApp <- ifelse( dataApp==2, 1, 0 )
      storeSingleSubjectMap_02[,counter] <- dataApp
      
      surfaceData_thickness <- read.table( surfaceFile_thickness[i], as.is=TRUE, header=TRUE )
      storeSingleSubjectMap_thickness[,counter] <- surfaceData_thickness[,2]
      
      surfaceData_curvature <- read.table( surfaceFile_curvature[i], as.is=TRUE, header=TRUE )
      storeSingleSubjectMap_curvature[,counter] <- surfaceData_curvature[,2]
      counter <- counter + 1
    }
    
    
    
  }
  
  averageLH <- round( apply( storeSingleSubjectMap, 1, mean ), 4 ) # average cluster [1,2]
  averageLH_01 <- round( apply( storeSingleSubjectMap_01, 1, mean ), 4 ) # proportion cluster 1
  averageLH_02 <- round( apply( storeSingleSubjectMap_02, 1, mean ), 4 ) # proportion cluster 2
  averageLH_thickness <- round( apply( storeSingleSubjectMap_thickness, 1, mean ), 4 ) # average thickness
  averageLH_curvature <- round( apply( storeSingleSubjectMap_curvature, 1, mean ), 4 ) # average curvature
  
  surfaceData[,7] <- averageLH
  surfaceData[,8] <- averageLH_01
  surfaceData[,9] <- averageLH_02
  surfaceData[,10] <- averageLH_thickness 
  surfaceData[,11] <- averageLH_curvature
  

  #write.table( surfaceData, file=sprintf( '%s/average_insulaClusters_%s_lh_selectedParticipants_corrected.1D.dset', atlasDir, filenameIn ), row.names=FALSE, col.names=FALSE )
  write.table( surfaceData, file=sprintf( '%s/average_insulaClusters_%s_lh_%s.1D.dset', atlasDir, filenameIn , group_name ), row.names=FALSE, col.names=FALSE )
  
  # right hemisphere
  setwd( saveDir )
  surfaceFile <- dir( pattern = '*_rh.1D.dset*')
  surfaceFile_thickness <- dir( pattern = '*rh.thickness_ppp_*')
  surfaceFile_curvature <- dir( pattern = '*rh.curv_ppp_*')
  storeSingleSubjectMap <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_01 <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_02 <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_thickness <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  storeSingleSubjectMap_curvature <- array( 9999, c( 198811, length( selectedParticipants ) ) )
  counter <- 1
  for ( i in 1:length(surfaceFile) ) { # i <- 1
    
    particiantLoop <- strsplit( surfaceFile[i], '_' )[[1]][1]
    if ( is.element( particiantLoop, selectedParticipants ) ) {
      print( sprintf('load subject file: %s', surfaceFile[i] ) )
      surfaceData <- read.table( surfaceFile[i], as.is=TRUE, header=TRUE )
      storeSingleSubjectMap[,counter] <- surfaceData[,7]
      dataApp <- surfaceData[,7]
      dataApp <- ifelse( dataApp==1, 1, 0 )
      storeSingleSubjectMap_01[,counter] <- dataApp
      dataApp <- surfaceData[,7]
      dataApp <- ifelse( dataApp==2, 1, 0 )
      storeSingleSubjectMap_02[,counter] <- dataApp
      
      surfaceData_thickness <- read.table( surfaceFile_thickness[i], as.is=TRUE, header=TRUE )
      storeSingleSubjectMap_thickness[,counter] <- surfaceData_thickness[,2]
      
      surfaceData_curvature <- read.table( surfaceFile_curvature[i], as.is=TRUE, header=TRUE )
      storeSingleSubjectMap_curvature[,counter] <- surfaceData_curvature[,2]
      counter <- counter + 1
    }
  }
  
  averageRH <- round( apply( storeSingleSubjectMap, 1, mean ), 4 ) # average cluster [1,2]
  averageRH_01 <- round( apply( storeSingleSubjectMap_01, 1, mean ), 4 ) # proportion cluster 1
  averageRH_02 <- round( apply( storeSingleSubjectMap_02, 1, mean ), 4 ) # proportion cluster 2
  averageRH_thickness <- round( apply( storeSingleSubjectMap_thickness, 1, mean ), 4 ) # average thickness
  averageRH_curvature <- round( apply( storeSingleSubjectMap_curvature, 1, mean ), 4 ) # average curvature
  
  surfaceData[,7] <- averageRH
  surfaceData[,8] <- averageRH_01
  surfaceData[,9] <- averageRH_02
  surfaceData[,10] <- averageRH_thickness 
  surfaceData[,11] <- averageRH_curvature
  

  #write.table( surfaceData, file=sprintf( '%s/average_insulaClusters_%s_rh_selectedParticipants_corrected.1D.dset', atlasDir, filenameIn ), row.names=FALSE, col.names=FALSE )
  write.table( surfaceData, file=sprintf( '%s/average_insulaClusters_%s_rh_%s.1D.dset', atlasDir, filenameIn , group_name ), row.names=FALSE, col.names=FALSE )

}


rm( list=ls() ); gc();
load("/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/analysis_15032022/09_data_vonEconomo1_lm_save.RData")
vnp1 = high_low
ppn1 = kUmap
load("/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/analysis_15032022/09_data_vonEconomo2_lm_save.RData")
high_low = rbind(vnp1,high_low)
ku = rbind(ppn1,kUmap)
# load("/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/analysis_15032022/09_data_lm_save.RData")
###### linear model #####
library(lme4)
high_low$clusterLocalization <- relevel(high_low$clusterLocalization, ref='high')
Age_group1_lm = high_low %>% filter(high_low$subjArray %in% Age_group1)
Age_group2_lm = high_low %>% filter(high_low$subjArray %in% Age_group2)
options(scipen=999)
summary( lmer( t1Intensity ~ corticalDepth * clusterLocalization + ( 1 | subjArray ), data=Age_group1_lm ) )
summary( lmer( t1Intensity ~ corticalDepth * clusterLocalization + ( 1 | subjArray ), data=Age_group2_lm ) )


