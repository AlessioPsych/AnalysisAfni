rm( list=ls() ); gc();
graphics.off();

freesurferDirPart01 <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/part1'
freesurferDirPart02 <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/part2'
mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insula'
saveName <- 'poly01'
setwd( mainDir )

set.seed(1);
library( plotly )
library(factoextra)
library(dplyr)
library(tidyr)
library(ggplot2)

#library( ica )
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniAtlasDir <- Sys.getenv(x='AFNI_ATLASDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )
source( sprintf( '%s/linearIndexFromCoordinate.r', afniSurfaces ) )

setwd( mainDir )

filesList <- dir( pattern = '*sub-*' )
filesList <- filesList[-c(15,16,99,100)] #to exclude subject 0007 and 0049 due to processing problems, participant 0050 was already excluded due to processing problems
hemispheres <- c('LH','RH')
storeProfiles01 <- array(0, c( length(filesList), 8 ) )
storeProfiles02 <- array(0, c( length(filesList), 8 ) )

column_names <- c("Participant_number", sapply(1:10, function(n) paste(n, "Clusters", sep = "_")))
LH_SH <- data.frame(matrix(ncol = length(column_names), nrow = 0))
RH_SH <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(LH_SH) <- column_names
colnames(RH_SH) <- column_names

for (nFile in 1:length(filesList) ) { #nFile <- 1; 
  for (nHemi in 1:length(hemispheres) ) { #nHemi <- 1;
    
    print('........................')
    print('.....main folder........')
    print('........................')
    setwd( mainDir ); print( getwd() ); print( sprintf('nFile: %d',nFile) )
    print('........................')
    print('........................')
    nFile <- nFile
    fileNameIn <- filesList[nFile] 
    participantName <- strsplit( fileNameIn, '_' )[[1]][1]
    partipantNumber <- as.numeric( strsplit( participantName, '-' )[[1]][2] )
    fileHemisphere <- strsplit( fileNameIn, '_' )[[1]][2]
    fileHemisphere <- strsplit( fileHemisphere, '[.]' )[[1]][1]
    selectedHemi <- fileHemisphere
    load( file=fileNameIn )
    
    
    dfOutSel <- dfOut[ dfOut$modality=='T1W', ]
    profilesArray <- dfOutSel[,8:15]
    roiCoords <- t( array( as.matrix( dfOutSel[,5:7] ), c( dim(dfOutSel)[1], 3 ) ) )
    
    cleanData <- function( idxIn ) { #store all the intercepts across depth, take the average and add the average intercept per profile
      #idxIn <- 1
      modCovar <- lm( profilesArray[,idxIn] ~ poly( dfOutSel$thickness, 1 ) * poly( dfOutSel$curvature, 1 ) )
      summary( modCovar )
      #plot( modCovar )
      dataClean <- round( modCovar$resid + coefficients( modCovar )[1], 2 )
      #dataClean <- round( modCovar$resid + 3000, 2 )
      #dataClean <- round( modCovar$resid, 2 )
      return( dataClean )
    }
    #profilesArray <- profilesArray
    profilesArray <- sapply( 1:8, FUN = cleanData, simplify = TRUE  )
    
    # check profiles:
    print('........................')
    print('.....check profiles.....')
    print('........................')
    #which( apply( abs( apply( profilesArray, 1, diff ) ), 2, sum ) < 0.00001 )
    absDifferences <- abs( apply( profilesArray, 1, diff ) )
    smallAbsDifferences <- apply( absDifferences < 0.0000001, 2, sum )
    if (selectedHemi=='RH') { idxKeepProfiles <- smallAbsDifferences < 2 & dfOutSel$x < 95 }
    if (selectedHemi=='LH') { idxKeepProfiles <- smallAbsDifferences < 2 & dfOutSel$x > 140 }
    #idxKeepProfiles <- smallAbsDifferences < 2
    table( smallAbsDifferences  )
    
    # filter profiles
    print('...............................')
    print('.....filter/scale profiles.....')
    print('...............................')
    roiCoords <- roiCoords[,idxKeepProfiles] 
    profilesArray <- profilesArray[idxKeepProfiles,]
    dfOutSel <- dfOutSel[idxKeepProfiles,] 
    
    # macro 'voxels'
    print('...............................')
    print('.....macro voxels..............')
    print('...............................')
    dataMat <- data.frame( t( roiCoords ) )
    names( dataMat ) <- c('xc','yc','zc')
    kMacroVoxels <- kmeans(x=dataMat, centers=800, nstart=10)
    #kMacroVoxels$cluster <- seq(1,dim(dataMat)[1])
    str( dataMat )
    p <- plot_ly( dataMat, x = ~xc, y = ~yc, z = ~zc, color = kMacroVoxels$cluster, size=1 ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    aggregateProfiles <- aggregate( profilesArray, by=list(kMacroVoxels$cluster), median ) # average profiles based on clustering
    aggregateProfilesMat <- as.matrix( aggregateProfiles[,2:dim(aggregateProfiles)[2]] )
    aggregateCoords <- aggregate( t(roiCoords), by=list(kMacroVoxels$cluster), median ) # average coordinates based on clustering
    aggregateCurvature <- aggregate( dfOutSel$curvature, by=list(kMacroVoxels$cluster), median ) # average curvature based on clustering
    aggregateThickness <- aggregate( dfOutSel$thickness, by=list(kMacroVoxels$cluster), median ) # average curvature based on clustering
    aggregateCoords$curvature <- aggregateCurvature[,2]
    aggregateCoords$thickness <- aggregateThickness[,2]
    
    names( aggregateCoords ) <- c('cluster','xc','yc','zc', 'curvature', 'thickness')
    str( aggregateCoords )
    p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~cluster, size=1 ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~curvature, size=1 ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = ~thickness, size=1 ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    # macro 'voxels'
    print('...............................')
    print('...model macro voxels..........')
    print('...............................')
    modelFunction <- function( idx, polyPar ) {
      singleProfile <- aggregateProfilesMat[idx, ]
      xPredictor <- seq(0,1,length.out = length(singleProfile))
      mod.lm <- lm( singleProfile ~ poly( xPredictor, polyPar ) )
      sum.mod.lm <- summary( mod.lm )
      return( as.numeric( c( coefficients( mod.lm ), sum.mod.lm$r.squared ) ) )
    }
    modelMacroVoxels <- sapply( 1:dim(aggregateProfilesMat)[1], modelFunction, polyPar=1 )
    
    str( aggregateCoords )
    p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = modelMacroVoxels[1,], size=1, colors='YlOrRd' ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    str( aggregateCoords )
    p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = modelMacroVoxels[2,], size=1, colors='YlOrRd' ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    colorMacroVoxel <- as.numeric( modelMacroVoxels[2,] < 0 )
    p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = colorMacroVoxel, size=1 ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    #####################################################################
    keepMacroVoxel <-  modelMacroVoxels[2,] < 0 #rep( TRUE, length( modelMacroVoxels[2,] ) ) # # 
    modelMacroVoxels_01 <- modelMacroVoxels[ , keepMacroVoxel ]
    aggregateProfiles_01 <- aggregateProfiles[ keepMacroVoxel, ]
    aggregateProfilesMat_01 <- aggregateProfilesMat[ keepMacroVoxel, ]
    aggregateCoords_01 <- aggregateCoords[ keepMacroVoxel,  ]
    #####################################################################
    
    p <- plot_ly( aggregateCoords_01, x = ~xc, y = ~yc, z = ~zc, size=1 ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    x11()
    par(mfrow=c(2,2))
    plot( modelMacroVoxels_01[2,] ~ modelMacroVoxels_01[1,], bty='n', las=1, xlab='intercept', ylab='slope', cex=1, pch=19, col=rgb(0.5,0.5,0.5,0.5) )
    abline( lm( modelMacroVoxels_01[2,] ~ modelMacroVoxels_01[1,] ), lwd=2, col='red' )
    hist( modelMacroVoxels_01[1,] )  
    hist( modelMacroVoxels_01[2,] )  
    
    # Determine number of clusters and some plots
    print('......................................')
    print('.....determine number of clusters.....')
    print('......................................')
    
    #kUmap <- kmeans( aggregateProfilesMat_01, centers = 2, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )
    kUmap <- kmeans( aggregateProfilesMat_01, centers = 2, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )
    
        #Silhouette Plot
    silhouette_scores <- fviz_nbclust(aggregateProfilesMat_01, kmeans, method = "silhouette", k.max = 10, nstart = 25) 
    sh_score = c(nFile, silhouette_scores$data$y)
    print(nHemi)
    if(nHemi == 1){LH_SH = rbind(LH_SH, setNames(as.data.frame(t(sh_score)), column_names))} else {RH_SH = rbind(RH_SH, setNames(as.data.frame(t(sh_score)), column_names))}
  
    
    #arrange clusters
    #pos01 <- apply( aggregateCoords_01[ kUmap$cluster==1, c(2:4) ], 2, mean )
    #pos02 <- apply( aggregateCoords_01[ kUmap$cluster==2, c(2:4) ], 2, mean )
    #if ( pos01[2] > pos02[2] ) {
    #  app <- kUmap$cluster
    #  kUmap$cluster[ app==1 ] <- 1
    #  kUmap$cluster[ app==2 ] <- 2
    #} else {
    #  app <- kUmap$cluster
    #  kUmap$cluster[ app==2 ] <- 1
    #  kUmap$cluster[ app==1 ] <- 2
    #}
    
    # arrange clusters
    meanProfile01 <- mean( apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean ) )
    meanProfile02 <- mean( apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean ) )
    if ( meanProfile01 < meanProfile02 ) {
      app <- kUmap$cluster
      kUmap$cluster[ app==1 ] <- 1
      kUmap$cluster[ app==2 ] <- 2
    } else {
      app <- kUmap$cluster
      kUmap$cluster[ app==2 ] <- 1
      kUmap$cluster[ app==1 ] <- 2
    }
    
    x11( width=4, height=4 )
    par(mfrow=c(1,1))
    profile01 <- apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean )
    profile01SE_P <- apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean ) + apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, sd )/2
    profile01SE_M <- apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, mean ) - apply( aggregateProfilesMat_01[ kUmap$cluster==1, ], 2, sd )/2
    profile02 <- apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean )
    profile02SE_P <- apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean ) + apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, sd )/2
    profile02SE_M <- apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, mean ) - apply( aggregateProfilesMat_01[ kUmap$cluster==2, ], 2, sd )/2
    plot( profile01~seq(0,1,length.out = 8), bty='n', las=1, col='blue', xlab='cortical depth', ylab='T1w intensity (A.U.)', type='l', lwd=4,
          ylim=c( min( c( profile01SE_M, profile02SE_M ) ), max( c( profile01SE_P, profile02SE_P ) ) ), axes=FALSE, cex.lab=1.5  )
    axis(1,seq(0,1,0.5), cex.axis=1.5); axis(2, seq( round(min( c( profile01SE_M, profile02SE_M ) ),0), round(max( c( profile01SE_P, profile02SE_P ) ),0), length.out = 3 ), cex.axis=1.5 )
    lines( seq(0,1,length.out = 8), profile02, col='darkorange', lwd=4 )
    lines( seq(0,1,length.out = 8), profile01SE_M, col='blue', lwd=2, lty=2 )
    lines( seq(0,1,length.out = 8), profile01SE_P, col='blue', lwd=2, lty=2 )
    lines( seq(0,1,length.out = 8), profile02SE_M, col='darkorange', lwd=2, lty=2 )
    lines( seq(0,1,length.out = 8), profile02SE_P, col='darkorange', lwd=2, lty=2 )
    
    p <- plot_ly( aggregateCoords_01, x = ~xc, y = ~yc, z = ~zc, color = kUmap$cluster, size=1, colors=c('blue','darkorange') ) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'xc'),
                          yaxis = list(title = 'yc'),
                          zaxis = list(title = 'zc')))
    p
    
    
    kmeansAIC = function(fit){
      m = ncol(fit$centers)
      n = length(fit$cluster)
      k = nrow(fit$centers)
      D = fit$tot.withinss
      return( data.frame(AIC = D + 2*m*k,
                         BIC = D + log(n)*m*k) )
    }
    
    nClusterToTest <- 8
    totwss <- rep(0, nClusterToTest)
    storeBIC <- array( 0, c( 1, nClusterToTest ) )
    storeAIC <- storeBIC 
    storeR2 <- storeBIC
    
    n <- dim( aggregateProfilesMat_01 )[1]
    k <- 1
    
    for (fitCounter in 1:nClusterToTest) {
      kFit <- kmeans( aggregateProfilesMat_01, centers = fitCounter, nstart = 25,  algorithm = c("Hartigan-Wong" ), iter.max=1500 )
      #kFit <- kmeans( selSf, fitCounter )
      storeBicTemp <- kmeansAIC( kFit )
      storeBIC[k,fitCounter] <- storeBicTemp$BIC
      storeAIC[k,fitCounter] <- storeBicTemp$AIC
      totwss[ fitCounter ] <- kFit$tot.withinss
    }
    
    rsq = 1-(totwss*(n-1))/(totwss[1]*(n-seq(1,nClusterToTest)))
    storeR2[k,] <- rsq  
    
    
    #matplot( t( storeBIC ) )
    #matplot( t( storeAIC ) )
    #matplot( t( storeR2 ) )
    
    x11(width=3.75,height=4)
    
    matplot( t( storeBIC ), pch=3, cex=1.25, las=1, bty='n', ylab='BIC score (sf)', xlab='n. clusters', cex.lab=1.25, cex.axis=1.15  )
    lines( apply( t(storeBIC), 1, median ), lwd=3, lty=2  )
    
    #### bring data back into anatomy standard space: ####
    
    length(kMacroVoxels$cluster); dim(dataMat)
    length(keepMacroVoxel); 
    length(kUmap$cluster)
    
    clusterDataFull <- rep( 9999, length(keepMacroVoxel) )
    clusterDataFull[ keepMacroVoxel==TRUE ] <- kUmap$cluster
    clusterDataFull[ keepMacroVoxel==FALSE ] <- 0
    clusterDataFull_toVoxels <- clusterDataFull[ kMacroVoxels$cluster ]
    
    if ( partipantNumber <= 50 ) {
      setwd( freesurferDirPart01 ); getwd()
    }else{
      setwd( freesurferDirPart02 ); getwd()
    }
    setwd( participantName ); getwd()
    setwd('ses-1')
    setwd('anat')
    setwd( 'AHEAD_test' ); getwd()
    setwd( 'SUMA' ); 
    
    print('...........................')
    print('.....current folder........')
    print('...........................')
    print( getwd() )
    print('...........................')
    print('...........................')
    referenceAnatomy <- dir( pattern='*_mod-t1w_orient-std_brain.nii.gz' )[1]
    anatFile <- read.AFNI( referenceAnatomy )
    anatVolume <- anatFile$brk[,,,1]
    anatVolumeDim <- dim( anatVolume )
    emptyAnatVolume <- array( 0, anatVolumeDim )
    coordsCluster01 <- dataMat[ clusterDataFull_toVoxels==1, ] 
    coordsCluster02 <- dataMat[ clusterDataFull_toVoxels==2, ] 
    linIndexCluster01 <- linearIndexFromCoordinate( t( as.matrix( coordsCluster01 ) ), anatVolumeDim )
    linIndexCluster02 <- linearIndexFromCoordinate( t( as.matrix( coordsCluster02 ) ), anatVolumeDim )
    emptyAnatVolume[ linIndexCluster01 ] <- 1
    emptyAnatVolume[ linIndexCluster02 ] <- 2
    
    if (file.exists('_ttt_clusters.nii.gz')) { system('rm _ttt_clusters.nii.gz') }
    write.AFNI( '_ttt_clusters.nii.gz', brk = emptyAnatVolume,
                origin = anatFile$origin, orient = anatFile$orient,
                defhead = anatFile$NI_head )
    
    if ( selectedHemi=='LH') {
      outName1D <- sprintf('insulaClusters_%s_lh.1D.dset',saveName)
      if (file.exists(outName1D)) { system( sprintf('rm %s', outName1D ) ) }
      filenameIn <- '_ttt_clusters.nii.gz'
      instr <- paste( '3dVol2Surf',                 
                      '-spec std.141.AHEAD_test_lh.spec',
                      '-surf_A std.141.lh.white.gii',
                      '-surf_B std.141.lh.pial.gii',
                      sprintf('-sv %s', filenameIn ),
                      sprintf('-grid_parent %s', filenameIn ),
                      '-map_func nzmode', 
                      '-f_steps 10',
                      sprintf('-out_1D %s', outName1D) ) #'-f_steps 10',
      system( instr )
    }
    if ( selectedHemi=='RH') {
      outName1D <- sprintf('insulaClusters_%s_rh.1D.dset',saveName)
      if (file.exists(outName1D)) { system( sprintf('rm %s', outName1D ) ) }
      filenameIn <- '_ttt_clusters.nii.gz'
      instr <- paste( '3dVol2Surf',                 
                      '-spec std.141.AHEAD_test_rh.spec',
                      '-surf_A std.141.rh.white.gii',
                      '-surf_B std.141.rh.pial.gii',
                      sprintf('-sv %s', filenameIn ),
                      sprintf('-grid_parent %s', filenameIn ),
                      '-map_func nzmode', 
                      '-f_steps 10',
                      sprintf('-out_1D %s', outName1D) ) #'-f_steps 10',
      system( instr )
    }
    if (file.exists('_ttt_clusters.nii.gz')) { system('rm _ttt_clusters.nii.gz') }
    graphics.off()
  }
}

## silhouette Plot ----
col_wise_stats = function(df){
  df = df[,-1]  
  # Column-wise mean
  column_means <- colMeans(df, na.rm = TRUE)
  # Column-wise SD
  column_sds <- apply(df, 2, sd, na.rm = TRUE)
  # Column-wise SEM
  column_sems <- column_sds / sqrt(nrow(df))
  # Combine all 
  summary_df <- data.frame(
    Mean = column_means,
    SD = column_sds,
    SEM = column_sems)
  return(summary_df)
}
stats_LH = col_wise_stats(LH_SH)
stats_RH = col_wise_stats(RH_SH)

#Combine
stats_LH$Hemisphere <- 'LH'
stats_RH$Hemisphere <- 'RH'
stats_LH$Cluster_no <- 1:nrow(stats_LH)
stats_RH$Cluster_no <- 1:nrow(stats_RH)
combined_stats <- rbind(stats_LH, stats_RH)
highest_point = combined_stats %>% arrange(desc(Mean)) %>% select(Mean,Cluster_no) 
highest_point = highest_point[1,]


#Plot for just mean
ggplot(combined_stats, aes(x = Cluster_no, y = Mean, color = Hemisphere, group = Hemisphere)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD, color = Hemisphere), width = 0.2) +  # Coloring error bars by Hemisphere
  theme_minimal() +
  labs(x = "Number of Clusters", y = "Mean Silhouette Score", color = "Hemisphere") +
  scale_x_continuous(breaks = 1:10) +
  theme(plot.title = element_blank()) +
  scale_y_continuous(breaks = c(0.0, 0.4)) +
  annotate("text", x = highest_point$Cluster_no, y = highest_point$Mean, 
           label = "Optimal No of Clusters", vjust = 35, hjust = -.05, color = "black") +
  geom_vline(xintercept = highest_point$Cluster_no, linetype = "dotted", color = "black", size = 1)

