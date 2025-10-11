rm( list=ls() ); gc();
graphics.off();

mainDir <- '/analyse/Project0165/Fracasso_Anatomy/AHEAD_database/10007840/results/insula'
setwd( mainDir )

set.seed(1);
library( plotly )
library( ica )
afniDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniAtlasDir <- Sys.getenv(x='AFNI_ATLASDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
source( sprintf( '%s/AFNIio.R', afniDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )

setwd( mainDir )

filesList <- dir( pattern = '*.RData' )
storeProfiles01 <- array(0, c( length(filesList), 8 ) )
storeProfiles02 <- array(0, c( length(filesList), 8 ) )

#for (nFile in 1:length(filesList) ) {
  
  nFile <- 3
  load( file=filesList[nFile] )
  
  dfOutSel <- dfOut[ dfOut$modality=='T1W', ]
  profilesArray <- dfOutSel[,8:15]
  roiCoords <- t( array( as.matrix( dfOutSel[,5:7] ), c( dim(dfOutSel)[1], 3 ) ) )
  
  cleanData <- function( idxIn ) {
    #idxIn <- 1
    modCovar <- lm( profilesArray[,idxIn] ~ poly( dfOutSel$thickness, 3 ) * poly( dfOutSel$curvature, 3 ) )
    summary( modCovar )
    dataClean <- round( modCovar$resid + coefficients( modCovar )[1], 2 )
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
  idxKeepProfiles <- smallAbsDifferences < 2
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
  kMacroVoxels <- kmeans(x=dataMat, centers=2000, nstart=10) 
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
  p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = modelMacroVoxels[2,], size=1 ) %>%
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
  keepMacroVoxel <-  modelMacroVoxels[2,] < 0 # rep( TRUE, length( modelMacroVoxels[2,] ) ) #
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
  
  #arrange clusters
  pos01 <- apply( aggregateCoords_01[ kUmap$cluster==1, c(2:4) ], 2, mean )
  pos02 <- apply( aggregateCoords_01[ kUmap$cluster==2, c(2:4) ], 2, mean )
  if ( pos01[2] < pos02[2] ) {
   app <- kUmap$cluster
   kUmap$cluster[ app==1 ] <- 1
   kUmap$cluster[ app==2 ] <- 2
  } else {
   app <- kUmap$cluster
   kUmap$cluster[ app==2 ] <- 1
   kUmap$cluster[ app==1 ] <- 2
  }
  
  p <- plot_ly( aggregateCoords_01, x = ~xc, y = ~yc, z = ~zc, color = kUmap$cluster, size=1, colors=c('blue','darkorange') ) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = 'xc'),
                        yaxis = list(title = 'yc'),
                        zaxis = list(title = 'zc')))
  p
  
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
  
  
  # save data on r file about clusterization
  # bring back data into single subject space and then into atlas for cross validation
  #kMacroVoxels$cluster == 1 
  
  
  # filename <- strsplit( filesList[nFile], '[.]' )[[1]][1]
  # dev.copy2pdf( file=sprintf('%s.pdf',filename), width=4, height=4)
  # dev.off()
  # 
  # storeProfiles01[ 
  #nFile, ] <- profile01
  # storeProfiles02[ nFile, ] <- profile02
#}


# x11( width=4, height=4 )
# par(mfrow=c(1,1))
# profile01 <- apply( storeProfiles01, 2, mean )
# profile01SE_P <- apply( storeProfiles01, 2, mean ) + apply( storeProfiles01, 2, sd )/2
# profile01SE_M <- apply( storeProfiles01, 2, mean ) - apply( storeProfiles01, 2, sd )/2
# profile02 <- apply( storeProfiles02, 2, mean )
# profile02SE_P <- apply( storeProfiles02, 2, mean ) + apply( storeProfiles02, 2, sd )/2
# profile02SE_M <- apply( storeProfiles02, 2, mean ) - apply( storeProfiles02, 2, sd )/2
# plot( profile01~seq(0,1,length.out = 8), bty='n', las=1, col='blue', xlab='cortical depth', ylab='T1w intensity (A.U.)', type='l', lwd=4,
#       ylim=c( min( c( profile01SE_M, profile02SE_M ) ), max( c( profile01SE_P, profile02SE_P ) ) ), axes=FALSE, cex.lab=1.5  )
# axis(1,seq(0,1,0.5), cex.axis=1.5); axis(2, seq( round(min( c( profile01SE_M, profile02SE_M ) ),0), round(max( c( profile01SE_P, profile02SE_P ) ),0), length.out = 3 ), cex.axis=1.5 )
# lines( seq(0,1,length.out = 8), profile02, col='darkorange', lwd=4 )
# lines( seq(0,1,length.out = 8), profile01SE_M, col='blue', lwd=2, lty=2 )
# lines( seq(0,1,length.out = 8), profile01SE_P, col='blue', lwd=2, lty=2 )
# lines( seq(0,1,length.out = 8), profile02SE_M, col='darkorange', lwd=2, lty=2 )
# lines( seq(0,1,length.out = 8), profile02SE_P, col='darkorange', lwd=2, lty=2 )
# 
# filename <- strsplit( filesList[nFile], '[.]' )[[1]][1]
# dev.copy2pdf( file=sprintf('averageProfiles_across_participants.pdf'), width=4, height=4)
# dev.off()
# 
# 

  #library( pvclust )
  #d <- dist( aggregateProfilesMat_01, method='euclidean' )
  #fit <- hclust(d, method='ward.D')
  #plot(fit)
  #groups <- cutree(fit, k=3)
  #rect.hclust(fit, k=3, border="red") 
  #kUmap <- 0
  #kUmap$cluster <- groups
  
  #fit <- pvclust( aggregateProfilesMat_01, method.hclust='ward', method.dist='euclidean' )
  #plot( fit )
  #kUmap <- kmeans( aggregateProfilesMat_01, centers = 2, nstart = 10 )
  

  #kfuz <- Fclust( aggregateProfilesMat, 3, type='polynomial')
  #par(mfrow=c(1,1))
  #p <- plot_ly( aggregateCoords, x = ~xc, y = ~yc, z = ~zc, color = kfuz$U[,3], sizes=1, colors=c('blue','darkorange') ) %>%
  #  add_markers() %>%
  #  layout(scene = list(xaxis = list(title = 'xc'),
  #                      yaxis = list(title = 'yc'),
  #                      zaxis = list(title = 'zc')))
  #p
  

