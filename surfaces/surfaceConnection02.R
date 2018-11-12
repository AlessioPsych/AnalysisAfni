rm( list=ls() )
library(abind)

participantArray <- c( 'SD', 'BH', 'VK', 'AF', 'AF01', 'RB' ) #SD BH VK AF AF01 RB

for ( nParticipant in 1:length(participantArray) ) {
#for ( nParticipant in 1:1 ) {
  
  graphics.off()
  
  #nParticipant <- 6
  
  participant <- participantArray[ nParticipant ]
  
  if (participant=='BH') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V5029leftright_Ben_ph/results_CBS'
    statDir <- 'statsSingleShot_mi+orig_interp_folder'
    roiDirRight <- 'rightV1_folder'
    roiDirLeft <- 'leftV1_folder'
    curvatureDir <- 'curvature_interp_folder_smoothed_8'
    thicknessDir <- 'thickness_interp_folder'
    amplitudeDir <- 'amplitudeSingleShot_mi+orig_interp_folder'
    anatomyDir <- 'anatomy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_1-2-1'
    angleShift <- 45
  }
  
  if (participant=='VK') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2699leftright_Vin/results_CBS'
    statDir <- 'statsSingleShot_mi+orig_interp_folder'
    roiDirRight <- 'rightV1_folder'
    roiDirLeft <- 'leftV1_folder'
    curvatureDir <- 'curvatureCopy_interp_folder_smoothed_8'
    thicknessDir <- 'thicknessCopy_interp_folder'
    amplitudeDir <- 'amplitudeSingleShot_mi+orig_interp_folder'
    anatomyDir <- 'anatomyCopy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_00-1'
    angleShift <- -10
  }
  
  if (participant=='AF') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2678leftright_Ale/results_CBS'
    statDir <- 'statsSingleShotEPI_add_lpc+orig_interp_folder'
    roiDirRight <- 'rightV1_folder'
    roiDirLeft <- 'leftV1_folder'
    curvatureDir <- 'curvature_interp_folder_smoothed_8'
    thicknessDir <- 'thickness_interp_folder'
    amplitudeDir <- 'amplitudeSingleShotEPI_add_lpc+orig_interp_folder'
    anatomyDir <- 'anatomy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_-1-1-1'
    angleShift <- 30
  }
  
  if (participant=='AF01') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V4847leftright_Ale_ph/results_CBS'
    statDir <- 'statsSingleShotEPI_lpc_left_add_interp_folder'
    roiDirRight <- 'rightV1_folder'
    roiDirLeft <- 'leftV1_folder'
    curvatureDir <- 'curvature_interp_folder_smoothed_8'
    thicknessDir <- 'thickness_interp_folder'
    amplitudeDir <- 'amplitudeSingleShotEPI_lpc_left_add_interp_folder'
    anatomyDir <- 'anatomy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_222'
    angleShift <- -10
  }
  
  
  if (participant=='SD') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2676leftright_Ser/results_CBS'  
    statDir <- 'statsSingleShotEPI_add_mi+orig_interp_folder'
    roiDirRight <- 'rightV1_folder'
    roiDirLeft <- 'leftV1_folder'  
    curvatureDir <- 'curvatureCopy_interp_folder_smoothed_8'
    thicknessDir <- 'thicknessCopy_interp_folder'  
    amplitudeDir <- 'amplitudeSingleShotEPI_add_mi+orig_interp_folder'
    anatomyDir <- 'anatomyCopy_interp_folder'
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_-10-1'
    angleShift <- -10
  }
  
  if (participant=='RB') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2677leftright_Ric/resultsCBS'  
    statDir <- 'statsSingleShotEPI_add_lpc_01_interp_folder'
    roiDirRight <- 'rightV1_folder'
    roiDirLeft <- 'leftV1_folder'  
    curvatureDir <- 'curvatureCopy_interp_folder_smoothed_8'
    thicknessDir <- 'thicknessCopy_interp_folder'
    amplitudeDir <- 'amplitudeSingleShotEPI_add_lpc_01_interp_folder'
    anatomyDir <- 'anatomyCopy_interp_folder'
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_-2-2-2'
    angleShift <- 80
  }
  
  
  source('/home/alessiofracasso/Dropbox/analysisAfni/load1DData.R')
  source('/home/alessiofracasso/Dropbox/analysisAfni/generateProfiles.R')
  source('/home/alessiofracasso/Dropbox/analysisAfni/scaleData.R')
  
  
  ## load ROI right files
  fileNumber <- c(1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20)
  directory <- sprintf('%s/%s', mainDir, roiDirRight )
  surfValListRight <- load1DData( directory, fileNumber, 0, '^roiSurface.*.1D.roi' )
  
  ## load ROI left files
  directory <- sprintf('%s/%s', mainDir, roiDirLeft )
  surfValListLeft <- load1DData( directory, fileNumber, 0, '^roiSurface.*.1D.roi' )
  
  ## load stat files
  fileNumber <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19)
  directory <- sprintf('%s/%s', mainDir, statDir )
  mapValListStat <- load1DData( directory, fileNumber, 5, '*.1D' )
  
  ## load curv files
  directory <- sprintf('%s/%s', mainDir, curvatureDir )
  mapValListCurv <- load1DData( directory, fileNumber, 5, '*.1D.dset'  )
  
  ## load thickness files
  directory <- sprintf('%s/%s', mainDir, thicknessDir )
  mapValListThick <- load1DData( directory, fileNumber, 5, '*surfval.1D'  )
  
  ## load amplitude files
  directory <- sprintf('%s/%s', mainDir, amplitudeDir )
  mapValListAmp <- load1DData( directory, fileNumber, 5, '*.1D'  )
  
  ## load anatomy files
  directory <- sprintf('%s/%s', mainDir, anatomyDir )
  mapValListAnat <- load1DData( directory, fileNumber, 5, '*surfval.1D'  )
  
  ## load B0angle files
  directory <- sprintf('%s/%s', mainDir, B0angleDir )
  mapValListB0Angle <- load1DData( directory, fileNumber, 5, '*.1D.dset'  )
  
  pRoiRightPos <- generateProfiles( surfValListRight, mapValListStat, 8 )
  pRoiRightNeg <- generateProfiles( surfValListRight, mapValListStat, 11 )
  pRoiRightCurv <- generateProfiles( surfValListRight, mapValListCurv, 7 )
  pRoiRightThick <- generateProfiles( surfValListRight, mapValListThick, 7 )
  pRoiRightAmp <- generateProfiles( surfValListRight, mapValListAmp, 7 )
  pRoiRightTstat <- generateProfiles( surfValListRight, mapValListStat, 15 )
  pRoiRightAnat <- generateProfiles( surfValListRight, mapValListAnat, 7 )
  pRoiRightTstatSinglePos <- generateProfiles( surfValListRight, mapValListStat, 9 )
  pRoiRightB0angle <- generateProfiles( surfValListRight, mapValListB0Angle, 2 )
  
  pRoiLeftPos <- generateProfiles( surfValListLeft, mapValListStat, 11 )
  pRoiLeftNeg <- generateProfiles( surfValListLeft, mapValListStat, 8 )
  pRoiLeftCurv <- generateProfiles( surfValListLeft, mapValListCurv, 7 )
  pRoiLeftThick <- generateProfiles( surfValListLeft, mapValListThick, 7 )
  pRoiLeftAmp <- generateProfiles( surfValListLeft, mapValListAmp, 7 )
  pRoiLeftTstat <- generateProfiles( surfValListLeft, mapValListStat, 15 )
  pRoiLeftAnat <- generateProfiles( surfValListLeft, mapValListAnat, 7 )
  pRoiLeftTstatSinglePos <- generateProfiles( surfValListLeft, mapValListStat, 12 )
  pRoiLeftB0angle <- generateProfiles( surfValListLeft, mapValListB0Angle, 2 )
  
  intensityRightStruct <- list(0)
  intensityLeftStruct <- list(0)
  
  intensityRightStruct$intensity <- abind( pRoiRightPos$intensity,
                                           pRoiRightNeg$intensity,
                                           pRoiRightCurv$intensity,
                                           pRoiRightThick$intensity,
                                           pRoiRightAmp$intensity,
                                           pRoiRightTstat$intensity,
                                           pRoiRightAnat$intensity,
                                           pRoiRightTstatSinglePos$intensity,
                                           pRoiRightB0angle$intensity,
                                           along=3)
  intensityRightStruct$unconnectedPoints <- pRoiRightCurv$unconnectedPoints
  
  intensityLeftStruct$intensity <- abind( pRoiLeftPos$intensity,
                                          pRoiLeftNeg$intensity,
                                          pRoiLeftCurv$intensity,
                                          pRoiLeftThick$intensity,
                                          pRoiLeftAmp$intensity,
                                          pRoiLeftTstat$intensity,
                                          pRoiLeftAnat$intensity,
                                          pRoiLeftTstatSinglePos$intensity,
                                          pRoiLeftB0angle$intensity,
                                          along=3)
  intensityLeftStruct$unconnectedPoints <- pRoiLeftCurv$unconnectedPoints
  
  #intensityRightStruct <- generateProfiles01( surfValListRight, mapValListStat, mapValListCurv, mapValListThick, mapValListAmp, 8, 11 ) 
  #intensityLeftStruct <- generateProfiles01( surfValListLeft, mapValListStat, mapValListCurv, mapValListThick, mapValListAmp, 11, 8 ) 
  
  intensityRight <- intensityRightStruct$intensity
  intensityLeft <- intensityLeftStruct$intensity
  intensity <- abind( intensityRight, intensityLeft, along=1 )
  
  plot( apply( intensity[,,1], 2, mean ) )
  
  storeCoeffs <- array(1,c(dim(intensity)[1],9))
  storeRes <- array( 1, dim(intensity[,,1]) )
  xVar <- seq(3,16)
  for (k in 1:dim(intensity)[1] ) {  
    y <- intensity[k,xVar,1]
    yNeg <- intensity[k,xVar,2]
    yAmp <- intensity[k,xVar,5]
    yCurv <- intensity[k,xVar,3]
    yThick <- intensity[k,xVar,4]
    yAnatInt <- intensity[k,xVar,7]
    #statProfile <- abs( min(  intensity[k,xVar,6]  ) )
    #meanProfilePositive <- mean( intensity[k,xVar,1] )
    yFilt <- y[ !is.na(y) ]
    xFilt <- seq(0,1,length.out=length( xVar[ !is.na(y) ] ) )
    yAmpFilt <- yAmp[  !is.na(y) ]
    yCurvFilt <- yCurv[  !is.na(y) ]
    yThickFilt <- yThick[  !is.na(y) ]
    if ( ( sum( y==rep(0,length(xVar)) ) < 4 ) | 
           ( sum( yNeg==rep(0,length(xVar)) ) < 4 ) ) {
      
      sCor <- cor( yFilt, xFilt, method = c('spearman') )
      mod <- lm( yFilt ~ xFilt )    
      sMod <- summary( mod )
      
      modAmp <- lm( yAmpFilt ~ xFilt )
      sCorAmp <- cor( yAmpFilt, xFilt, method = c('spearman') )
      storeCoeffs[k,] <- c( mod$coefficients[1:2], 0, 0, sMod$r.squared, modAmp$coefficients[1:2], sCor, sCorAmp )
      #mod <- lm( yFilt ~ xFilt )
      #sMod <- summary( mod )
      #storeCoeffs[k,] <- c( mod$coefficients, sMod$r.squared, 0, 0, 0, 0 )
      
      #mod <- lm( yFilt ~ xFilt + yCurvFilt )
      #sMod <- summary( mod )
      #storeCoeffs[k,] <- c( mod$coefficients, sMod$r.squared, 0, 0, 0 )
      storeRes[k,xVar] <- mod$residuals
    }
    
  }
  
  plot( storeCoeffs[,2]~storeCoeffs[,1], cex=0.1, ylim=c(-10,10) ) 
  plot( storeCoeffs[,7]~storeCoeffs[,6], cex=0.1 ) 
  #rankProfiles <- storeCoeffs[,6]
  #rankProfiles <- storeCoeffs[,5] * storeCoeffs[,2]
  rankProfiles <- storeCoeffs[,8]
  sum(is.na(rankProfiles))
  
  
  delRows <- c( intensityRightStruct$unconnectedPoints, intensityLeftStruct$unconnectedPoints )
  filtRows <- which( delRows==0 )
  if ( length( filtRows ) > 0 ) { delRows <- delRows[-filtRows] }
  
  rm( naEl )
  #if (length( delRows ) > 0) { naEl <- unique( c( which(is.na(rankProfiles)), which( storeCoeffs[,7]==1 ), delRows ) ) }
  #if (length( delRows ) == 0) { naEl <- unique( c( which(is.na(rankProfiles)), which( storeCoeffs[,7]==1 ) ) ) }
  if (length( delRows ) > 0) { naEl <- unique( c( which(is.na(rankProfiles)), delRows ) ) }
  if (length( delRows ) == 0) { naEl <- unique( c( which(is.na(rankProfiles)) ) ) }
  #naEl <- 0
  
  surfValRight <- data.frame( surfValListRight[[10]] )
  surfValLeft <- data.frame( surfValListLeft[[10]] )
  surfVal <- c( surfValRight[,1], surfValLeft[,1] )
  
  dim(intensity)[1]
  length(rankProfiles)
  length(surfVal)
  
  boldMedian <- apply( intensity[,xVar,1], 1, median )
  boldMin <- apply( intensity[,xVar,1], 1, min )
  #boldFilt <- which( boldMedian < quantile( boldMedian, c( 0.025 ) ) | boldMedian > quantile( boldMedian, c( 0.975 ) ) )
  boldFilt <- which( boldMin <0 | boldMedian > quantile( boldMedian, c( 0.975 ) ) )
  curvMedian <- apply( intensity[ ,xVar, 3 ], 1, median  )
  curvFilt <- which( curvMedian < -0.7 | curvMedian > 0.7 )
  thickMedian <- apply( intensity[ ,xVar, 4 ], 1, median  )
  thickFilt <- which( thickMedian<0.5 | thickMedian>quantile(thickMedian,0.975)  )
  anatMin <- apply( intensity[ ,xVar, 7 ], 1, min  )
  anatFilt <- which( anatMin <= 10000 | anatMin > 45000 )
  
  tZero <- apply( intensity[,xVar,8], 1, median )
  tMeanFilt <- which( tZero<1.945 )
  
  #tMean <- apply( intensity[ , seq(8,12), 8 ], 1, min )
  #tMeanFilt2 <- which( tMean<1 )
  
  ampBeg <- apply( intensity[, xVar[ 1 : 3 ], 5 ], 1, mean ) 
  ampEnd <- apply( intensity[, xVar[ ( length(xVar)-2 ) : length(xVar) ], 5 ], 1, mean )
  #ampLinearFitFilt <- which( ampBeg<ampEnd*0.7 
  #                           | ampBeg < quantile(ampBeg, c(0.025) ) | ampBeg > quantile(ampBeg, c(0.975) ) 
  #                           | ampEnd < quantile(ampEnd, c(0.025) ) | ampEnd > quantile(ampEnd, c(0.975) ) ) 
  
  ampLinearFitFilt <- which( storeCoeffs[,9] < 0 | is.na( storeCoeffs[,9] ) ) 
  
  
  
  
  #naElOut <- unique( c( ampLinearFitFilt, naEl, anatFilt, tMeanFilt, thickFilt, curvFilt ) )
  naElOut <- unique( c( naEl, anatFilt, thickFilt, curvFilt, tMeanFilt ) )
  length(naElOut)
  
  intensityFilt <- intensity[-naElOut,,]
  rankProfilesFilt <- rankProfiles[-naElOut]
  surfValFilt <- surfVal[-naElOut]
  storeCoeffsFilt <- storeCoeffs[-naElOut,]
  storeResFilt <- storeRes[-naElOut,]
  
  plot( apply( storeResFilt[,xVar], 2, median ) )
  
  plot( apply( intensityFilt[,xVar,1], 2, mean ) )
  plot( apply( intensityFilt[,xVar,2], 2, mean ) )
  
  #plot
  
  rankProfilesBin <- cut( rankProfilesFilt, breaks=quantile( rankProfilesFilt, probs=seq(0,1,0.33) ), include.lowest=TRUE )
  rankProfilesBinNumeric <- as.numeric( rankProfilesBin )
  
  x11(height=2.5, width=7)
  par(mfrow=c(1,3))
  aggData <- aggregate( intensityFilt[,xVar,1], by=list(rankProfilesBinNumeric), FUN=median)
  matplot( t( aggData[ , 2:dim(aggData)[2] ] ) )
  aggData <- aggregate( intensityFilt[,xVar,2], by=list(rankProfilesBinNumeric), FUN=median)
  matplot( t( aggData[ , 2:dim(aggData)[2] ] ) )
  aggData <- aggregate( intensityFilt[,xVar,5], by=list(rankProfilesBinNumeric), FUN=median)
  matplot( t( aggData[ , 2:dim(aggData)[2] ] ) )
  
  perThr <- 0.45
  minRank <- 0
  rankThr <- quantile( rankProfilesFilt, probs=c(minRank, perThr) )
  splitVar <- ifelse( rankProfilesFilt<rankThr[2], 1, 2 )
  splitVar[ rankProfilesFilt<rankThr[1] ] <- 0
  table(splitVar)
  
  x11(height=5, width=10)
  xValues <- seq(0,1,length.out=length(xVar))
  par(mfrow=c(2,4))
  plot( apply( intensityFilt[ splitVar==1, xVar, 1 ], 2, mean )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5 )
  plot( apply( intensityFilt[ splitVar==1, xVar, 2 ], 2, mean )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5  )
  plot( apply( intensityFilt[ splitVar==1, xVar, 5 ], 2, mean )~xValues, bty='n', las=1, ylab='raw BOLD', xlab='cortical depth', pch=17, cex=1.5  )
  plot( apply( intensityFilt[ splitVar==1, xVar, 7 ], 2, mean )~xValues, bty='n', las=1, ylab='T1 intensity', xlab='cortical depth', pch=17, cex=1.5  )
  
  plot( apply( intensityFilt[ splitVar==2, xVar, 1 ], 2, mean )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5   )
  plot( apply( intensityFilt[ splitVar==2, xVar, 2 ], 2, mean )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5   )
  plot( apply( intensityFilt[ splitVar==2, xVar, 5 ], 2, mean )~xValues, bty='n', las=1, ylab='raw BOLD', xlab='cortical depth', pch=17, cex=1.5   )
  plot( apply( intensityFilt[ splitVar==2, xVar, 7 ], 2, mean )~xValues, bty='n', las=1, ylab='T1 intensity', xlab='cortical depth', pch=17, cex=1.5   )
  
  pSteps <- seq(0.3,0.95,0.01)
  posProfileMatrix <- array( 0, c( length( pSteps ), 25 ) )
  for ( k in 1:length( pSteps ) ) {
    rankThr <- quantile( rankProfilesFilt, probs=pSteps[k] )
    splitVar <- ifelse( rankProfilesFilt<=rankThr, 1, 2 )
    
    posProfileTemp <- apply( intensityFilt[ splitVar==1, xVar, 1 ], 2, mean )
    sp <- smooth.spline( posProfileTemp )
    xOut <- seq(1,14,length.out=25)
    spInt <- predict(sp,xOut)
    
    posProfileMatrix[k,] <- scaleData( spInt$y, 1, 0 )  
  }
  x11(height=4, width=3.5)
  image( t( posProfileMatrix ), col=rainbow(1000, start=0.2, end=0.9), las=1 )
  
  
  
  ## write out ROI as a map (to double check the mapping)
  nBins <- length( unique( rankProfilesBinNumeric ) )
  rankProfilesBinNumericFlag <- array( 0, c( length( rankProfilesBinNumeric ), nBins ) )
  for ( k in 1:nBins ) {
    rankProfilesBinNumericFlag[ rankProfilesBinNumeric==k, k ] <- 1
  }
  emptyMap <- array( 0, c( dim( mapValListStat[[10]] )[1], 3+nBins ) )
  emptyMap[,1] <- seq( 0, dim(emptyMap)[1]-1 )
  emptyMap[ surfValFilt+1, 2 ] <- 1
  emptyMap[ surfValFilt+1, 3 ] <- rankProfilesBinNumeric
  emptyMap[ surfValFilt+1, seq(4,3+nBins) ] <- rankProfilesBinNumericFlag
  setwd(mainDir)
  write.table( emptyMap, file='storeMap01.1D', row.names=FALSE, col.names=FALSE )
  
  
  #tapply( intensityFilt[,xVar,1], list( rankProfilesBinNumeric ), median  )
  #apply( intensityFilt, 2, median )
  
  # x11(width=6, height=4)
  # par(mfrow=c(2,3))
  # for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(1,12), bty='n', col=k, type='l' ) }
  #   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), col=k ) }  
  # }
  # for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(-5,0.5), bty='n', col=k, type='l' ) }
  #   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), col=k ) }  
  # }
  # for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 5 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(50000,110000), bty='n', col=k, type='l' ) }
  #   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 5 ], 2, median ), col=k ) }  
  # }
  # for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #   if (k==1) { plot( xVar, apply( abs( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ] ), 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(0,5), bty='n', col=k, type='l' ) }
  #   if (k>1)  { lines( xVar, apply( abs( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ] ), 2, median ), col=k ) }  
  # }
  # for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 7 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(15000,40000), bty='n', col=k, type='l' ) }
  #   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 7 ], 2, median ), col=k ) }  
  # }
  
  
  
  ## shift data accordingly
  angleShift <- 45
  angles <- intensityFilt[,,9] + angleShift
  angles <- ifelse( angles>180, angles-180, angles )
  angles <- ifelse( angles<0, angles+180, angles )
  
  x11(height=6, width=9)
  par(mfrow=c(3,5))
  xValues <- seq(0,1,length.out=length(xVar))
  storeValues <- array(0,c(length(xVar),1))
  for (k in 1:length(xVar)) {
    surfNum <- k
    B0angleBreaks <- cut( angles[,surfNum], breaks = quantile( angles[,surfNum], seq(0,1,0.1) ), include.lowest=TRUE )
    #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
    bold_angle <- tapply( intensityFilt[,surfNum,1], list( B0angleBreaks ), median  )
    B0_angle <- tapply( angles[,surfNum], list( B0angleBreaks ), median  )
    B0_angle2 <- B0_angle^2 
    cosPred <- cos( (B0_angle)*pi/180 )^2
    linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
    modCos <- lm( bold_angle ~ cosPred ) 
    sModCos <- summary( modCos )
    predMod <- predict( modCos, data.frame( cosPred ) )
    
    #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
    #predMod01 <- predict( modCos, data.frame( x01 ) )  
    
    plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1 )    
    interpData <- spline( B0_angle, predMod, n=50, xout=seq(40,140,10) )
    
    #lines( B0_angle, predMod, lwd=2, col='red' )
    lines( interpData$x, interpData$y, lwd=2, col='red' )
    
    #lines( seq(40,170,length.out=length(cosPred)), predMod01, lwd=2, col='blue' )
    storeValues[k] <- sModCos$r.squared 
  }
  
  
  angles <- intensityFilt[,,9] + angleShift
  angles <- ifelse( angles>180, angles-180, angles )
  angles <- ifelse( angles<0, angles+180, angles )
  
  x11(height=6, width=9)
  par(mfrow=c(3,5))
  xValues <- seq(0,1,length.out=length(xVar))
  storeValues <- array(0,c(length(xVar),1))
  for (k in 1:length(xVar)) {
    surfNum <- k
    B0angleBreaks <- cut( angles[,surfNum], breaks = quantile( angles[,surfNum], seq(0,1,0.1) ), include.lowest=TRUE )
    #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
    bold_angle <- tapply( intensityFilt[,surfNum,2], list( B0angleBreaks ), median  )
    B0_angle <- tapply( angles[,surfNum], list( B0angleBreaks ), median  )
    B0_angle2 <- B0_angle^2 
    cosPred <- cos( (B0_angle)*pi/180 )^2
    linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
    modCos <- lm( bold_angle ~ cosPred ) 
    sModCos <- summary( modCos )
    predMod <- predict( modCos, data.frame( cosPred ) )
    
    #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
    #predMod01 <- predict( modCos, data.frame( x01 ) )  
    
    plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1 )    
    interpData <- spline( B0_angle, predMod, n=50, xout=seq(40,140,10) )
    
    #lines( B0_angle, predMod, lwd=2, col='red' )
    lines( interpData$x, interpData$y, lwd=2, col='red' )
    
    #lines( seq(40,170,length.out=length(cosPred)), predMod01, lwd=2, col='blue' )
    storeValues[k] <- sModCos$r.squared 
  }
  
  angles <- intensityFilt[,,9] + angleShift
  angles <- ifelse( angles>180, angles-180, angles )
  angles <- ifelse( angles<0, angles+180, angles )
  
  x11(height=6, width=9)
  par(mfrow=c(3,5))
  xValues <- seq(0,1,length.out=length(xVar))
  storeValues <- array(0,c(length(xVar),1))
  for (k in 1:length(xVar)) {
    surfNum <- k
    B0angleBreaks <- cut( angles[,surfNum], breaks = quantile( angles[,surfNum], seq(0,1,0.1) ), include.lowest=TRUE )
    #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
    bold_angle <- tapply( intensityFilt[,surfNum,5], list( B0angleBreaks ), median  )
    B0_angle <- tapply( angles[,surfNum], list( B0angleBreaks ), median  )
    B0_angle2 <- B0_angle^2 
    cosPred <- cos( (B0_angle)*pi/180 )^2
    linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
    modCos <- lm( bold_angle ~ cosPred ) 
    sModCos <- summary( modCos )
    predMod <- predict( modCos, data.frame( cosPred ) )
    
    #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
    #predMod01 <- predict( modCos, data.frame( x01 ) )  
    
    plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1 )    
    interpData <- spline( B0_angle, predMod, n=50, xout=seq(40,140,10) )
    
    #lines( B0_angle, predMod, lwd=2, col='red' )
    lines( interpData$x, interpData$y, lwd=2, col='red' )
    
    #lines( seq(40,170,length.out=length(cosPred)), predMod01, lwd=2, col='blue' )
    storeValues[k] <- sModCos$r.squared 
  }
  
  setwd(mainDir)
  fname <-sprintf('%s_image_V1_no_smoothing_allData.RData', participant)
  save.image(file=fname)
  
}



#dev.copy2pdf(file='t1_intensity.pdf', height=2.5, width=9)
#dev.off()





# setwd(mainDir)
# B0map <- array( 0 , c( dim( mapValListAmp[[10]] )[1], 1 ) )
# B0mapApp <- mapValListB0Angle[[10]]
# B0map <- B0mapApp[,2] + angleShift
# B0map <- ifelse( B0map>180, B0map-180, B0map )
# B0map <- ifelse( B0map<0, B0map+180, B0map )
# write.table( B0map, file='B0mapCorr.1D.dset', row.names=FALSE, col.names=FALSE )
# 
# 
# summary(B0map)
# 
# 
# 
# medianB0Angle <- apply( angles, 1, median )
# rankProfilesBin <- cut( medianB0Angle, breaks=seq(0,180,45), include.lowest=TRUE )
# rankProfilesBinNumeric <- as.numeric( rankProfilesBin )
# 
# medianBold <- apply( intensityFilt[,,1], 1, median  )
# tapply( medianBold, list(rankProfilesBinNumeric), median )
# 
# x11(width=6, height=4)
# par(mfrow=c(2,3))
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
#   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(1,8), bty='n', col=k, type='l' ) }
#   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
#   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(-5,0.5), bty='n', col=k, type='l' ) }
#   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
#   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 5 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(200000,380000), bty='n', col=k, type='l' ) }
#   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 5 ], 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
#   if (k==1) { plot( xVar, apply( abs( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ] ), 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(0,5), bty='n', col=k, type='l' ) }
#   if (k>1)  { lines( xVar, apply( abs( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ] ), 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
#   if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 7 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(15000,40000), bty='n', col=k, type='l' ) }
#   if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 7 ], 2, median ), col=k ) }  
# }
# 
# 
# 
# #plot( intensityFilt[,9,1]~intensityFilt[,9,9] )
# 
# #boldArray <- array( t( intensityFilt[,xVar,1] ) )
# #curvArray <- rep( apply( intensityFilt[,xVar,3], 1, mean ), rep( length(xVar), dim(intensityFilt)[1] )  )
# #thickArray <- rep( apply( intensityFilt[,xVar,4], 1, mean ), rep( length(xVar), dim(intensityFilt)[1] )  )
# #ampArray <- rep( apply( intensityFilt[,xVar,5], 1, mean ), rep( length(xVar), dim(intensityFilt)[1] )  )
# #depthArray <- rep( scaleData( xVar, 1, 0 ), rep( dim(intensityFilt)[1] )  )
# #profileId <- as.factor( rep( seq(1,dim(intensityFilt)[1]),  rep( length(xVar), dim(intensityFilt)[1] ) ) )
# 
# #library(lme4)
# 
# #mod01 <- lmer( boldArray ~ curvArray + thickArray + depthArray  + (1|profileId) )
# #mod02 <- lmer( boldArray ~ curvArray + thickArray + depthArray + depthArray*thickArray  + (1|profileId) )
# 
# #summary(mod01)
# #summary(mod02)
# 
# #anova( mod01, mod02 )
# 
# #pMod02 <- predict( mod02 )
# #qPMod <- quantile( pMod02, seq(0.2,1,0.02) )
# #qBold <- quantile( boldArray, seq(0.2,1,0.02) )
# 
# 
# 
# #library(quantreg)
# 
# meanProfileBold <- apply( intensityFilt[,xVar,8], 1, median )
# meanProfileAmp <- apply( intensityFilt[,xVar,5], 1, median )
# meanProfileThick <- apply( intensityFilt[,xVar,4], 1, median )
# meanProfileCurv <- apply( intensityFilt[,xVar,3], 1, median )
# meanProfileT <- apply( intensityFilt[,xVar,8], 1, median )
# ampPred <- as.numeric( cut( meanProfileAmp, breaks = quantile( meanProfileAmp, seq(0,1,0.5) ), include.lowest=TRUE ) )
# thickPred <- as.numeric( cut( meanProfileThick, breaks = quantile( meanProfileThick, seq(0,1,0.25) ), include.lowest=TRUE ) )
# curvPred <- as.numeric( cut( meanProfileCurv, breaks = quantile( meanProfileCurv, seq(0,1,0.5) ), include.lowest=TRUE ) )
# 
# library( quantreg )
# mod <- rq( meanProfileBold ~ as.factor( curvPred ) * as.factor( thickPred ) )
# summary(mod)
# pmod <- predict(mod)
# 
# tMat <- tapply( meanProfileT,  list( as.numeric( curvPred ) , as.numeric( thickPred ) ), median )
# image( tMat )
# 
# plot( xVar, apply( intensityFilt[  curvPred==1 & thickPred==1, xVar , 1], 2, median ), type='l', col=gray(0), lty=1, lwd=2, ylim=c(0,6) )
# lines( xVar, apply( intensityFilt[  curvPred==1 & thickPred==2, xVar , 1], 2, median ), col=gray(0.2), lty=1, lwd=2 )
# lines( xVar, apply( intensityFilt[  curvPred==1 & thickPred==3, xVar , 1], 2, median ), col=gray(0.4), lty=1, lwd=2 )
# lines( xVar, apply( intensityFilt[  curvPred==1 & thickPred==4, xVar , 1], 2, median ), col=gray(0.8), lty=1, lwd=2 )
# 
# plot( xVar, apply( intensityFilt[  curvPred==2 & thickPred==1, xVar , 1], 2, median ), type='l', col=gray(0), lty=1, lwd=2, ylim=c(0,6) )
# lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==2, xVar , 1], 2, median ), col=gray(0.2), lty=1, lwd=2 )
# lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==3, xVar , 1], 2, median ), col=gray(0.4), lty=1, lwd=2 )
# lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==4, xVar , 1], 2, median ), col=gray(0.8), lty=1, lwd=2 )
# 
# plot( xVar, apply( intensityFilt[  thickPred==1, xVar , 1], 2, median ), type='l', col=gray(0), lty=1, lwd=2, ylim=c(0,6) )
# lines( xVar, apply( intensityFilt[  thickPred==2, xVar , 1], 2, median ), col=gray(0.2), lty=1, lwd=2 )
# lines( xVar, apply( intensityFilt[  thickPred==3, xVar , 1], 2, median ), col=gray(0.4), lty=1, lwd=2 )
# lines( xVar, apply( intensityFilt[  thickPred==4, xVar , 1], 2, median ), col=gray(0.8), lty=1, lwd=2 )
# 
# 
# 
# #plot( xVar, apply( intensityFilt[  curvPred==1 & thickPred==2 & ampPred==1, xVar , 1], 2, median ), type='l', col=gray(0.9), lty=1, lwd=2, ylim=c(1,9) )
# #lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==2 & ampPred==1, xVar , 1], 2, median ), col=gray(0.9), lty=2, lwd=2 )
# #lines( xVar, apply( intensityFilt[  curvPred==1 & thickPred==2 & ampPred==2, xVar , 1], 2, median ), col=gray(0.4), lty=1, lwd=2 )
# #lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==2 & ampPred==2, xVar , 1], 2, median ), col=gray(0.4), lty=2, lwd=2 )
# 
# 
# #plot( xVar, apply( intensityFilt[  curvPred==1 & thickPred==1, xVar , 8], 2, median ), type='l', col=gray(0.6), lty=1, lwd=2, ylim=c(0.5,2.5) )
# #lines( xVar, apply( intensityFilt[  curvPred==1 & thickPred==2, xVar , 8], 2, median ), col=gray(0), lty=1, lwd=2 )
# #lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==1, xVar , 8], 2, median ), col=gray(0.6), lty=2, lwd=2 )
# #lines( xVar, apply( intensityFilt[  curvPred==2 & thickPred==2, xVar , 8], 2, median ), col=gray(0), lty=2, lwd=2 )
# 
# 
# #plot( tapply( apply( intensityFilt[ , xVar, 1 ], 1, median ), list( ampPred ), median ) )
# 
# 
# 
# #plot( xVar, apply( intensityFilt[  ampPred==1 & thickPred==1 & curvPred==1, xVar , 1], 2, mean ), type='l', col=1, ylim=c(0,5) )
# #lines( xVar, apply( intensityFilt[  ampPred==1 & thickPred==2 & curvPred==1, xVar , 1], 2, mean ), col=2 )
# #lines( xVar, apply( intensityFilt[  ampPred==2 & thickPred==1 & curvPred==1, xVar , 1], 2, mean ), col=1, lty=2 )
# #lines( xVar, apply( intensityFilt[  ampPred==2 & thickPred==2 & curvPred==1, xVar , 1], 2, mean ), col=2, lty=2 )
# 
# #plot( xVar, apply( intensityFilt[  ampPred==1 & thickPred==1 & curvPred==2, xVar , 1], 2, mean ), type='l', col=1, ylim=c(0,5) )
# #lines( xVar, apply( intensityFilt[  ampPred==1 & thickPred==2 & curvPred==2, xVar , 1], 2, mean ), col=2 )
# #lines( xVar, apply( intensityFilt[  ampPred==2 & thickPred==1 & curvPred==2, xVar , 1], 2, mean ), col=1, lty=2 )
# #lines( xVar, apply( intensityFilt[  ampPred==2 & thickPred==2 & curvPred==2, xVar , 1], 2, mean ), col=2, lty=2 )
# 
# # scaledInt <- scale( t( intensityFilt[,xVar,1] ) )
# # mydata <- t( as.matrix( scaledInt ) )
# # sum( is.na( mydata ) )
# # naflag <- apply( mydata, 1, is.na )
# # sumnaflag <- apply( naflag, 2, sum )
# # sumnaflagProfile <- sumnaflag == 0
# # mydata <- mydata[ sumnaflagProfile, ]
# # 
# # library(cluster)
# # clust.out <- pam(mydata, 3)
# # prop.table( table( clust.out$clustering ) )
# # aggData <- aggregate(intensityFilt[,xVar,1],by=list(clust.out$clustering),FUN=median)
# # aggDataNeg <- aggregate(intensityFilt[,xVar,2],by=list(clust.out$clustering),FUN=median)
# # 
# # # plot
# # x11(height=4.5, width=4)
# # xValues <- seq(0,1,length.out=length(xVar))
# # par(mfrow=c(2,2))
# # plot( as.numeric( aggData[ 1, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# # plot( as.numeric( aggData[ 2, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# # plot( as.numeric( aggDataNeg[ 1, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# # plot( as.numeric( aggDataNeg[ 2, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# # 
# # 
# # distMat <- dist(  mydata, method='manhattan' )
# 
# 
# 
# 
# 
# rankProfilesBin <- cut( rankProfilesFilt, breaks=quantile( rankProfilesFilt, probs=seq(0,1,0.25) ), include.lowest=TRUE )
# rankProfilesBinNumeric <- as.numeric( rankProfilesBin )
# 
# 
# x11(width=6, height=4)
# par(mfrow=c(2,3))
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
# if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(1,12), bty='n', col=k, type='l' ) }
# if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
# if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(-5,0.5), bty='n', col=k, type='l' ) }
# if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
# if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 5 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(50000,110000), bty='n', col=k, type='l' ) }
# if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 5 ], 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
# if (k==1) { plot( xVar, apply( abs( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ] ), 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(0,5), bty='n', col=k, type='l' ) }
# if (k>1)  { lines( xVar, apply( abs( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ] ), 2, median ), col=k ) }  
# }
# for (k in 1:length(unique(rankProfilesBinNumeric))) {
# if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 7 ], 2, median ), xlim=c(min(xVar),max(xVar)), ylim=c(15000,40000), bty='n', col=k, type='l' ) }
# if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 7 ], 2, median ), col=k ) }  
# }
# 
# 
# meanCurvature <- apply( intensityFilt[,xVar,3], 1, median )
# meanThickness <- apply( intensityFilt[,xVar,4], 1, median )
# x11(height=3.5, width=5.5)
# par(mfrow=c(1,2))
# plot( rankProfilesFilt~meanCurvature, cex=0.1 ) 
# plot( rankProfilesFilt~meanThickness, cex=0.1 ) 
# dev.copy2pdf(file='curvature_thickness.pdf', height=3.5, width=5.5)
# dev.off()
# 
# dev.copy2pdf(file='general_plot.pdf', height=5, width=10)
# dev.off()
# 
# 
# perThr <- 0.45
# minRank <- 0.025
# rankThr <- quantile( rankProfilesFilt, probs=c(minRank, perThr) )
# splitVar <- ifelse( rankProfilesFilt<rankThr[2], 1, 2 )
# splitVar[ rankProfilesFilt<rankThr[1] ] <- 0
# table(splitVar)
# 
# x11(height=5, width=10)
# xValues <- seq(0,1,length.out=length(xVar))
# par(mfrow=c(2,4))
# plot( apply( intensityFilt[ splitVar==1, xVar, 1 ], 2, median )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5 )
# plot( apply( intensityFilt[ splitVar==1, xVar, 2 ], 2, median )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5  )
# plot( apply( intensityFilt[ splitVar==1, xVar, 5 ], 2, median )~xValues, bty='n', las=1, ylab='raw BOLD', xlab='cortical depth', pch=17, cex=1.5  )
# plot( apply( intensityFilt[ splitVar==1, xVar, 7 ], 2, median )~xValues, bty='n', las=1, ylab='T1 intensity', xlab='cortical depth', pch=17, cex=1.5  )
# 
# plot( apply( intensityFilt[ splitVar==2, xVar, 1 ], 2, median )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5   )
# plot( apply( intensityFilt[ splitVar==2, xVar, 2 ], 2, median )~xValues, bty='n', las=1, ylab='%BOLD', xlab='cortical depth', pch=17, cex=1.5   )
# plot( apply( intensityFilt[ splitVar==2, xVar, 5 ], 2, median )~xValues, bty='n', las=1, ylab='raw BOLD', xlab='cortical depth', pch=17, cex=1.5   )
# plot( apply( intensityFilt[ splitVar==2, xVar, 7 ], 2, median )~xValues, bty='n', las=1, ylab='T1 intensity', xlab='cortical depth', pch=17, cex=1.5   )
# 
# pSteps <- seq(0.3,0.95,0.01)
# posProfileMatrix <- array( 0, c( length( pSteps ), 25 ) )
# for ( k in 1:length( pSteps ) ) {
# rankThr <- quantile( rankProfilesFilt, probs=pSteps[k] )
# splitVar <- ifelse( rankProfilesFilt<=rankThr, 1, 2 )
# 
# posProfileTemp <- apply( intensityFilt[ splitVar==1, xVar, 1 ], 2, mean )
# sp <- smooth.spline( posProfileTemp )
# xOut <- seq(1,14,length.out=25)
# spInt <- predict(sp,xOut)
# 
# posProfileMatrix[k,] <- scaleData( spInt$y, 1, 0 )  
# }
# x11(height=4, width=3.5)
# image( t( posProfileMatrix ), col=rainbow(1000, start=0.2, end=0.9), las=1 )
# 
# dev.copy2pdf(file='heat_map.pdf', height=4, width=3.5)
# dev.off()
# 
# 
# ## write out ROI as a map (to double check the mapping)
# nBins <- length( unique( rankProfilesBinNumeric ) )
# rankProfilesBinNumericFlag <- array( 0, c( length( rankProfilesBinNumeric ), nBins ) )
# for ( k in 1:nBins ) {
# rankProfilesBinNumericFlag[ rankProfilesBinNumeric==k, k ] <- 1
# }
# emptyMap <- array( 0, c( dim( mapValListStat[[10]] )[1], 3+nBins ) )
# emptyMap[,1] <- seq( 0, dim(emptyMap)[1]-1 )
# emptyMap[ surfValFilt+1, 2 ] <- 1
# emptyMap[ surfValFilt+1, 3 ] <- rankProfilesBinNumeric
# emptyMap[ surfValFilt+1, seq(4,3+nBins) ] <- rankProfilesBinNumericFlag
# setwd(mainDir)
# write.table( emptyMap, file='storeMap01.1D', row.names=FALSE, col.names=FALSE )
# 
# 
# scaledInt <- scale( t( intensityFilt[,xVar,1] ) )
# mydata <- t( as.matrix( scaledInt ) )
# sum( is.na( mydata ) )
# naflag <- apply( mydata, 1, is.na )
# sumnaflag <- apply( naflag, 2, sum )
# sumnaflagProfile <- sumnaflag == 0
# mydata <- mydata[ sumnaflagProfile, ]
# 
# 
# #mydatacopy <- intensityFilt[,xVar,1]
# #mydata <- apply( mydatacopy, 1, scaleData, newMax=1, newMin=0  )
# #sum( is.na( mydata ) )
# mydataCurv <- intensityFilt[,xVar,3]
# meanCurv <- apply( mydataCurv, 1, mean  )
# plot( apply( mydataCurv[meanCurv>0,], 2, mean  ) )
# plot( apply( mydataCurv[meanCurv<0,], 2, mean  ) )
# 
# 
# mydata <- intensityFilt[,xVar,7]
# #mydataPos <- intensityFilt[,xVar,1]
# #mydataNeg <- intensityFilt[,xVar,2]
# 
# # Determine number of clusters
# #wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# #for (i in 2:15) wss[i] <- sum(kmeans(mydata,
# #                                     centers=i)$withinss)
# #plot(1:15, wss, type="b", xlab="Number of Clusters",
# #     ylab="Within groups sum of squares")
# 
# # K-Means Cluster Analysis
# #fit <- kmeans(mydata, 5, iter.max=50 ) # 5 cluster solution
# #prop.table( table( fit$cluster ) )
# # get cluster means
# #aggData <- aggregate(mydata,by=list(fit$cluster),FUN=mean)
# #aggDataPos <- aggregate(mydataPos,by=list(fit$cluster),FUN=mean)
# #aggDataNeg <- aggregate(mydataNeg,by=list(fit$cluster),FUN=mean)
# 
# library(cluster)
# library(randomForest)
# clust.out <- pam(mydata, 4)
# prop.table( table( clust.out$clustering ) )
# aggData <- aggregate(mydata,by=list(clust.out$clustering),FUN=median)
# #plot( clust.out )
# out <- prcomp(mydata, scale=FALSE)
# 
# # plot
# x11(height=2.5, width=9)
# xValues <- seq(0,1,length.out=length(xVar))
# par(mfrow=c(1,4))
# plot( as.numeric( aggData[ 1, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# plot( as.numeric( aggData[ 2, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# plot( as.numeric( aggData[ 3, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# plot( as.numeric( aggData[ 4, 2:dim(aggData)[2] ] )~xValues, bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=17, cex=1.5 )
# 
# dev.copy2pdf(file='t1_intensity.pdf', height=2.5, width=9)
# dev.off()
# 
# 
# mydataShort <- mydata[1:300,]
# iris.urf <- randomForest( mydataShort, ntree=5000 )
# 
# iris.mds <- cmdscale(1 - iris.urf$proximity, eig=TRUE)
# plot( iris.mds$points[,1],iris.mds$points[,2] )
# 
# 
# 
# 
# library(randomForest)
# dim( t( mydata ) )
# myDataSel <- mydataPos[ sample( seq( 1, dim(mydata)[1] ), size=500  ), ]
# p.urf <- randomForest( myDataSel, ntree=5000 )
# p.scale <- cmdscale(1 - p.urf$proximity, eig=TRUE)
# plot( p.scale$points )
# 
# 
# peakInfo <- array( 0, c( dim(intensityValuesMat)[1], 2 ) )
# for (k in 1:dim(intensityValuesMat)[1] ) {
# intArray <- as.numeric( scale( intensityValuesMat[k,1:10] ) )
# a <- findpeaks( intArray, minpeakheight=0.2 )
# if ( is.numeric(a) ) {
# peakLocation <- a[,2]
# flagPeak <- peakLocation > 2 & peakLocation <= 7
# whichPeakRow <- which(flagPeak) ## separate this from integer 0 !!
# if is.integer( whichPeakRow[1] ) { peakInfo[k,1] <- 1 }
# }  
# }
# sum( peakInfo[,1]==0 )/length(peakInfo[,1])
# 
# 
# ## write out ROI as a map (to double check the mapping)
# emptyMap <- array( 0, c( dim( mapValList[[5]] )[1], 3 ) )
# surfVal <- as.data.frame( surfValList[[5]] )
# emptyMap[ surfVal[,1]+1, 2 ] <- 1
# emptyMap[ surfVal[,1]+1, 3 ] <- peakInfo[,1]
# emptyMap[,1] <- seq( 0, dim(emptyMap)[1]-1 )
# write.table( emptyMap, file='storeRoi03.1D.dset', row.names=FALSE, col.names=FALSE )
# 
# 
# 
# 
# 
# 
# 
# mydata <- intensity[,,1]
# # Determine number of clusters
# wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
# for (i in 2:15) wss[i] <- sum(kmeans(mydata,
# centers=i)$withinss)
# plot(1:15, wss, type="b", xlab="Number of Clusters",
# ylab="Within groups sum of squares")
# 
# # K-Means Cluster Analysis
# fit <- kmeans(mydata, 5) # 5 cluster solution
# prop.table( table( fit$cluster ) )
# # get cluster means
# aggData <- aggregate(mydata,by=list(fit$cluster),FUN=mean)
# # append cluster assignment
# plot( as.numeric( aggData[5,2:10] ) )
# 
# 
# # Ward Hierarchical Clustering
# d <- dist(mydata, method = "manhattan") # distance matrix
# fit <- hclust(d, method="ward")
# plot(fit) # display dendogram
# groups <- cutree(fit, k=5) # cut tree into 5 clusters
# # draw dendogram with red borders around the 5 clusters
# rect.hclust(fit, k=5, border="red") 
# 
# 
# 
# 
# 
# 
# ## load ROI files
# setwd( sprintf('%s/%s', mainDir, roiDir ) )
# fileRoiSurface <- dir( pattern='^roiSurface.*.1D.roi')
# 
# surfValList <- list(0)
# fileNumber <- c(1,2,3,4,5,6,8,9,10,11,12) #depending on which surface you use as baseline, this line has to change!!
# for ( k in 1:length(fileNumber) ) {
# fileIndex <- fileNumber[k]
# surfValList[[k]] <- read.table( fileRoiSurface[ fileIndex ], as.is=TRUE )
# }
# dfFiles <- data.frame(fileRoiSurface)
# data.frame( dfFiles[fileNumber,] )
# 
# 
# ## load stat files
# setwd( sprintf('%s/%s', mainDir, statDir ) )
# fileSurfacesMapping <- dir( pattern='*surfval.1D')
# 
# mapValListStat <- list(0)
# fileNumber <- c(1,2,3,4,5,6,7,8,9,10,11)
# for ( k in 1:length(fileNumber) ) {
# fileIndex <- fileNumber[k]
# mapValListStat[[k]] <- read.table( fileSurfacesMapping[ fileIndex ], skip=5, as.is=TRUE )  
# }
# data.frame(fileSurfacesMapping)
# 
# ## load curvature files
# setwd( sprintf('%s/%s', mainDir, curvatureDir ) )
# fileSurfacesMapping <- dir( pattern='*surfval.1D')
# 
# mapValListCurv <- list(0)
# fileNumber <- c(1,2,3,4,5,6,7,8,9,10,11)
# for ( k in 1:length(fileNumber) ) {
# fileIndex <- fileNumber[k]
# mapValListCurv[[k]] <- read.table( fileSurfacesMapping[ fileIndex ], skip=5, as.is=TRUE )  
# }
# data.frame(fileSurfacesMapping)
# 
# ## load thickness files
# setwd( sprintf('%s/%s', mainDir, thicknessDir ) )
# fileSurfacesMapping <- dir( pattern='*surfval.1D')
# 
# mapValListThick <- list(0)
# fileNumber <- c(1,2,3,4,5,6,7,8,9,10,11)
# for ( k in 1:length(fileNumber) ) {
# fileIndex <- fileNumber[k]
# mapValListThick[[k]] <- read.table( fileSurfacesMapping[ fileIndex ], skip=5, as.is=TRUE )  
# }
# data.frame(fileSurfacesMapping)
# 
# ## load amplitude files
# setwd( sprintf('%s/%s', mainDir, amplitudeDir ) )
# fileSurfacesMapping <- dir( pattern='*surfval.1D')
# 
# mapValListAmp <- list(0)
# fileNumber <- c(1,2,3,4,5,6,7,8,9,10,11)
# for ( k in 1:length(fileNumber) ) {
# fileIndex <- fileNumber[k]
# mapValListAmp[[k]] <- read.table( fileSurfacesMapping[ fileIndex ], skip=5, as.is=TRUE )  
# }
# data.frame(fileSurfacesMapping)
# 
# 
# ## generate profiles from intensity files and roi files
# intensityValuesMat <- array(0, c( dim(surfValList[[1]])[1], 11, 5 ) )
# for (k in 1:11) {
# map <- as.data.frame( mapValListStat[[k]] )
# curv <- as.data.frame( mapValListCurv[[k]] )
# thick <- as.data.frame( mapValListThick[[k]] )
# amp <- as.data.frame( mapValListAmp[[k]] )
# surfVal <- as.data.frame( surfValList[[k]] )
# 
# unconnectedPoints <- which( surfVal[,1]==-1 | surfVal[,1]==1 )
# print( unconnectedPoints )
# surfValClean <- ifelse( surfVal[,1]==-1, 1, surfVal[,1] )
# 
# if ( length( unconnectedPoints )>0 && exists('storeUnconnectedPoints') ) {
# storeUnconnectedPoints <- c( storeUnconnectedPoints, unconnectedPoints )
# }    
# if ( length( unconnectedPoints )>0 && !exists('storeUnconnectedPoints') ) {
# storeUnconnectedPoints <- unconnectedPoints
# }
# 
# intValue <- map[ surfValClean, c(8,11) ]  
# curvValue <- curv[ surfValClean, 7 ]
# thickValue <- thick[ surfValClean, 7 ]
# ampValue <- amp[ surfValClean, 7 ]
# 
# intensityValuesMat[,k,1] <- intValue[,1]
# intensityValuesMat[,k,2] <- intValue[,2]
# intensityValuesMat[,k,3] <- curvValue
# intensityValuesMat[,k,4] <- thickValue
# intensityValuesMat[,k,5] <- ampValue
# }
# intensityValuesMat <- intensityValuesMat[ -1*storeUnconnectedPoints, , ]
# 
# 
# load1DData <- function(directory, fileNumber, skipLines, filePattern) {  
# setwd( directory )
# fileSurfacesMapping <- dir( pattern=filePattern)
# mapValList <- list(0)  
# print( sprintf('%d files to load', length(fileNumber) ) )
# 
# for ( k in 1:length(fileNumber) ) {
# fileIndex <- fileNumber[k]
# print( sprintf('reading file n. %d: %s ...', k, fileSurfacesMapping[ fileIndex ] ) )
# mapValList[[k]] <- read.table( fileSurfacesMapping[ fileIndex ], skip=skipLines, as.is=TRUE )        
# }
# print( 'done' )
# return( mapValList )  
# }
# 
# 
# generateProfiles01 <- function( surfValList, mapValListStat, mapValListCurv, mapValListThick, mapValListAmp, statIndexPositive, statIndexNegative ) {  
# ## generate profiles from intensity files and roi files
# intensityValuesMat <- array(0, c( dim(surfValList[[1]])[1], 11, 5 ) )
# 
# for (k in 1:11) {
# 
# print( sprintf('reading data from boundary %d', k ) )
# 
# map <- as.data.frame( mapValListStat[[k]] )
# curv <- as.data.frame( mapValListCurv[[k]] )
# thick <- as.data.frame( mapValListThick[[k]] )
# amp <- as.data.frame( mapValListAmp[[k]] )
# surfVal <- as.data.frame( surfValList[[k]] )
# 
# unconnectedPoints <- which( surfVal[,1]==-1 | surfVal[,1]==1 )
# print( unconnectedPoints )
# surfValClean <- ifelse( surfVal[,1]==-1, 1, surfVal[,1] )
# 
# if ( length( unconnectedPoints )>0 && exists('storeUnconnectedPoints') ) {
# storeUnconnectedPoints <- c( storeUnconnectedPoints, unconnectedPoints )
# }    
# if ( length( unconnectedPoints )>0 && !exists('storeUnconnectedPoints') ) {
# storeUnconnectedPoints <- unconnectedPoints
# }
# 
# intValue <- map[ surfValClean, c(statIndexPositive,statIndexNegative) ]  
# curvValue <- curv[ surfValClean, 7 ]
# thickValue <- thick[ surfValClean, 7 ]
# ampValue <- amp[ surfValClean, 7 ]
# 
# intensityValuesMat[,k,1] <- intValue[,1]
# intensityValuesMat[,k,2] <- intValue[,2]
# intensityValuesMat[,k,3] <- curvValue
# intensityValuesMat[,k,4] <- thickValue
# intensityValuesMat[,k,5] <- ampValue
# }
# 
# 
# if (exists('storeUnconnectedPoints')) { 
# print( storeUnconnectedPoints ) 
# #intensityValuesMat <- intensityValuesMat[ -1*storeUnconnectedPoints, , ]
# out <- list( intensity=intensityValuesMat,
# unconnectedPoints=storeUnconnectedPoints )  
# }
# if (!exists('storeUnconnectedPoints')) { 
# out <- list( intensity=intensityValuesMat,
# unconnectedPoints=0 )
# }
# 
# return( out )
# 
# }
