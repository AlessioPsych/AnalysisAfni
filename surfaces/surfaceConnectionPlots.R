rm( list=ls() )
library(abind)

savePlots <- 0
participantArray <- c('SD', 'BH', 'VK', 'AF', 'AF01', 'RB' ) #SD BH VK AF AF01 RB
#participantArray <- c('SD', 'BH' ) #SD BH VK AF AF01 RB
boldAngleArray <- array( 0, c( length(participantArray), 10, 16, 7 ) )
depthArrayAngle <- array( 0, c( length(participantArray), 16, 3, 4  ) )
depthArray <- array( 0, c( length(participantArray), 16, 3  ) )
storeBold <- array( 0, c( length(participantArray), 16, 3, 3  ) )
storeBoldWhole <- array( 0, c( length(participantArray), 16, 3  ) )
storeT1 <- array( 0, c( length(participantArray), 15  ) )
storeBoldEst <-  array( 0, c( length(participantArray), 3  ) )
storeT1Est <-  array( 0, c( length(participantArray), 3  ) )
boldEx <- rep(0, length(participantArray))
overallEx <- rep(0, length(participantArray))
storeMotion <- rep(0, length(participantArray))
storeNprofiles <- rep(0, length(participantArray))
storeAngleCor <-  rep(0, length(participantArray))
mainDirSavePlots <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/plots_all'


for ( nPart in 1:length(participantArray) ) {
  #for ( nPart in 1:1 ) {
  
  roiFlag <- 3
  
  #nPart <- 5
  participant <- participantArray[ nPart ]    
  
  graphics.off()
  if (participant=='BH') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V5029leftright_Ben_ph/results_CBS'
    statDir <- 'statsSingleShot+orig_interp_folder_smoothed_5'
    roiDirRight <- 'rightRoi_folder'
    roiDirLeft <- 'leftRoi_folder'
    curvatureDir <- 'curvature_interp_folder_smoothed_8'
    thicknessDir <- 'thickness_interp_folder'
    amplitudeDir <- 'amplitudeSingleShot+orig_interp_folder_smoothed_5'
    anatomyDir <- 'anatomy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_1-2-1'  
    setwd(mainDir)
    if ( roiFlag==0 ) { load( 'BH_image_roi.RData' ) }
    if ( roiFlag==1 ) { load( 'BH_image_V1.RData' ) }
    if ( roiFlag==2 ) { load( 'BH_image_V1_no_smoothing.RData' ) }
    if ( roiFlag==3 ) { load( 'BH_image_V1_no_smoothing_allData.RData' ) }
    
    angleShift <- 20
  }
  
  if (participant=='VK') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2699leftright_Vin/results_CBS'
    statDir <- 'statsSingleShotEPI_add_interp_folder_smoothed_5'
    roiDirRight <- 'rightRoi_folder'
    roiDirLeft <- 'leftRoi_folder'
    curvatureDir <- 'curvatureCopy_interp_folder_smoothed_8'
    thicknessDir <- 'thicknessCopy_interp_folder'
    amplitudeDir <- 'amplitudeAnatomyEPI_add_interp_folder_smoothed_5'
    anatomyDir <- 'anatomyCopy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_00-1'  
    setwd(mainDir)
    if ( roiFlag==0 ) { load( 'VK_image_roi.RData' )  }
    if ( roiFlag==1 ) { load( 'VK_image_V1.RData' )  }
    if ( roiFlag==2 ) { load( 'VK_image_V1_no_smoothing.RData' )  }
    if ( roiFlag==3 ) { load( 'VK_image_V1_no_smoothing_allData.RData' ) }
    angleShift <- -20
  }
  
  if (participant=='AF') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2678leftright_Ale/results_CBS'
    statDir <- 'statsSingleShotEPI_add+orig_interp_folder_smoothed_5'
    roiDirRight <- 'rightRoi_folder'
    roiDirLeft <- 'leftRoi_folder'
    curvatureDir <- 'curvature_interp_folder_smoothed_8'
    thicknessDir <- 'thickness_interp_folder'
    amplitudeDir <- 'amplitudeSingleShotEPI_add+orig_interp_folder_smoothed_5'
    anatomyDir <- 'anatomy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_-1-1-1'  
    setwd(mainDir)
    if ( roiFlag==0 ) { load( 'AF_image_roi.RData' ) }
    if ( roiFlag==1 ) { load( 'AF_image_V1.RData' ) }
    if ( roiFlag==2 ) { load( 'AF_image_V1_no_smoothing.RData' ) }
    if ( roiFlag==3 ) { load( 'AF_image_V1_no_smoothing_allData.RData' ) }
    
    angleShift <- -30
  }
  
  if (participant=='AF01') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V4847leftright_Ale_ph/results_CBS'
    statDir <- 'statsSingleShotEPI_add_interp_folder_smoothed_5'
    roiDirRight <- 'rightRoi_folder'
    roiDirLeft <- 'leftRoi_folder'
    curvatureDir <- 'curvature_interp_folder_smoothed_8'
    thicknessDir <- 'thickness_interp_folder'
    amplitudeDir <- 'amplitudeSingleShotEPI_add_interp_folder_smoothed_5'
    anatomyDir <- 'anatomy_interp_folder' 
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_222'  
    setwd(mainDir)
    if ( roiFlag==0 ) { load( 'AF01_image_roi.RData' ) }
    if ( roiFlag==1 ) { load( 'AF01_image_V1.RData' ) }
    if ( roiFlag==2 ) { load( 'AF01_image_V1_no_smoothing.RData' ) }
    if ( roiFlag==3 ) { load( 'AF01_image_V1_no_smoothing_allData.RData' ) }
    angleShift <- -30
  }
  
  
  if (participant=='SD') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2676leftright_Ser/results_CBS'  
    statDir <- 'statsSingleShotEPI_add+orig_interp_folder_smoothed_5'
    roiDirRight <- 'rightRoi_folder'
    roiDirLeft <- 'leftRoi_folder'  
    curvatureDir <- 'curvatureCopy_interp_folder_smoothed_8'
    thicknessDir <- 'thicknessCopy_interp_folder'  
    amplitudeDir <- 'amplitudeSingleShotEPI_add+orig_interp_folder_smoothed_5'
    anatomyDir <- 'anatomyCopy_interp_folder'
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_-10-1'  
    setwd(mainDir)
    if ( roiFlag==0 ) { load( 'SD_image_roi.RData' )  }
    if ( roiFlag==1 ) { load( 'SD_image_V1.RData' )  }
    if ( roiFlag==2 ) { load( 'SD_image_V1_no_smoothing.RData' )  }
    if ( roiFlag==3 ) { load( 'SD_image_V1_no_smoothing_allData.RData' ) }
    angleShift <- 0
  }
  
  if (participant=='RB') {
    mainDir <- '/media/alessiofracasso/storage2/dataLaminar_leftVSrightNEW/V2677leftright_Ric/resultsCBS'  
    statDir <- 'statsSingleShotEPI_add_interp_folder_smoothed_5'
    roiDirRight <- 'rightRoi_folder'
    roiDirLeft <- 'leftRoi_folder'  
    curvatureDir <- 'curvatureCopy_interp_folder_smoothed_8'
    thicknessDir <- 'thicknessCopy_interp_folder'
    amplitudeDir <- 'amplitudeSingleShotEPI_add_interp_folder_smoothed_5'
    anatomyDir <- 'anatomyCopy_interp_folder'
    B0angleDir <- 'B0Angles_multiple/B0_angle_smoothed_5_-2-2-2'  
    setwd(mainDir)
    if ( roiFlag==0 ) { load( 'RB_image_roi.RData' ) }
    if ( roiFlag==1 ) { load( 'RB_image_V1.RData' ) }
    if ( roiFlag==2 ) { load( 'RB_image_V1_no_smoothing.RData' ) }
    if ( roiFlag==3 ) { load( 'RB_image_V1_no_smoothing_allData.RData' ) }
    angleShift <- 60
  }
  
  
  source('/home/alessiofracasso/Dropbox/analysisAfni/load1DData.R')
  source('/home/alessiofracasso/Dropbox/analysisAfni/generateProfiles.R')
  source('/home/alessiofracasso/Dropbox/analysisAfni/scaleData.R')
  
  motionData <- read.table( file='../motionCorrect.results/dfile_rall.1D' )
  
  storeMotion[nPart] <- mean( apply( abs( as.matrix(motionData) ), c(2), max ) )
  
  
  delRows <- c( intensityRightStruct$unconnectedPoints, intensityLeftStruct$unconnectedPoints )
  filtRows <- which( delRows==0 )
  if ( length( filtRows ) > 0 ) { delRows <- delRows[-filtRows] }
  
  rm( naEl )
  if (length( delRows ) > 0) { naEl <- unique( c( which(is.na(rankProfiles)), delRows ) ) }
  if (length( delRows ) == 0) { naEl <- unique( c( which( is.na( rankProfiles ) ) ) )  }
  naEl
  
  surfValRight <- data.frame( surfValListRight[[10]] )
  surfValLeft <- data.frame( surfValListLeft[[10]] )
  surfVal <- c( surfValRight[,1], surfValLeft[,1] )
  
  dim(intensity)[1]
  length(rankProfiles)
  length(surfVal)
  
  xVar <- seq(2,17)
  
  boldMedian <- apply( intensity[,xVar,1], 1, median )  
  boldFilt <- which( boldMedian < quantile( boldMedian, c( 0.025 ) ) | boldMedian > quantile( boldMedian, c( 0.975 ) ) )
  curvMedian <- apply( intensity[ ,xVar, 3 ], 1, median  )
  curvFilt <- which( curvMedian < -0.7 | curvMedian > 0.7 )
  thickMedian <- apply( intensity[ ,xVar, 4 ], 1, median  )
  thickFilt <- which( thickMedian < 1 | thickMedian > quantile(thickMedian,0.975)   )
  anatMin <- apply( intensity[ ,xVar, 7 ], 1, min  )
  anatFilt <- which( anatMin <= 12000 | anatMin > 40000 )  
  tZero <- apply( intensity[,xVar,8], 1, median )
  tMeanFilt <- which( tZero<=2 )
  
  storeAmp <- rep( 0, dim(intensity)[1] )
  for ( k in 1:dim(intensity)[1] ) {
    sCorPearAll <- cor( intensity[ k, xVar, 5 ], xVar, method = c('pearson') )
    storeAmp[ k ] <- sCorPearAll
    #sCorPearAll <- lm( intensity[ k, xVar, 5 ] ~ xVar )
    #storeAmp[ k ] <- sCorPearAll
  }
  ampFilt <- which( storeAmp < 0 )
  
  
  
  #naElOut <- unique( c( naEl, anatFilt, thickFilt, curvFilt, tMeanFilt, ampFilt ) )
  naElOut <- unique( c( naEl, anatFilt, thickFilt, curvFilt, tMeanFilt ) )
  length(naElOut)
  
  boldEx[nPart] <- length(tMeanFilt)/dim(intensity)[1]
  overallEx[nPart] <- length(naElOut)/dim(intensity)[1]
  
  
  
  tZeroFilt <- tZero[ -naElOut ]
  intensityFilt <- intensity[-naElOut,,]
  rankProfilesFilt <- rankProfiles[-naElOut]
  surfValFilt <- surfVal[-naElOut]
  
  storeNprofiles[nPart] <-  dim(intensityFilt)[1]
  
  angles <- intensityFilt[,,9] + angleShift
  
  storeCoeff <- array( 0, c( dim(intensityFilt)[1], 6 ) )  
  storeResAngle <- array( 0, c( dim( angles )[1], length(xVar) ) )
  for ( k in 1:dim(intensityFilt)[1] ) {
    modBoldAngle <- summary( lm( intensityFilt[ k, xVar, 1 ] ~ angles[k,xVar]  ) )
    storeResAngle[k,] <- residuals( modBoldAngle ) 
    
    modBold <- summary( lm( intensityFilt[ k, xVar, 1 ] ~ xVar  ) )
    
    sCorSpear <- cor( intensityFilt[ k, xVar, 1 ], xVar, method = c('spearman') )
    sCorPear <- cor( intensityFilt[ k, xVar, 1 ], xVar, method = c('pearson') )
    
    sCorPearAmp <- cor( intensityFilt[ k, xVar, 5 ], xVar, method = c('spearman') )
    
    storeCoeff[k,] <- c( modBold$coefficients[1:2,1], modBold$r.squared, sCorSpear, sCorPear, sCorPearAmp )   
  }
  
  rankProfiles <- storeCoeff[,5]
  #rankProfiles <- storeCoeff[,2]
  rankProfilesBin <- cut( rankProfiles, breaks=quantile( rankProfiles, probs=c(0,0.40,0.70,1) ), include.lowest=TRUE )
  rankProfilesBinNumeric <- as.numeric( rankProfilesBin )
  
  #plot(  )
  
  #aggData <- aggregate( storeResAngle, by=list(rankProfilesBinNumeric), FUN=mean)  
  #matplot( t( aggData[ , 2:dim(aggData)[2] ] ), bty='n', xlab='cortical depth', ylab='% BOLD (contra-lateral)', axes=FALSE, pch=19 )
  #axis(1, seq(1,length(xVar),length.out=3), seq(0,1,0.5) )
  #axis(2, round( seq(min(aggData[ , 2:dim(aggData)[2] ]),max(aggData[ , 2:dim(aggData)[2] ]),length.out=3), 1 ), 
  #     round( seq(min(aggData[ , 2:dim(aggData)[2] ]),max(aggData[ , 2:dim(aggData)[2] ]),length.out=3), 1 ), las=1 )
  
  library( boot )
  nRepetitions <- 2000
  
  boot_fun <- function( x, indices, splitVar ) {
    mVal <- aggregate( x[indices,], by=list(splitVar[indices]), FUN=mean)
    mVal <- t( as.matrix( mVal[,2:dim(mVal)[2]] ) )
    return( mVal )
  }
#  out <- boot_fun( intensityFilt[,xVar,1], sample( dim( intensityFilt[,xVar,1] )[1], replace=TRUE ), rankProfilesBinNumeric )
#  plot( out[,1])
    
  
  x11(height=7.5, width=7)
  par(mfrow=c(3,3))
  
  bootMod <- boot( intensityFilt[,xVar,1], boot_fun, splitVar=rep(1,dim(intensityFilt[,xVar,7])[1]), R=nRepetitions )  
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )      
  aggData <- bootCi[2,]
  seMin <- bootCi[1,]
  seMax <- bootCi[3,]
  #aggData <- apply( intensityFilt[,xVar,1], 2, FUN=mean)  
  #seData <- ( apply( intensityFilt[,xVar,1], 2, FUN=sd) / sqrt( dim( intensityFilt[,xVar,1] )[1] ) ) 
  plot( aggData, bty='n', xlab='cortical depth', ylab='% BOLD (contra-lateral)', axes=FALSE, pch=19, type='l', lwd=2, lty=1, cex.lab=1.5, ylim=c( min( seMin ), max( seMax )  ) )
  xAxis <- seq(1,length(aggData))
  lines( xAxis, seMin, col = 'grey')
  lines( xAxis, seMax, col = 'grey')  
  polygon(c(xAxis, rev(xAxis)), c(seMax, rev(seMin)),
          col = rgb(0.5,0.5,0.5,0.8), border = NA)  
  #segments( seq(1,length(aggData)), seMin, seq(1,length(aggData)), seMax )  
  axis(1, seq(1,length(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
  axis(2, round( seq(min(aggData),max(aggData),length.out=3), 1 ), 
       round( seq(min(aggData),max(aggData),length.out=3), 1 ), las=1, cex.axis=1.5 )
  storeBoldWhole[ nPart, ,1 ] <- aggData
  
  bootMod <- boot( intensityFilt[,xVar,2], boot_fun, splitVar=rep(1,dim(intensityFilt[,xVar,7])[1]), R=nRepetitions )  
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )
  aggData <- bootCi[2,]
  seMin <- bootCi[1,]
  seMax <- bootCi[3,]
  #aggData <- apply( intensityFilt[,xVar,1], 2, FUN=mean)  
  #seData <- ( apply( intensityFilt[,xVar,1], 2, FUN=sd) / sqrt( dim( intensityFilt[,xVar,1] )[1] ) ) 
  plot( aggData, bty='n', xlab='cortical depth', ylab='% BOLD (contra-lateral)', axes=FALSE, pch=19, type='l', lwd=2, lty=1, cex.lab=1.5, ylim=c( min( seMin ), max( seMax )  ) )
  xAxis <- seq(1,length(aggData))
  lines( xAxis, seMin, col = 'grey')
  lines( xAxis, seMax, col = 'grey')  
  polygon(c(xAxis, rev(xAxis)), c(seMax, rev(seMin)),
          col = rgb(0.5,0.5,0.5,0.8), border = NA)  
  #segments( seq(1,length(aggData)), seMin, seq(1,length(aggData)), seMax )  
  axis(1, seq(1,length(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
  axis(2, round( seq(min(aggData),max(aggData),length.out=3), 1 ), 
       round( seq(min(aggData),max(aggData),length.out=3), 1 ), las=1, cex.axis=1.5 )
  storeBoldWhole[ nPart, ,2 ] <- aggData
  
  bootMod <- boot( intensityFilt[,xVar,5], boot_fun, splitVar=rep(1,dim(intensityFilt[,xVar,7])[1]), R=nRepetitions )  
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )      
  aggData <- bootCi[2,]
  seMin <- bootCi[1,]
  seMax <- bootCi[3,]
  #aggData <- apply( intensityFilt[,xVar,1], 2, FUN=mean)  
  #seData <- ( apply( intensityFilt[,xVar,1], 2, FUN=sd) / sqrt( dim( intensityFilt[,xVar,1] )[1] ) ) 
  plot( aggData, bty='n', xlab='cortical depth', ylab='% BOLD (contra-lateral)', axes=FALSE, pch=19, type='l', lwd=2, lty=1, cex.lab=1.5, ylim=c( min( seMin ), max( seMax )  ) )
  xAxis <- seq(1,length(aggData))
  lines( xAxis, seMin, col = 'grey')
  lines( xAxis, seMax, col = 'grey')  
  polygon(c(xAxis, rev(xAxis)), c(seMax, rev(seMin)),
          col = rgb(0.5,0.5,0.5,0.8), border = NA)  
  #segments( seq(1,length(aggData)), seMin, seq(1,length(aggData)), seMax )  
  axis(1, seq(1,length(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
  axis(2, round( seq(min(aggData),max(aggData),length.out=3), 1 ), 
       round( seq(min(aggData),max(aggData),length.out=3), 0 ), las=0, cex.axis=1.5 )
  storeBoldWhole[ nPart, ,3 ] <- aggData
  
  
  
  library(quantmod)
  bootPeakFun <- function( x, index, sparParam ) {
    anatData <- x[index,]  
    mAnatData <- apply( anatData, 2, median )
    xVar <- seq( 0, 1, length.out=length(mAnatData) )
    xVarInterp <- seq( 0, 1, length.out=length(mAnatData)*2 )
    spAnat <- spline( xVar, mAnatData, xout=xVarInterp  )
    sm_spAnat <- smooth.spline( spAnat$x, spAnat$y, spar = sparParam )  
    #indexPeak <- which( diff( sm_spAnat$y ) < 0 )
    #indexPeak <- findPeaks( sm_spAnat$y, thresh = 0.005 ) 
    indexPeak <- which( sm_spAnat$y == max( sm_spAnat$y ) )
    if ( length(indexPeak)>0 ) { estPeak <- sm_spAnat$x[ indexPeak[1] ]  }
    if ( length(indexPeak)==0 ) { estPeak <- 999  }
    return( estPeak )
   }
  funData <- intensityFilt[rankProfilesBinNumeric==1,xVar,1]
  bootPeakFunData <- boot( funData, bootPeakFun, sparParam = 0.2, R=nRepetitions )
  funBoot <- quantile( bootPeakFunData$t[ bootPeakFunData$t!=999 ], probs=c( 0.025, 0.5, 0.975 ) )
  funBoot
#funBoot[1] <- funBoot[2]-sd(bootPeakFunData$t[ bootPeakFunData$t!=999 ])/2
  #funBoot[3] <- funBoot[2]+sd(bootPeakFunData$t[ bootPeakFunData$t!=999 ])/2

  out <- boot_fun( intensityFilt[,xVar,1], sample( dim( intensityFilt[,xVar,1] )[1] ), rankProfilesBinNumeric )
  bootMod <- boot( intensityFilt[,xVar,1], boot_fun, splitVar=rankProfilesBinNumeric, R=nRepetitions )
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )  
  
  #bootCi01 <- apply( bootMod$t, 2 , mean ) -  apply( bootMod$t, 2 , sd )
  #bootCi02 <- apply( bootMod$t, 2 , mean )
  #bootCi03 <- apply( bootMod$t, 2 , mean ) +  apply( bootMod$t, 2 , sd )
  #bootCi <- rbind( bootCi01, bootCi02, bootCi03 )  

  mVal <- cbind( bootCi[2,1:16], bootCi[2,17:32], bootCi[2,33:48] )
  seValmin <- rbind( bootCi[1,1:16], bootCi[1,17:32], bootCi[1,33:48] )
  seValmax <- rbind( bootCi[3,1:16], bootCi[3,17:32], bootCi[3,33:48] )
  xVar <- seq( 1, dim(seValmin)[2] )
  matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, col=c('blue','darkgreen','red'), xlim=c(min(xVar),max(xVar)), ylim=c( round( c(min(seValmin),max(seValmax)), 0 ) ) )
  lines( xVar, seValmin[1,], col = rgb(0,0,0.8,0.2) )
  lines( xVar, seValmax[1, ], col = rgb(0,0,0.8,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[1,], rev(seValmin[1,])),
          col = rgb(0,0,0.8,0.2), border = NA)  
  lines( xVar, seValmin[2,], col = rgb(0,0.8,0,0.2) )
  lines( xVar, seValmax[2, ], col = rgb(0,0.8,0,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[2,], rev(seValmin[2,])),
          col = rgb(0,0.8,0,0.2), border = NA)  
  lines( xVar, seValmin[3,], col = rgb(0.8,0,0,0.2) )
  lines( xVar, seValmax[3, ], col = rgb(0.8,0,0,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[3,], rev(seValmin[3,])),
          col = rgb(0.8,0,0,0.2), border = NA)    
  #  segments( xVar+0.15, seValmin[1,], xVar+0.15, seValmax[1,], col=c('blue')  )
  #  segments( xVar, seValmin[2,], xVar, seValmax[2,], col=c('darkgreen')  )
  #  segments( xVar-0.15, seValmin[3,], xVar-0.15, seValmax[3,], col=c('red')  )
  axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
  axis(2, round( seq(min(mVal),max(mVal),length.out=3), 0 ), 
       round( seq(min(mVal),max(mVal),length.out=3), 0 ), las=1, cex.axis=1.5 )
  abline( v=funBoot[1]*max(xVar), lwd=2, lty=2, col='gray50' )
  abline( v=funBoot[2]*max(xVar), lwd=2, lty=1, col='gray50' )
  abline( v=funBoot[3]*max(xVar), lwd=2, lty=2, col='gray50' )  
  storeBold[nPart, ,1,] <- mVal 
  storeBoldEst[nPart,] <- c(  funBoot[1], funBoot[2], funBoot[3] )
  
  out <- boot_fun( intensityFilt[,xVar,2], sample( dim( intensityFilt[,xVar,1] )[1] ), rankProfilesBinNumeric )
  bootMod <- boot( intensityFilt[,xVar,2], boot_fun, splitVar=rankProfilesBinNumeric, R=nRepetitions )
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )  
  mVal <- cbind( bootCi[2,1:16], bootCi[2,17:32], bootCi[2,33:48] )
  seValmin <- rbind( bootCi[1,1:16], bootCi[1,17:32], bootCi[1,33:48] )
  seValmax <- rbind( bootCi[3,1:16], bootCi[3,17:32], bootCi[3,33:48] )
  xVar <- seq( 1, dim(seValmin)[2] )
  matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, col=c('blue','darkgreen','red'), ylim=c( round( c(min(seValmin),max(seValmax)), 2 ) ) )
  lines( xVar, seValmin[1,], col = rgb(0,0,0.8,0.2) )
  lines( xVar, seValmax[1, ], col = rgb(0,0,0.8,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[1,], rev(seValmin[1,])),
          col = rgb(0,0,0.8,0.2), border = NA)  
  lines( xVar, seValmin[2,], col = rgb(0,0.8,0,0.2) )
  lines( xVar, seValmax[2, ], col = rgb(0,0.8,0,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[2,], rev(seValmin[2,])),
          col = rgb(0,0.8,0,0.2), border = NA)  
  lines( xVar, seValmin[3,], col = rgb(0.8,0,0,0.2) )
  lines( xVar, seValmax[3, ], col = rgb(0.8,0,0,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[3,], rev(seValmin[3,])),
          col = rgb(0.8,0,0,0.2), border = NA)    
  #  segments( xVar+0.15, seValmin[1,], xVar+0.15, seValmax[1,], col=c('blue')  )
  #  segments( xVar, seValmin[2,], xVar, seValmax[2,], col=c('darkgreen')  )
  #  segments( xVar-0.15, seValmin[3,], xVar-0.15, seValmax[3,], col=c('red')  )
  axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
  axis(2, round( seq(min(seValmin),max(seValmax),length.out=3), 2 ), 
       round( seq(min(seValmin),max(seValmax),length.out=3), 2 ), las=1, cex.axis=1.5 )
  storeBold[nPart, ,2,] <- mVal 
  
  out <- boot_fun( intensityFilt[,xVar,5], sample( dim( intensityFilt[,xVar,1] )[1] ), rankProfilesBinNumeric )
  bootMod <- boot( intensityFilt[,xVar,5], boot_fun, splitVar=rankProfilesBinNumeric, R=nRepetitions )
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )  
  mVal <- cbind( bootCi[2,1:16], bootCi[2,17:32], bootCi[2,33:48] )
  seValmin <- rbind( bootCi[1,1:16], bootCi[1,17:32], bootCi[1,33:48] )
  seValmax <- rbind( bootCi[3,1:16], bootCi[3,17:32], bootCi[3,33:48] )
  xVar <- seq( 1, dim(seValmin)[2] )
  matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, col=c('blue','darkgreen','red'), ylim=c( round( c(min(seValmin),max(seValmax)), 0 ) ) )
  lines( xVar, seValmin[1,], col = rgb(0,0,0.8,0.2) )
  lines( xVar, seValmax[1, ], col = rgb(0,0,0.8,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[1,], rev(seValmin[1,])),
          col = rgb(0,0,0.8,0.2), border = NA)  
  lines( xVar, seValmin[2,], col = rgb(0,0.8,0,0.2) )
  lines( xVar, seValmax[2, ], col = rgb(0,0.8,0,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[2,], rev(seValmin[2,])),
          col = rgb(0,0.8,0,0.2), border = NA)  
  lines( xVar, seValmin[3,], col = rgb(0.8,0,0,0.2) )
  lines( xVar, seValmax[3, ], col = rgb(0.8,0,0,0.2) )  
  polygon(c(xAxis, rev(xAxis)), c(seValmax[3,], rev(seValmin[3,])),
          col = rgb(0.8,0,0,0.2), border = NA)    
  #  segments( xVar+0.15, seValmin[1,], xVar+0.15, seValmax[1,], col=c('blue')  )
  #  segments( xVar, seValmin[2,], xVar, seValmax[2,], col=c('darkgreen')  )
  #  segments( xVar-0.15, seValmin[3,], xVar-0.15, seValmax[3,], col=c('red')  )
  axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
  axis(2, round( seq(min(mVal),max(mVal),length.out=3), 0 ), 
       round( seq(min(mVal),max(mVal),length.out=3), 0 ), las=0, cex.axis=1.5 )
  storeBold[nPart, ,3,] <- mVal 

  
  mydata <- intensityFilt[,xVar,7]
  library(cluster)  
  clust.out <- pam(mydata, 2)
  prop.table( table( clust.out$clustering ) )  
  boot_fun <- function( x, indices, splitVar ) {
    mVal <- aggregate( x[indices,], by=list(splitVar[indices]), FUN=median)    
    mVal <- t( as.matrix( mVal[,2:dim(mVal)[2]] ) )
    
    xVar <- seq( 0, 1, length.out=length(mVal) )
    sp <- smooth.spline( xVar, mVal, spar=0.1 )
    xOut <- seq(0,1,length.out=50)
    spInt <- predict(sp,xOut)
    
    indexMin <- which( diff( spInt$y )>0 )
    if ( length( indexMin ) > 0 ) { peakVal <- spInt$x[ indexMin[ length( indexMin ) ] ]  }
    if ( length( indexMin ) == 0 ) { peakVal <- 999  }
    
    return( c( mVal, peakVal ) )
  }
  bootMod <- boot( intensityFilt[,xVar,7], boot_fun, splitVar=rep(1,dim(intensityFilt[,xVar,7])[1]), R=nRepetitions )  
  bootCi <- apply( bootMod$t, 2, quantile, probs=c(0.025, 0.5, 0.975) )    
  bigAverage <- bootCi[2,1:16]
  #aggData <- aggregate(mydata,by=list(clust.out$clustering),FUN=mean)    
  #bigAverage <- apply( mydata, 2, median )
  bigAverage <- bigAverage[2:length(bigAverage)]
  xValues <- seq(0,1,length.out=length(bigAverage))
  plot( bigAverage~xValues, 
        bty='n', las=1, xlab='cortical depth',
        ylab='T1 intensity', pch=20, cex=1.5, type='l', lwd=2, lty=1,
        axes=FALSE, xlim=c(0,1), ylim=c(min( bigAverage ), max( bigAverage )), cex.lab=1.5   )
  lines( xValues, bootCi[1,2:16], col = rgb(0.5,0.5,0.5,0.8) )
  lines( xValues, bootCi[3,2:16], col = rgb(0.5,0.5,0.5,0.8) )  
  polygon(c(xValues, rev(xValues)), c(bootCi[3,2:16], rev(bootCi[1,2:16])),
          col = rgb(0.5,0.5,0.5,0.8), border = NA)  
  
  #  segments( xValues, bootCi[1,2:16], xValues, bootCi[3,2:16] )  
  axis( 1, c(0,0.5,1), c(0,0.5,1), cex.axis=1.5 )
  axis( 2, round( seq( min( bigAverage ), max( bigAverage ),length.out=3 ), 0), cex.axis=1.5 )
  storeT1[nPart,] <- bigAverage
  abline(  v=bootCi[1,17], lwd=2, lty=2, col='gray50')
  abline(  v=bootCi[2,17], lwd=2, lty=1, col='gray50')
  abline(  v=bootCi[3,17], lwd=2, lty=2, col='gray50')
  storeT1Est[nPart,] <- c( bootCi[1:3,17]  ) 
  #plot( as.numeric( aggData[ 1, 3:dim(aggData)[2] ] )~xValues[2:length(xValues)], bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=20, cex=1.5, type='l', lwd=2, lty=1  )
  #plot( as.numeric( aggData[ 2, 3:dim(aggData)[2] ] )~xValues[2:length(xValues)], bty='n', las=1, xlab='cortical depth', ylab='T1 intensity', pch=20, cex=1.5, type='l', lwd=2, lty=1  )
  
  pSteps <- seq(0.40,1,0.01)
  posProfileMatrix <- array( 0, c( length( pSteps ), 50 ) )
  for ( k in 1:length( pSteps ) ) {
    rankThr <- quantile( rankProfiles, probs=pSteps[k] )
    splitVar <- ifelse( rankProfiles<=rankThr, 1, 2 )
    
    posProfileTemp <- apply( intensityFilt[ splitVar==1, xVar, 1 ], 2, mean )
    sp <- smooth.spline( posProfileTemp, spar=0.25 )
    xOut <- seq(1,length(xVar),length.out=50)
    spInt <- predict(sp,xOut)
    
    posProfileMatrix[k,] <- scaleData( spInt$y, 1, 0 )  
  }  
  image( t( posProfileMatrix[,4:50] ), col=rainbow(1000, start=0.2, end=0.85), las=1, axes=FALSE, ylab='percentile', xlab='cortical depth', cex.lab=1.5 )
  axis( 1, c(0,0.5,1), c(0,0.5,1), cex.axis=1.5 )
  axis( 2, c(0,0.5,1), round( seq(min(pSteps),1,length.out=3), 2)*100, cex.axis=1.5, las=1 )
  
  if ( savePlots==1 ) {
    dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'boldCorticalThicknessRank_all' ), height=7.5, width=7)
    system( sprintf( 'convert -density 500 %s %s', 
                     sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'boldCorticalThicknessRank_all' ), 
                     sprintf('%s/%s_%s.png', mainDirSavePlots, participant ,'boldCorticalThicknessRank_all' )  ) )
  }
  dev.off()
  
  
  ## write out ROI as a map (to double check the mapping)
  nBins <- length( unique( rankProfilesBinNumeric ) )
  rankProfilesBinNumericFlag <- array( 0, c( length( rankProfilesBinNumeric ), nBins ) )
  for ( k in 1:nBins ) {
    rankProfilesBinNumericFlag[ rankProfilesBinNumeric==k, k ] <- 1
  }
  emptyMap <- array( 0, c( dim( mapValListStat[[10]] )[1], 5+nBins ) )
  emptyMap[,1] <- seq( 0, dim(emptyMap)[1]-1 )
  emptyMap[ surfValFilt+1, 2 ] <- 1
  emptyMap[ surfValFilt+1, 3 ] <- rankProfilesBinNumeric
  emptyMap[ surfValFilt+1, 4 ] <- rankProfiles
  emptyMap[ surfValFilt+1, 5 ] <- tZeroFilt
  emptyMap[ surfValFilt+1, seq(6,5+nBins) ] <- rankProfilesBinNumericFlag
  setwd(mainDir)
  write.table( emptyMap, file='storeMap01.1D', row.names=FALSE, col.names=FALSE )
  
  
  
     ## shift data accordingly
     #nangleShift <- 0
     angles <- intensityFilt[,,9] + angleShift
     angles <- ifelse( angles>180, angles-180, angles )
     angles <- ifelse( angles<0, angles+180, angles )  

     meanAngles <- apply(angles, 1, mean)
     storeAngleCor[nPart] <- cor( rankProfiles, meanAngles )     
     
      ## plot
     x11(height=7.5, width=7.5)
     par(mfrow=c(4,4))
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
       
       boldAngleArray[ nPart, , k, 1 ] <- bold_angle
       boldAngleArray[ nPart, , k, 2 ] <- rep( sModCos$r.squared, length(bold_angle)  )
       
     }
     #dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'raw_bold_positive_angle' ), height=7.5, width=7.5)
     #dev.off()
  #   
  #   
  #   
  #   x11(height=7.5, width=7.5)
  #   par(mfrow=c(4,4))
  #   xValues <- seq(0,1,length.out=length(xVar))
  #   storeValues <- array(0,c(length(xVar),1))
  #   for (k in 1:length(xVar)) {
  #     surfNum <- k
  #     B0angleBreaks <- cut( angles[,surfNum], breaks = quantile( angles[,surfNum], seq(0,1,0.1) ), include.lowest=TRUE )
  #     #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
  #     bold_angle <- tapply( intensityFilt[,surfNum,2], list( B0angleBreaks ), median  )
  #     B0_angle <- tapply( angles[,surfNum], list( B0angleBreaks ), median  )
  #     B0_angle2 <- B0_angle^2 
  #     cosPred <- cos( (B0_angle)*pi/180 )^2
  #     linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
  #     modCos <- lm( bold_angle ~ cosPred ) 
  #     sModCos <- summary( modCos )
  #     predMod <- predict( modCos, data.frame( cosPred ) )
  #     
  #     #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
  #     #predMod01 <- predict( modCos, data.frame( x01 ) )  
  #     
  #     plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1 )    
  #     interpData <- spline( B0_angle, predMod, n=50, xout=seq(40,140,10) )
  #     
  #     #lines( B0_angle, predMod, lwd=2, col='red' )
  #     lines( interpData$x, interpData$y, lwd=2, col='red' )
  #     
  #     #lines( seq(40,170,length.out=length(cosPred)), predMod01, lwd=2, col='blue' )
  #     storeValues[k] <- sModCos$r.squared 
  #     
  #     boldAngleArray[ nPart, , k, 3 ] <- bold_angle    
  #     boldAngleArray[ nPart, , k, 4 ] <- rep( sModCos$r.squared, length(bold_angle)  )
  #   }
  #   dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'raw_bold_negative_angle' ), height=7.5, width=7.5)
  #   #dev.off()
  #   
  #   
  #   
  #   x11(height=7.5, width=7.5)
  #   par(mfrow=c(4,4))
  #   xValues <- seq(0,1,length.out=length(xVar))
  #   storeValues <- array(0,c(length(xVar),1))
  #   for (k in 1:length(xVar)) {
  #     surfNum <- k
  #     B0angleBreaks <- cut( angles[,surfNum], breaks = quantile( angles[,surfNum], seq(0,1,0.1) ), include.lowest=TRUE )
  #     #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
  #     bold_angle <- tapply( intensityFilt[,surfNum,5], list( B0angleBreaks ), median  )
  #     B0_angle <- tapply( angles[,surfNum], list( B0angleBreaks ), median  )
  #     B0_angle2 <- B0_angle^2 
  #     cosPred <- cos( (B0_angle)*pi/180 )^2
  #     linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
  #     modCos <- lm( bold_angle ~ cosPred ) 
  #     sModCos <- summary( modCos )
  #     predMod <- predict( modCos, data.frame( cosPred ) )
  #     
  #     #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
  #     #predMod01 <- predict( modCos, data.frame( x01 ) )  
  #     
  #     plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1 )    
  #     interpData <- spline( B0_angle, predMod, n=50, xout=seq(40,140,10) )
  #     
  #     #lines( B0_angle, predMod, lwd=2, col='red' )
  #     lines( interpData$x, interpData$y, lwd=2, col='red' )
  #     
  #     #lines( seq(40,170,length.out=length(cosPred)), predMod01, lwd=2, col='blue' )
  #     storeValues[k] <- sModCos$r.squared 
  #     
  #     boldAngleArray[ nPart, , k, 5 ] <- bold_angle
  #     boldAngleArray[ nPart, , k, 6 ] <- scale( bold_angle, center=TRUE, scale=TRUE )
  #     boldAngleArray[ nPart, , k, 7 ] <- rep( sModCos$r.squared, length(bold_angle)  )
  #   }
  #   dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'raw_bold_angle' ), height=7.5, width=7.5)
  #   #dev.off()
  #   
  #   
  
  #   medianB0Angle <- apply( angles, 1, median )
  #   rankProfilesBin <- cut( medianB0Angle, breaks=seq(0,180,45), include.lowest=TRUE )
  #   rankProfilesBinNumeric <- as.numeric( rankProfilesBin )
  #   
  #   medianBold <- apply( intensityFilt[,,1], 1, median  )
  #   tapply( medianBold, list(rankProfilesBinNumeric), median )
  #   colList <- c('red','green','lightblue','orange')
  #   
  #   x11(width=6.5, height=2.5)
  #   par(mfrow=c(1,2))
  #   for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #     if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ),
  #                       xlim=c( min(xVar),max(xVar) ), ylim=c(1,8), bty='n', col=colList[k], type='l', axes=FALSE,
  #                       ylab='% BOLD (contralateral)', xlab='cortical depth')
  #               axis( 1, seq( min(xVar),max(xVar), length.out=3 ),  seq(0,1,0.5) )
  #               axis( 2,seq( 1, 8, length.out=3 ) ) }
  #     if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), col=colList[k] ) }
  #     depthArrayAngle[nPart,,1,k] <- apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median )
  #   }
  #   for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #     if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ),
  #                       xlim=c(min(xVar),max(xVar)), ylim=c(-2,0.5), bty='n', col=colList[k], type='l', axes=FALSE,
  #                       ylab='% BOLD (ipsilateral)', xlab='cortical depth')
  #                 axis( 1, seq( min(xVar),max(xVar), length.out=3 ),  seq(0,1,0.5) )
  #                 axis( 2,seq( -2, 0.5, length.out=3 ) ) }
  #     if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), col=colList[k] ) }  
  #     depthArrayAngle[nPart,,2,k] <- apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median )
  #   }
  #   for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #     if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ], 2, median ),
  #                       xlim=c(min(xVar),max(xVar)), ylim=c(0,4), bty='n', col=colList[k], type='l', axes=FALSE,
  #                       ylab='T stat (ipsilateral)', xlab='cortical depth')
  #                 axis( 1, seq( min(xVar),max(xVar), length.out=3 ),  seq(0,1,0.5) )
  #                 axis( 2,seq( 0, 4, length.out=3 ) ) }
  #     if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ], 2, median ), col=colList[k] ) }  
  #     depthArrayAngle[nPart,,3,k] <- apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ], 2, median )
  #   }
  #   dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'bold_corticalThickness_angle' ), height=2.5, width=4.5)
  #   dev.off()
  #   
  #   
  #   rankProfilesBinNumeric <- rep(1,length(rankProfilesBinNumeric))
  #   
  #   x11(width=4.5, height=2.5)
  #   par(mfrow=c(1,2))
  #   for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #     if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ),
  #                       xlim=c( min(xVar),max(xVar) ), ylim=c(1,8), bty='n', col=colList[k], type='l', axes=FALSE,
  #                       ylab='% BOLD (contralateral)', xlab='cortical depth')
  #                 axis( 1, seq( min(xVar),max(xVar), length.out=3 ),  seq(0,1,0.5) )
  #                 axis( 2,seq( 1, 8, length.out=3 ) ) }
  #     if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median ), col=colList[k] ) }  
  #   }
  #   for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #     if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ),
  #                       xlim=c(min(xVar),max(xVar)), ylim=c(-2,0.5), bty='n', col=colList[k], type='l', axes=FALSE,
  #                       ylab='% BOLD (ipsilateral)', xlab='cortical depth')
  #                 axis( 1, seq( min(xVar),max(xVar), length.out=3 ),  seq(0,1,0.5) )
  #                 axis( 2,seq( -2, 0.5, length.out=3 ) ) }
  #     if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median ), col=colList[k] ) }  
  #   }
  #   for (k in 1:length(unique(rankProfilesBinNumeric))) {
  #     if (k==1) { plot( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ], 2, median ),
  #                       xlim=c(min(xVar),max(xVar)), ylim=c(0,4), bty='n', col=colList[k], type='l', axes=FALSE,
  #                       ylab='T stat (ipsilateral)', xlab='cortical depth')
  #                 axis( 1, seq( min(xVar),max(xVar), length.out=3 ),  seq(0,1,0.5) )
  #                 axis( 2,seq( 0, 4, length.out=3 ) ) }
  #     if (k>1)  { lines( xVar, apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ], 2, median ), col=colList[k] ) }  
  #   }
  #   depthArray[nPart,,1] <- apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 1 ], 2, median )
  #   depthArray[nPart,,2] <- apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 2 ], 2, median )
  #   depthArray[nPart,,3] <- apply( intensityFilt[ rankProfilesBinNumeric==k, xVar, 6 ], 2, median )
  #   
  #   dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, participant ,'bold_corticalThickness' ), height=2.5, width=4.5)
  #   dev.off()
  
  
  
  
}

mean( storeMotion )
sd( storeMotion )




x11( height=7.5, width=7 )
par(mfrow=c(3,3))

mVal <- apply( storeBoldWhole[,,1], 2, mean )
#mVal <- mVal[,c(3,2,1)]
seVal <- ( apply( storeBoldWhole[,,1], 2, sd ) / sqrt( length(participantArray) ) ) / 2 
seValmin <- mVal - seVal
seValmax <- mVal + seVal
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, ylim=round( c(min(seValmin),max(seValmax)), 0 ) ) 
#segments( xVar+0.15, mVal-seVal, xVar+0.18, mVal+seVal, col=c('black')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), las=1, cex.axis=1.5 )
lines( xVar, seValmin, col = rgb(0.5,0.5,0.5,0.8) )
lines( xVar, seValmax, col = rgb(0.5,0.5,0.5,0.8) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax, rev(seValmin)),
        col = rgb(0.5,0.5,0.5,0.2), border = NA)    

mVal <- apply( storeBoldWhole[,,2], 2, mean )
#mVal <- mVal[,c(3,2,1)]
seVal <- ( apply( storeBoldWhole[,,1], 2, sd ) / sqrt( length(participantArray) ) ) / 2 
seValmin <- mVal - seVal
seValmax <- mVal + seVal
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, ylim=round( c(min(seValmin),max(seValmax)), 0 ) ) 
#segments( xVar+0.15, mVal-seVal, xVar+0.18, mVal+seVal, col=c('black')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), las=1, cex.axis=1.5 )
lines( xVar, seValmin, col = rgb(0.5,0.5,0.5,0.8) )
lines( xVar, seValmax, col = rgb(0.5,0.5,0.5,0.8) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax, rev(seValmin)),
        col = rgb(0.5,0.5,0.5,0.2), border = NA)    

mVal <- apply( storeBoldWhole[,,3], 2, mean )
#mVal <- mVal[,c(3,2,1)]
seVal <- ( apply( storeBoldWhole[,,3], 2, sd ) / sqrt( length(participantArray) ) ) / 2 
seValmin <- mVal - seVal
seValmax <- mVal + seVal
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='raw BOLD', type='l', lwd=2, lty=1, cex.lab=1.5, ylim=round( c(min(seValmin),max(seValmax)), 0 ) ) 
#segments( xVar+0.15, mVal-seVal, xVar+0.18, mVal+seVal, col=c('black')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), las=0, cex.axis=1.5 )
lines( xVar, seValmin, col = rgb(0.5,0.5,0.5,0.8) )
lines( xVar, seValmax, col = rgb(0.5,0.5,0.5,0.8) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax, rev(seValmin)),
        col = rgb(0.5,0.5,0.5,0.2), border = NA)    

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

bootPeakFun <- function( x, index, sparParam ) {
  anatData <- x[index,]  
  mAnatData <- apply( anatData, 2, mean )
  xVar <- seq( 0, 1, length.out=length(mAnatData) )
  xVarInterp <- seq( 0, 1, length.out=length(mAnatData)*5 )
  spAnat <- spline( xVar, mAnatData, xout=xVarInterp  )
  sm_spAnat <- smooth.spline( spAnat$x, spAnat$y, spar = sparParam )  
  #indexPeak <- which( diff( sm_spAnat$y ) < 0 )
  indexPeak <- which(  sm_spAnat$y == max( sm_spAnat$y )  )
  if ( length(indexPeak)>0 ) { estPeak <- sm_spAnat$x[ indexPeak[1] ]  }
  if ( length(indexPeak)==0 ) { estPeak <- 999  }
  return( estPeak )
}

funData <- storeBold[,,1,1]
bootPeakFun <- boot( funData, bootPeakFun, sparParam = 0.01, R=nRepetitions )
funBoot <- quantile( bootPeakFun$t, probs=c( 0.025, 0.5, 0.975 ) )

mVal <- apply( storeBold[,,1,], c(2,3), mean )
seValmin <- t( mVal - ( apply( storeBold[,,1,], c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 )
seValmax <- t( mVal + ( apply( storeBold[,,1,], c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 )
xVar <- seq( 1, dim(seValmin)[2] )
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, col=c('blue','darkgreen','red'), ylim=c( round( c(min(seValmin),max(seValmax)), 0 ) ) )
lines( xVar, seValmin[1,], col = rgb(0,0,0.8,0.2) )
lines( xVar, seValmax[1,], col = rgb(0,0,0.8,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[1,], rev(seValmin[1,])),
        col = rgb(0,0,0.8,0.2), border = NA)  
lines( xVar, seValmin[2,], col = rgb(0,0.8,0,0.2) )
lines( xVar, seValmax[2, ], col = rgb(0,0.8,0,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[2,], rev(seValmin[2,])),
        col = rgb(0,0.8,0,0.2), border = NA)  
lines( xVar, seValmin[3,], col = rgb(0.8,0,0,0.2) )
lines( xVar, seValmax[3, ], col = rgb(0.8,0,0,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[3,], rev(seValmin[3,])),
        col = rgb(0.8,0,0,0.2), border = NA)    
abline( v=mean(storeBoldEst[,1])*max(xVar), lwd=3, lty=2, col='gray50' )
abline( v=mean(storeBoldEst[,2])*max(xVar), lwd=3, lty=1, col='gray50' )
abline( v=mean(storeBoldEst[,3])*max(xVar), lwd=3, lty=2, col='gray50' )
#  segments( xVar+0.15, seValmin[1,], xVar+0.15, seValmax[1,], col=c('blue')  )
#  segments( xVar, seValmin[2,], xVar, seValmax[2,], col=c('darkgreen')  )
#  segments( xVar-0.15, seValmin[3,], xVar-0.15, seValmax[3,], col=c('red')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq(min(mVal),max(mVal),length.out=3), 0 ), 
     round( seq(min(mVal),max(mVal),length.out=3), 0 ), las=1, cex.axis=1.5 )


mVal <- apply( storeBold[,,2,], c(2,3), mean )
seValmin <- t( mVal - ( apply( storeBold[,,2,], c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 )
seValmax <- t( mVal + ( apply( storeBold[,,2,], c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 )
xVar <- seq( 1, dim(seValmin)[2] )
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, col=c('blue','darkgreen','red'), ylim=c( round( c(min(seValmin),max(seValmax)), 0 ) ) )
lines( xVar, seValmin[1,], col = rgb(0,0,0.8,0.2) )
lines( xVar, seValmax[1,], col = rgb(0,0,0.8,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[1,], rev(seValmin[1,])),
        col = rgb(0,0,0.8,0.2), border = NA)  
lines( xVar, seValmin[2,], col = rgb(0,0.8,0,0.2) )
lines( xVar, seValmax[2, ], col = rgb(0,0.8,0,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[2,], rev(seValmin[2,])),
        col = rgb(0,0.8,0,0.2), border = NA)  
lines( xVar, seValmin[3,], col = rgb(0.8,0,0,0.2) )
lines( xVar, seValmax[3, ], col = rgb(0.8,0,0,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[3,], rev(seValmin[3,])),
        col = rgb(0.8,0,0,0.2), border = NA)    
#  segments( xVar+0.15, seValmin[1,], xVar+0.15, seValmax[1,], col=c('blue')  )
#  segments( xVar, seValmin[2,], xVar, seValmax[2,], col=c('darkgreen')  )
#  segments( xVar-0.15, seValmin[3,], xVar-0.15, seValmax[3,], col=c('red')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq(min(mVal),max(mVal),length.out=3), 0 ), 
     round( seq(min(mVal),max(mVal),length.out=3), 0 ), las=1, cex.axis=1.5 )


storeZscore <- ( storeBold[,,3,] - mean( storeBold[,,3,] ) ) / sd( storeBold[,,3,] )
mVal <-  apply( storeZscore, c(2,3), mean )
seValmin <- t( mVal - ( apply( storeZscore, c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 )
seValmax <- t( mVal + ( apply( storeZscore, c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 )
xVar <- seq( 1, dim(seValmin)[2] )
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, col=c('blue','darkgreen','red'), ylim=c( round( c(min(seValmin),max(seValmax)), 2 ) ) )
lines( xVar, seValmin[1,], col = rgb(0,0,0.8,0.2) )
lines( xVar, seValmax[1,], col = rgb(0,0,0.8,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[1,], rev(seValmin[1,])),
        col = rgb(0,0,0.8,0.2), border = NA)  
lines( xVar, seValmin[2,], col = rgb(0,0.8,0,0.2) )
lines( xVar, seValmax[2, ], col = rgb(0,0.8,0,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[2,], rev(seValmin[2,])),
        col = rgb(0,0.8,0,0.2), border = NA)  
lines( xVar, seValmin[3,], col = rgb(0.8,0,0,0.2) )
lines( xVar, seValmax[3, ], col = rgb(0.8,0,0,0.2) )  
polygon(c(xAxis, rev(xAxis)), c(seValmax[3,], rev(seValmin[3,])),
        col = rgb(0.8,0,0,0.2), border = NA)    
#  segments( xVar+0.15, seValmin[1,], xVar+0.15, seValmax[1,], col=c('blue')  )
#  segments( xVar, seValmin[2,], xVar, seValmax[2,], col=c('darkgreen')  )
#  segments( xVar-0.15, seValmin[3,], xVar-0.15, seValmax[3,], col=c('red')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq(min(seValmin),max(seValmax),length.out=3), 2 ), 
     round( seq(min(seValmin),max(seValmax),length.out=3), 2 ), las=1, cex.axis=1.5 )

bootPeak <- function( x, index, sparParam ) {
  anatData <- x[index,]  
  mAnatData <- apply( anatData, 2, mean )
  xVar <- seq( 0, 1, length.out=length(mAnatData) )
  xVarInterp <- seq( 0, 1, length.out=length(mAnatData)*5 )
  spAnat <- spline( xVar, mAnatData, xout=xVarInterp  )
  sm_spAnat <- smooth.spline( spAnat$x, spAnat$y, spar = sparParam )  
  indexPeak_0 <- which( diff( sm_spAnat$y ) > 0 )
  indexPeak_1 <- which( diff( sm_spAnat$y ) < 0 )
  #indexPeak <- which( indexPeak_0 > indexPeak_1 )
  if ( length(indexPeak_0)>0 ) { estPeak <- sm_spAnat$x[ indexPeak_0[length(indexPeak_0)] ]  }
  if ( length(indexPeak_0)==0 ) { estPeak <- 999  }
  return( estPeak )
}

bootPeakAnat <- boot( storeT1, bootPeak, sparParam = 0.01, R=nRepetitions )
anatBoot <- quantile( bootPeakAnat$t, probs=c( 0.025, 0.5, 0.975 ) )

mVal <- apply( storeT1, 2, mean )
#mVal <- mVal[,c(3,2,1)]
seVal <- ( apply( storeT1, 2, sd ) / sqrt( length(participantArray) ) ) / 2 
seValmin <- mVal - seVal
seValmax <- mVal + seVal
xVarAnat <- xVar[1:length(mVal)]
matplot( xVarAnat, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='T1 intensity', type='l', lwd=2, lty=1, cex.lab=1.5, ylim=round( c(min(seValmin),max(seValmax)), 0 ) ) 
#segments( xVar+0.15, mVal-seVal, xVar+0.18, mVal+seVal, col=c('black')  )
axis(1, seq(min(xVarAnat),max(xVarAnat),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), round( seq( min(seValmin),max(seValmax), length.out=3 ), 0), las=0, cex.axis=1.5 )
lines( xVarAnat, seValmin, col = rgb(0.5,0.5,0.5,0.2) )
lines( xVarAnat, seValmax, col = rgb(0.5,0.5,0.5,0.2) )  
polygon(c(xVarAnat, rev(xVarAnat)), c(seValmax, rev(seValmin)),
        col = rgb(0.5,0.5,0.5,0.2), border = NA)    
abline( v=mean( storeT1Est[,1] )*max(xVarAnat), lwd=3, lty=2, col='gray50' )
abline( v=mean( storeT1Est[,2] )*max(xVarAnat), lwd=3, lty=1, col='gray50' )
abline( v=mean( storeT1Est[,3] )*max(xVarAnat), lwd=3, lty=2, col='gray50' )

if (savePlots==1) {
  dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, 'average' ,'bold_corticalThickness_all' ), height=7.5, width=7)
}
dev.off()







mVal <- apply( storeBold[,,1,], c(2,3), mean )
#mVal <- mVal[,c(3,2,1)]
seVal <- ( apply( storeBold[,,1,], c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2 
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5 )
segments( xVar+0.15, mVal[,1]-seVal[,1], xVar+0.18, mVal[,1]+seVal[,1], col=c('black')  )
segments( xVar, mVal[,2]-seVal[,2], xVar, mVal[,2]+seVal[,2], col=c('red')  )
segments( xVar-0.15, mVal[,3]-seVal[,3], xVar-0.18, mVal[,3]+seVal[,3], col=c('darkgreen')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, round( seq(min(mVal),max(mVal),length.out=3), 0 ), 
     round( seq(min(mVal),max(mVal),length.out=3), 0 ), las=1, cex.axis=1.5 )

predQuad <- seq(-1,1,length.out=length(xVar))^2
predLin <- seq(-1,1,length.out=length(xVar))
summary( lm( mVal[,1] ~ predQuad ) )
summary( lm( mVal[,2] ~ predQuad ) )
summary( lm( mVal[,3] ~ predQuad ) )
summary( lm( mVal[,1] ~ predLin ) )
summary( lm( mVal[,2] ~ predLin ) )
summary( lm( mVal[,3] ~ predLin ) )


mVal <- apply( storeBold[,,2,], c(2,3), mean ) 
seVal <- ( apply( storeBold[,,2,], c(2,3), sd ) / sqrt( length(participantArray) ) ) / 2
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (ipsi-lateral)', type='l', lwd=2, lty=1, ylim=c(-1, 0), cex.lab=1.5 )
segments( xVar+0.18, mVal[,1]-seVal[,1], xVar+0.18, mVal[,1]+seVal[,1], col=c('black')  )
segments( xVar, mVal[,2]-seVal[,2], xVar, mVal[,2]+seVal[,2], col=c('red')  )
segments( xVar-0.18, mVal[,3]-seVal[,3], xVar-0.18, mVal[,3]+seVal[,3], col=c('darkgreen')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, c(0, -0.5, -1), c(0, -0.5, -1), las=1, cex.axis=1.5 )

predQuad <- seq(-1,1,length.out=length(xVar))^2
predLin <- seq(-1,1,length.out=length(xVar))
summary( lm( mVal[,1] ~ predQuad ) )
summary( lm( mVal[,2] ~ predQuad ) )
summary( lm( mVal[,3] ~ predQuad ) )
summary( lm( mVal[,1] ~ predLin ) )
summary( lm( mVal[,2] ~ predLin ) )
summary( lm( mVal[,3] ~ predLin ) )

summary( lm( apply( mVal, 1, mean) ~ predLin ) )


mVal <- scale( apply( storeBold[,,3,], c(2,3), mean ) ) 
seVal <- ( scale( apply( storeBold[,,3,], c(2,3), sd ) ) / sqrt( length(participantArray) ) ) / 2
matplot( xVar, mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='raw BOLD (z-score)', type='l', lwd=2, lty=1, ylim=c(-1.7, 1.7), cex.lab=1.5 )
segments( xVar+0.15, mVal[,1]-seVal[,1], xVar+0.18, mVal[,1]+seVal[,1], col=c('black')  )
segments( xVar, mVal[,2]-seVal[,2], xVar, mVal[,2]+seVal[,2], col=c('red')  )
segments( xVar-0.15, mVal[,3]-seVal[,3], xVar-0.18, mVal[,3]+seVal[,3], col=c('darkgreen')  )
axis(1, seq(min(xVar),max(xVar),length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, c(-1.5, 0, 1.5), c(-1.5, 0, 1.5), las=1, cex.axis=1.5 )

predQuad <- seq(-1,1,length.out=length(xVar))^2
predLin <- seq(-1,1,length.out=length(xVar))
summary( lm( mVal[,1] ~ predQuad + predLin ) )
summary( lm( mVal[,2] ~ predQuad + predLin ) )
summary( lm( mVal[,3] ~ predQuad + predLin ) )

mVal <- apply( storeT1, 2, mean )
seVal <- ( apply( storeT1, 2, sd ) / sqrt( length(participantArray) ) ) / 2 
plot( xVar[1:15], mVal, bty='n', pch=19, axes=FALSE, xlab='cortical depth', ylab='%BOLD (contra-lateral)', type='l', lwd=2, lty=1, cex.lab=1.5, xlim=c(1,15), ylim=c( min(mVal), max(mVal) ) )
segments( xVar+0.15, mVal-seVal, xVar+0.18, mVal+seVal, col=c('black')  )
axis(1, seq(1,15,length.out=3), seq(0,1,0.5), cex.axis=1.5 )
axis(2, seq( min(mVal), max(mVal), length.out=3 ), round( seq( min(mVal), max(mVal), length.out=3 ), 0), las=0, cex.axis=1.5 )

# 
# mVal <- apply( storeT1, 2, mean )
# dev.copy2pdf(file = sprintf('%s/%s_%s.pdf', mainDirSavePlots, 'average' ,'bold_corticalThickness' ), height=7.5, width=7.5)
# dev.off()
# 


x11(height=7.5, width=7.5)
par(mfrow=c(4,4))
xValues <- seq(0,1,length.out=length(xVar))
bold_angleStoreVar <- array( 0, c(length(xVar)) ) 
storeValues <- array(0,c(length(xVar),1))
for (k in 1:length(xVar)) {
  surfNum <- k
  #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
  bold_angleStoreVar[k] <- mean( apply( boldAngleArray[,,k,2], c(2), mean ) )
  bold_angle <- apply( boldAngleArray[,,k,1], c(2), mean )
  bold_angle_se <- apply( boldAngleArray[,,k,1], c(2), sd ) / sqrt( length(participantArray) )
  B0_angle <- seq( 0, 180, length.out=length(bold_angle) ) 
  B0_angle2 <- B0_angle^2 
  cosPred <- cos( (B0_angle)*pi/180 )^2
  linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
  modCos <- lm( bold_angle ~ cosPred ) 
  sModCos <- summary( modCos )
  predMod <- predict( modCos, data.frame( cosPred ) )
  
  #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
  #predMod01 <- predict( modCos, data.frame( x01 ) )  
  
  plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1, ylab='% BOLD', xlab=expression(paste(phi)), axes=FALSE,
        ylim=c( round( min(bold_angle), 1 ), round( max(bold_angle), 1 ) ), cex.lab=1.25 )    
  axis( 1, seq(0,180,length.out=3), cex.axis=1.25 )
  axis( 2, seq( round( min(bold_angle), 1 ), round( max(bold_angle), 1 ), length.out=3), cex.axis=1.25 )
  segments( B0_angle, bold_angle-bold_angle_se/2, B0_angle, bold_angle+bold_angle_se/2 )
  
  interpData <- spline( B0_angle, predMod, n=50, xout=seq(20,160,10) )
  
  #lines( B0_angle, predMod, lwd=2, col='red' )
  lines( interpData$x+5, interpData$y, lwd=2, col='red' )
    
}
if (savePlots==1) {
  dev.copy2pdf(file = sprintf('%s/%s.pdf', mainDirSavePlots, 'average_B0Angle_percent_BOLD_all' ), height=7.5, width=7.5)
}
dev.off()






x11(height=7.5, width=7.5)
par(mfrow=c(4,4))
xValues <- seq(0,1,length.out=length(xVar))
storeValues <- array(0,c(length(xVar),1))
for (k in 1:length(xVar)) {
  surfNum <- k
  #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
  bold_angle <- apply( boldAngleArray[,,k,5], c(2), mean )
  
  B0_angle <- seq( 0, 180, length.out=length(bold_angle) )
  B0_angle2 <- B0_angle^2 
  cosPred <- cos( (B0_angle)*pi/180 )^2
  linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
  modCos <- lm( bold_angle ~ cosPred ) 
  sModCos <- summary( modCos )
  predMod <- predict( modCos, data.frame( cosPred ) )
  
  #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
  #predMod01 <- predict( modCos, data.frame( x01 ) )  
  
  plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1, ylab='raw BOLD', xlab=expression(paste(phi)), axes=FALSE,
  ylim=c( floor( min(bold_angle) ), floor( max(bold_angle) ) ), cex.lab=1.25 )    
  axis( 1, seq(0,180,length.out=3), cex.axis=1.25  )
  axis( 2, seq( floor( min(bold_angle)  ), floor( max(bold_angle) ), length.out=2), cex.axis=1.25  )
  
  interpData <- spline( B0_angle, predMod, n=50, xout=seq(20,160,10) )
  
  #lines( B0_angle, predMod, lwd=2, col='red' )
  lines( interpData$x, interpData$y, lwd=2, col='red' )
  
}
dev.copy2pdf(file = sprintf('%s/%s.pdf', mainDirSavePlots, 'average_B0Angle_raw_BOLD_all' ), height=7.5, width=7.5)
dev.off()

x11(height=7.5, width=7.5)
par(mfrow=c(4,4))
xValues <- seq(0,1,length.out=length(xVar))
storeValues <- array(0,c(length(xVar),1))
for (k in 1:length(xVar)) {
  surfNum <- k
  #B0angleBreaks <- cut( intensityFilt[,surfNum,9], breaks = seq(0,180,30), include.lowest=TRUE )
  bold_angle <- apply( boldAngleArray[,,k,5], c(2), mean )
  bold_angle_se <- apply( boldAngleArray[,,k,5], c(2), sd ) / sqrt( length(participantArray) )
  B0_angle <- seq( 0, 180, length.out=length(bold_angle) )
  B0_angle2 <- B0_angle^2 
  cosPred <- cos( (B0_angle)*pi/180 )^2
  linPred <- seq( min(cosPred), max(cosPred), length.out=length(cosPred) )  
  modCos <- lm( bold_angle ~ cosPred ) 
  sModCos <- summary( modCos )
  predMod <- predict( modCos, data.frame( cosPred ) )
  
  #x01 <- cos( seq(40,170,length.out=length(cosPred))*pi/180 )^2
  #predMod01 <- predict( modCos, data.frame( x01 ) )  
  
  plot( bold_angle ~ B0_angle, xlim=c(0,180), bty='n', pch=15, cex=1, ylab='normalized BOLD', xlab=expression(paste(phi)), axes=FALSE,
        ylim=c( round( min(bold_angle), 1 ), round( max(bold_angle), 1 ) ), cex.lab=1.25  )    
  axis( 1, seq(0,180,length.out=3), cex.axis=1.25  )
  axis( 2, seq( round( min(bold_angle), 1 ), round( max(bold_angle), 1 ), length.out=3), cex.axis=1.25  )
  segments( B0_angle, bold_angle-bold_angle_se/2, B0_angle, bold_angle+bold_angle_se/2 )
  interpData <- spline( B0_angle, predMod, n=50, xout=seq(20,160,10) )
  
  #lines( B0_angle, predMod, lwd=2, col='red' )
  lines( interpData$x, interpData$y, lwd=2, col='red' )
  
}
dev.copy2pdf(file = sprintf('%s/%s.pdf', mainDirSavePlots, 'average_B0Angle_raw_BOLD_zscore_all' ), height=7.5, width=7.5)
dev.off()






x11(height=3.5, width=6)
par(mfrow=c(1,2))
xValues <- seq(0,1,length.out=length(xVar))
y <- apply( depthArray[,,1], 2, mean )
y_se <- apply( depthArray[,,1], 2, sd ) / sqrt( length(participantArray)  )
plot( y~xValues, axes=FALSE, ylab='% BOLD (contra-lateral)', xlab='cortical depth', las=1,
      ylim=c( min(y-y_se)-0.1*min(y-y_se), max(y+y_se)+0.1*max(y+y_se) ), cex.lab=1.25, pch=15  )
segments( xValues, y-y_se, xValues, y+y_se )
axis( 1, seq(0,1,0.5), cex.axis=1.25 )
axis( 2, seq( round( min(y), 1 ), round( max(y), 0 ), length.out=3 ), cex.axis=1.25 )

y <- apply( depthArray[,,2], 2, mean )
y_se <- ( apply( depthArray[,,2], 2, sd ) / sqrt( length(participantArray)  ) ) 
plot( y~xValues, axes=FALSE, ylab='% BOLD (ipsi-lateral)', xlab='cortical depth', las=1,
      ylim=c( -1, 0 ), cex.lab=1.25, pch=15 )
segments( xValues, y-y_se, xValues, y+y_se )
axis( 1, seq(0,1,0.5), cex.axis=1.25 )
axis( 2, seq( -1, 0, length.out=3 ), cex.axis=1.25 )

# y <- apply( depthArray[,,3], 2, mean )
# y_se <- apply( depthArray[,,3], 2, sd ) / sqrt( length(participantArray)  )
# plot( y~xValues, axes=FALSE, ylab='T stat (ipsi-lateral)', xlab='cortical depth', las=1,
#       ylim=c( 0, 1.5 ), cex.lab=1.25, pch=15 )
# segments( xValues, y-y_se, xValues, y+y_se )
# axis( 1, seq(0,1,0.5), cex.axis=1.25 )
# axis( 2, seq( 0, 1.5, length.out=3 ), cex.axis=1.25 )

dev.copy2pdf(file = sprintf('%s/%s.pdf', mainDirSavePlots, 'average_acrossDepth_percent_BOLD_all' ), height=3.5, width=6)
dev.off()



x11(height=3.5, width=6)
par(mfrow=c(1,2))
xValues <- seq(0,1,length.out=length(xVar))
colArray <- c('orange','blue','lightblue','red')
dataPlot <- t( apply( depthArrayAngle[,,1,], c(3,2), mean ) )
dataPlot_sd <- t( apply( depthArrayAngle[,,1,], c(3,2), sd ) ) / sqrt( length( participantArray ) )
seq(0,180,45)
matplot( xValues, dataPlot, axes=FALSE, ylab='% BOLD (contra-lateral)',
         xlab='cortical depth', las=1, cex.lab=1.25, pch=15, ylim=c(2,7), col=colArray  )
segments( xValues, dataPlot[,1]-dataPlot_sd[,1], xValues, dataPlot[,1]+dataPlot_sd[,1],
          col=colArray[1] )
segments( xValues, dataPlot[,2]-dataPlot_sd[,2], xValues, dataPlot[,2]+dataPlot_sd[,2],
          col=colArray[2] )
segments( xValues, dataPlot[,3]-dataPlot_sd[,3], xValues, dataPlot[,3]+dataPlot_sd[,3],
          col=colArray[3] )
segments( xValues, dataPlot[,4]-dataPlot_sd[,4], xValues, dataPlot[,4]+dataPlot_sd[,4],
          col=colArray[4] )
legend("bottomright", legend = c('(0,45]','(45,90]','(90,135]','(135,180]'), col=colArray, pch=15, bty='n')
axis( 1, seq(0,1,0.5), cex.axis=1.25 )
axis( 2, seq( 2, 7, length.out=3 ), cex.axis=1.25 )

dataPlot <- t( apply( depthArrayAngle[,,2,], c(3,2), mean ) )
dataPlot_sd <- t( apply( depthArrayAngle[,,2,], c(3,2), sd ) ) / sqrt( length( participantArray ) )
matplot( xValues, dataPlot, axes=FALSE, ylab='% BOLD (contra-lateral)',
         xlab='cortical depth', las=1, cex.lab=1.25, pch=15, ylim=c(-1.1,0), col=colArray  )
segments( xValues, dataPlot[,1]-dataPlot_sd[,1], xValues, dataPlot[,1]+dataPlot_sd[,1],
          col=colArray[1] )
segments( xValues, dataPlot[,2]-dataPlot_sd[,2], xValues, dataPlot[,2]+dataPlot_sd[,2],
          col=colArray[2] )
segments( xValues, dataPlot[,3]-dataPlot_sd[,3], xValues, dataPlot[,3]+dataPlot_sd[,3],
          col=colArray[3] )
segments( xValues, dataPlot[,4]-dataPlot_sd[,4], xValues, dataPlot[,4]+dataPlot_sd[,4],
          col=colArray[4] )
axis( 1, seq(0,1,0.5), cex.axis=1.25 )
axis( 2, seq( -1.1, 0, length.out=3 ), cex.axis=1.25 )

dev.copy2pdf(file = sprintf('%s/%s.pdf', mainDirSavePlots, 'average_acrossDepth_angle_percent_BOLD_all' ), height=3.5, width=6)
dev.off()


# dataPlot <- t( apply( depthArrayAngle[,,3,], c(3,2), mean ) )
# dataPlot_sd <- t( apply( depthArrayAngle[,,3,], c(3,2), sd ) ) / sqrt( length( participantArray ) )
# matplot( xValues, dataPlot, axes=FALSE, ylab='% BOLD (contra-lateral)',
#          xlab='cortical depth', las=1, cex.lab=1.25, pch=15, ylim=c(0,2), col=colArray  )
# segments( xValues, dataPlot[,1]-dataPlot_sd[,1], xValues, dataPlot[,1]+dataPlot_sd[,1],
#           col=colArray[1] )
# segments( xValues, dataPlot[,2]-dataPlot_sd[,2], xValues, dataPlot[,2]+dataPlot_sd[,2],
#           col=colArray[2] )
# segments( xValues, dataPlot[,3]-dataPlot_sd[,3], xValues, dataPlot[,3]+dataPlot_sd[,3],
#           col=colArray[3] )
# segments( xValues, dataPlot[,4]-dataPlot_sd[,4], xValues, dataPlot[,4]+dataPlot_sd[,4],
#           col=colArray[4] )
# axis( 1, seq(0,1,0.5), cex.axis=1.25 )
# axis( 2, seq( 0, 2, length.out=3 ), cex.axis=1.25 )
 

posBoldMatrix <- rbind( storeBold[,,1,1], storeBold[,,1,2], storeBold[,,1,3] )
filenameMat <- sprintf('%s/%s.txt',mainDirSavePlots,'posBold')
write.table( posBoldMatrix, file=filenameMat, row.names=FALSE, col.names=FALSE )

negBoldMatrix <- rbind( storeBold[,,2,1], storeBold[,,2,2], storeBold[,,2,3] )
filenameMat <- sprintf('%s/%s.txt',mainDirSavePlots,'negBold')
write.table( negBoldMatrix, file=filenameMat, row.names=FALSE, col.names=FALSE )

plot( apply( negBoldMatrix[13:18,], 2, mean ) )

