debug <- 0

if (debug==0) {
  args <- commandArgs(T)
  print( args )
}
if (debug==1) {
  setwd('/analyse/Project0165/Fracasso_Anatomy/tests/testMRDenoisePhase/dirTest')
  args <- c('EPIFOLDER/','PHASEFOLDER/','*.nii','4','EPI_denoised_test/')
}

mainDir <- getwd()

set.seed(1);
afniInstallDir <- Sys.getenv(x='AFNI_INSTALLDIR')
afniAtlasDir <- Sys.getenv(x='AFNI_ATLASDIR')
afniSurfaces <- Sys.getenv(x='AFNI_TOOLBOXDIRSURFACES')
generalPurposeDir <- Sys.getenv(x='AFNI_TOOLBOXDIRGENERALPURPOSE')
source( sprintf( '%s/scaleData.R', generalPurposeDir ) )
source( sprintf( '%s/AFNIio.R', afniInstallDir ) )
source( sprintf( '%s/coordinateFromLinearIndex.r', afniSurfaces ) )
source( sprintf( '%s/linearIndexFromCoordinate.r', afniSurfaces ) )
library( pracma )

#arrange input
inputEpiDir <- args[1]
inputEpiPhaseDir <- args[2]
epiFormat <- args[3]
ncores <- args[4]
outputNameDir <- args[5]

print( sprintf('input EPI directory = %s', inputEpiDir ) )
print( sprintf('input EPI phase directory = %s', inputEpiPhaseDir ) )
print( sprintf('EPI file format = %s', epiFormat ) )
print( sprintf('ncores = %s', ncores ) )
print( sprintf('output directory = %s', outputNameDir ) )

# create output directory
instr <- sprintf('mkdir %s', outputNameDir ); system( instr )

# get epi filenames
setwd( inputEpiDir )
filesWorkDir <- dir(pattern=epiFormat);
setwd( mainDir )

# get epi phase filenames
setwd( inputEpiPhaseDir )
filesPhaseWorkDir <- dir(pattern=epiFormat);
setwd( mainDir )

# perform data denoising correction
for (k in 1:length(filesWorkDir)) { # k <- 1

  setwd( mainDir )
  
  # convert phase into radians, from -pi to pi
  print( sprintf('convert phase into radians, file %s', filesPhaseWorkDir[k] ) )
  phaseFile <- read.AFNI( sprintf( '%s/%s', inputEpiPhaseDir, filesPhaseWorkDir[k] ) )
  phaseBrk <- phaseFile$brk
  oldMax <- max( phaseBrk )
  oldMin <- min( phaseBrk )
  newMax <- round( pi, 4 )
  newMin <- round( -pi, 4 )
  delta <- (newMax-newMin)/(oldMax-oldMin);
  
#  instr <- sprintf('3dcalc -a %s%s -expr \u0027 %1.6f*(a-%1.2f)+%1.2f \u0027 -datum short -prefix %s/phaseScaled_pipi.nii ', 
#                 inputEpiPhaseDir, filesPhaseWorkDir[k], delta, oldMin, newMin, outputNameDir )
  instr <- sprintf('3dcalc -a %s%s -expr \u0027 %1.6f*(a-%1.2f)+%1.2f \u0027 -prefix %s/phaseScaled_pipi.nii ', 
                   inputEpiPhaseDir, filesPhaseWorkDir[k], delta, oldMin, newMin, outputNameDir )
  print( instr )
  system( instr )
  
  print( sprintf('convert phase files into mif format' ) )
  instr <- sprintf( 'mrconvert %s/phaseScaled_pipi.nii %s/phase_pipi.mif', outputNameDir, outputNameDir  )
  print( instr )
  system( instr )

  print( sprintf('convert amplitude into mif format' ) )
  instr <- sprintf( 'mrconvert %s/%s %s/magnitude_pipi.mif', inputEpiDir, filesWorkDir[k], outputNameDir  )
  print( instr )
  system( instr )

  print( sprintf('create complex file' ) )
  instr <- sprintf( 'mrcalc %s/magnitude_pipi.mif %s/phase_pipi.mif -polar %s/complex_pipi.mif', outputNameDir, outputNameDir, outputNameDir  )
  print( instr )
  system( instr )
  
  print( sprintf('processing complex epi file: complex_pipi.mif, merged from %s and %s', filesWorkDir[k], filesPhaseWorkDir[k] ) )
  
  instr <- sprintf( 'dwidenoise -nthreads %s %s/complex_pipi.mif %s/out_pipi.mif', ncores, outputNameDir, outputNameDir )
  print( instr )
  system( instr )
  
  instr <- sprintf( 'mrcalc %s/out_pipi.mif -abs %s/out_magnitude_pipi.mif', outputNameDir, outputNameDir )
  print( instr )
  system( instr )
  
  instr <- sprintf( 'mrconvert %s/out_magnitude_pipi.mif %s/%s', outputNameDir, outputNameDir, filesWorkDir[k] )
  print( instr )
  system( instr )
  
  print( sprintf('clean up') )
  setwd( outputNameDir )
  system('rm *.mif')
  system('rm phaseScaled_pipi.nii')
  
  #instr <- sprintf('dwidenoise -nthreads %s %s/%s%s %s/%s%s ', ncores, mainDir, inputEpiDir, filesWorkDir[k], mainDir, outputNameDir, filesWorkDir[k]  ); system(instr)
  # does it work with datum long?
  
}
