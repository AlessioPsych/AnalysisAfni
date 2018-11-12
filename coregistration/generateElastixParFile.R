args <- commandArgs(T)
print( args )

setwd('/home/alessiofracasso/data/elastix_example_v4.8')

args <- c( 'FixedInternalImagePixelType', 'float',
           'MovingInternalImagePixelType', 'float', 
           'UseDirectionCosines', 'true', 
           
           'ResampleInterpolator', 'FinalBSplineInterpolator',
           'FixedImagePyramid', 'FixedRecursiveImagePyramid',
           'MovingImagePyramid', 'MovingRecursiveImagePyramid',
           'Registration', 'MultiResolutionRegistration',
           'Interpolator', 'BSplineInterpolator',
           'Metric', 'AdvancedMattesMutualInformation',
           'BSplineInterpolationOrder', '3',

           'Resampler', 'DefaultResampler',
           'Optimizer', 'AdaptiveStochasticGradientDescent',
           'Transform', 'EulerTransform',
                      
           'ErodeMask', 'false',
           'ErodeFixedMask', 'false',
           
           'NumberOfResolutions', '4',
           'MaximumNumberOfIterations', '1500',
           'AutomaticScalesEstimation', 'true',
           'AutomaticTransformInitialization', 'true',
           
           'HowToCombineTransforms', 'Compose',
           
           'NumberOfHistogramBins', '32',
           
           'NumberOfSpatialSamples', '10000',
           'ImageSampler', 'RandomCoordinate',
           'CheckNumberOfSamples', 'false',
           'NewSamplesEveryIteration', 'true',
           'MaximumNumberOfSamplingAttempts', '100',
           'FinalBSplineInterpolationOrder', '1',
           
           'DefaultPixelValue', '0',
           'WriteTransformParametersEachIteration', 'false',
           'WriteResultImage', 'true',
           'ResultImageFormat', 'nii',
           'ResultImagePixelType', 'float' )

## add the parameters for bspline:
#(FinalGridSpacingInPhysicalUnits 35 35 35)
#(MovingImageDerivativeScales 1 1 1)


for (k in seq(1,(length(args)-1),2) ) {
  argSplit01 <- strsplit( args[k], split='-' )
  argSplit02 <- strsplit( args[k+1], split='-' )
  arg1 <- argSplit01[[1]]
  arg2 <- argSplit02[[1]]
  
  ## if the argument is numeric, then 
  nArg2 <- as.numeric( arg2 )
  if (is.na( nArg2 )[1])  {
    stringOut <- sprintf( '(%s %s)', arg1, arg2  )
  }
  if (!is.na( nArg2 )[1]) {
    for ( nEl in 1:length( nArg2 ) ) {
      if (nEl==1) { stringOut <- sprintf( '%s %d', arg1, nArg2[ nEl ]  ) }
      if (nEl>1) { stringOut <- sprintf( '%s %d', stringOut, nArg2[ nEl ]  ) }
    }
    stringOut <- sprintf( '(%s)', stringOut )
  }
    
  print( stringOut )
  
  if (k==1) { stringOutFull <- c( stringOut ) }
  if (k>1) { stringOutFull <- c( stringOutFull, stringOut ) }
  
}

fileConn<-file("output.txt")
writeLines( stringOutFull, fileConn)
close(fileConn)
