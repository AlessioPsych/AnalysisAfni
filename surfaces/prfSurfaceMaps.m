params.inVol = prfModelFile;
params.distVol = 'anatomy/distance_WM+orig.BRIK';
params.depth = linspace(0,2.5,6);
params.saveDir = 'anatomy/';
params.names{1} = 'coherence';
params.names{2} = 'sigmaUpdatedCenter';
params.names{3} = 'sigmaUpdatedSurround';
params.names{4} = 'sigma';
params.names{5} = 'eccentricity';
params.names{6} = 'betaInterceptArray';
params.names{7} = 'betaSlopeArray1';
params.names{8} = 'betaSlopeArray2';
generateSurfaceMaps( params )

cd anatomy

system('rm prfCoherence+orig.BRIK')
system('rm prfCoherence+orig.HEAD')
system('rm prfCoherence.nii')

[a b] = fileparts( params.inVol );
[ functionalVolume infoVol ] = BrikLoad( b );
prfCoherence = functionalVolume(:,:,:,1);
Opt.Prefix = strcat( 'prfCoherence' );
WriteBrik( prfCoherence, infoVol, Opt )
system('3dAFNItoNIFTI prfCoherence+orig');

cd ..