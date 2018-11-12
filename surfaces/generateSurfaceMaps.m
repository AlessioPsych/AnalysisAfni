function generateSurfaceMaps( params )

[ functionalVolume infoVol ] = BrikLoad( params.inVol );
%[ distVolume ] = BrikLoad( params.distVol );
nSteps = length( params.depth );
mapValues = zeros( [ size( params.v_wm, 1 ) ...
    length( params.depth ) size( functionalVolume, 4 ) ] );



%[errFlag v] = AFNI_XYZcontinuous2Index( params.v_wm, params.info );
%f = params.f_wm;

convMat = inv( [ reshape( params.info.IJK_TO_DICOM_REAL, [4 3] )'; 0 0 0 1 ] );
v = convMat * [ params.v_wm ones( size(params.v_wm, 1), 1) ]';
v = v(1:3,:)';
v = v(:,[2 1 3]);
f = params.f_wm;

for k=1:nSteps
    
    step = params.depth( k );
    new_v_wm = perform_normal_displacement( v, f, step );
    for j=1:size( functionalVolume, 4 )
        mapValues(:,k,j) = interp3( functionalVolume(:,:,:,j), new_v_wm(:,1), ...
            new_v_wm(:,2), new_v_wm(:,3), 'linear', 0 );
    end
    
    progressbar(k,nSteps)

end

for k = 1:size( functionalVolume, 4 )
    filename = strcat( params.saveDir, 'surfaceMap_', cell2mat( params.names(k) ), '.1D' );
    dlmwrite( filename, [ mapValues(:,:,k) mapValues(:,:,1) ], 'delimiter', '\t', 'precision', 4 );
end