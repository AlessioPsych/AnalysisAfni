%% load roi
%roiFilename = strcat('anatomy/',roiName);
%distanceMapFile = 'anatomy/distance_WM+orig.BRIK';
%coordsMeshFile = 'anatomy/WM_Surface_smooth_consistent.1D.coord';
%faceMeshFile = 'anatomy/WM_Surface_smooth_consistent.1D.topo';
roiName = 'leftV1.niml.roi';
roiFilename = strcat('anatomy_CBS/',roiName);
segFile = 'anatomy_CBS/SEG_al_epi.nii.gz';
distanceMapFile = 'anatomy_CBS/DEPTH_al_epi.nii.gz';
coordsMeshFile = 'anatomy_CBS/boundary00_sm.1D.coord';
faceMeshFile = 'anatomy_CBS/boundary00_or.1D.topo';



[ roiFileParts1 roiFileParts2 ] = fileparts( roiFilename );
v_wm = dlmread( coordsMeshFile );
f_wm = dlmread( faceMeshFile );
f_wm = f_wm + 1;

S = afni_niml_readsimpleroi( roiFilename );
emptyNodesFilt = [ S{1}.region{1}; S{1}.edge{1} ];
options.face_vertex_color = zeros( size(v_wm,1), 3 );
options.face_vertex_color( emptyNodesFilt, : ) = repmat( [ 250 0 0 ], length( emptyNodesFilt ), 1 );
figure(10), plot_mesh( v_wm, f_wm, options )

[ errFlag segVol infoVolume ] = BrikLoad( segFile );

nSteps = 8;
storeVertex = zeros( [ size( v_wm ), nSteps ] );
storeIndxRoi = [];
stepArray=linspace(-0.5,3.5,nSteps);
for k=1:nSteps
    new_v_wm = perform_normal_displacement( v_wm, f_wm, stepArray(k) );
    storeVertex(:,:,k) = new_v_wm;
    vertexList = new_v_wm( emptyNodesFilt, : );
    [err, coordsRoi] = AFNI_XYZcontinuous2Index( vertexList, infoVolume );
    coordsRoiUnique = unique( coordsRoi, 'rows' );
    indexNull = find( coordsRoiUnique==0 );
    if ~isempty( indexNull )
        coordsRoiUnique( indexNull ) = 1;
    end
    indxRoi = coords2Indices( coordsRoiUnique', size( segVol ) );
    indexGray = segVol( indxRoi );
    storeIndxRoi = [storeIndxRoi indxRoi(indexGray==1)]; 
end

emptyVolume = zeros( size( segVol ) );
emptyVolume( storeIndxRoi ) = 1;
emptyVolume = smooth3( emptyVolume );
emptyVolume = emptyVolume >= 0.4;
Opt.Prefix = [roiFileParts1, '/', roiFileParts2];
Opt.OverWrite = 'y';
WriteBrik( emptyVolume, infoVolume, Opt )
system( sprintf( '3dAFNItoNIFTI -prefix %s  %s+orig.BRIK', Opt.Prefix(1:end-5), Opt.Prefix ) );


%[err, coordsRoi] = AFNI_XYZcontinuous2Index( v_wm( emptyNodesFilt, : ), infoVolume );
%coordsRoiUnique = unique( coordsRoi, 'rows' );
%indxRoi = coords2Indices( coordsRoiUnique', size( distanceMap ) );

%indxDistanceMap = find( distanceMap > -0.5 & distanceMap < 3 );
%nearpoints()

%roiCoordsRAS = [ v_wm( emptyNodesFilt, : ), ones( length(emptyNodesFilt), 1 ) ];
%tMat = [ reshape( infoVolume.IJK_TO_DICOM_REAL , [4 3] )'; 0 0 0 1 ] ;
%roiCoordsVOX = inv(tMat)*roiCoordsRAS';
%roiCoordsVOX = roiCoordsVOX(1:3,:);

%round(roiCoordsVOX)

%index_GM_WM = find( distanceMap > -0.5 & distanceMap < 0.5 );
%coords_GM_WM = indices2Coords( index_GM_WM, size( distanceMap ) );

%nearpoints( roiCoordsVOX )


%disp = 0.25;
%startGrow = -0.5;
%stopGrow = 3;
%arrayGrow = linspace( startGrow, stopGrow, disp );


