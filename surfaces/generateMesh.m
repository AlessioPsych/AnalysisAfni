cd anatomy

system('rm *.spec')
system('rm *.dset')
system('rm *.coord')
system('rm *.topo')
system('rm t1MyelinSegmentation_bool_WM+orig.BRIK')
system('rm t1MyelinSegmentation_bool_WM+orig.HEAD')
system('rm t1MyelinSegmentation_bool_CSF+orig.BRIK')
system('rm t1MyelinSegmentation_bool_CSF+orig.HEAD')
system('rm distance_WM+orig.BRIK')
system('rm distance_WM+orig.HEAD')
system('rm *.dset')
system('rm *.M2M')
system('rm *.1D')

PATH = getenv('PATH');
setenv( 'PATH', [PATH pathToAFNICommands] );
setenv( 'PATH', [PATH pathToSHCommands] );

myStartup('+','matFileAddOn')
myStartup('+','afni_matlab');
myStartup('+','vistasoft')

growNormalFlag = 0;
methodMeshSmooth = 1;
smoothIterationsWm = 150;
smoothIterations = 3;

anatomyFile = anatomyFileInput;
segmentationFile = segmentationFileInput;
[volAnatomy infoAnatomy] = BrikLoad( anatomyFile );
[segmentation info] = BrikLoad( segmentationFile );
segmentationBool = double( zeros( size( segmentation ) ) );
segmentationBool( segmentation == 3 | segmentation == 4 ) = 1;
Opt.Prefix = 't1MyelinSegmentation_bool_WM';
WriteBrik( segmentationBool, infoAnatomy, Opt )

[segmentation info] = BrikLoad( segmentationFile );
segmentationBool = double( zeros( size( segmentation ) ) );
segmentationBool( segmentation == 3 | segmentation == 4 |...
    segmentation == 5 | segmentation == 6 ) = 1;
Opt.Prefix = 't1MyelinSegmentation_bool_CSF';
WriteBrik( segmentationBool, infoAnatomy, Opt )

system('. surfaceCommands_1')
 
%% load surfaces
v_wm = dlmread( 'WM_Surface_smooth_consistent.1D.coord' );
f_wm = dlmread( 'WM_Surface_smooth_consistent.1D.topo' );
f_wm = f_wm + 1;
normals_white = compute_normal( v_wm, f_wm )';

[v_wm_s, f_wm_s] = smoothMesh( v_wm, f_wm, smoothIterationsWm);
dlmwrite( 'WM_Surface_inflated.1D.coord', v_wm_s, 'delimiter', ' ' );
dlmwrite( 'WM_Surface_inflated.1D.topo', f_wm_s-1, 'delimiter', ' ' );

system('. surfaceCommands_3')
%system('. ../computeSurfaceDistance')

v_csf = dlmread( 'CSF_Surface_smooth_consistent.1D.coord' );
f_csf = dlmread( 'CSF_Surface_smooth_consistent.1D.topo' );
f_csf = f_csf + 1;

[segmentation info] = BrikLoad( segmentationFile );
params.volume = segmentation;
params.info = info;
params.v_wm = v_wm;
params.f_wm = f_wm;
params.v_csf = v_csf;
params.f_csf = f_csf;
params.gaussianSize = 3;
params.gaussianStd = 0.6; 
params.whiteThr = 0.5;
params.csfThr = 0.5;
[dist distCsf v f] = computeDistanceMap_afni( params );
Opt.Prefix = 'distance_WM';
WriteBrik( dist, infoAnatomy, Opt )

cd ..




% clean up generateMesh file
% build a file that generates maps based on
% distance from gm/wm border
% work on rois

% if growNormalFlag == 1
%     
%     %% grow normals, dompute distance with csf mesh and take the first minima of the grown normal
%     dDisp = 0.20;
%     startSampling = -0.20;
%     displacement = startSampling;
%     nSteps = 20;
%     storeDist = NaN( size( v_wm, 1 ), nSteps );
%     storeCoord = NaN( size( v_wm, 1 ), 3, nSteps );
%     storeCsfVertex = NaN( size( v_wm, 1 ), nSteps );
%     
%     for n = 1:nSteps
%         
%         displacement = displacement + dDisp;
%         vertex1 = perform_normal_displacement( v_wm, f_wm, displacement );
%         
%         [vertex_id d] = dsearchn( v_csf, vertex1 );
%         
%         storeDist(:,n) = d;
%         storeCoord(:,:,n) = vertex1;
%         storeCsfVertex(:,n) = vertex_id;
%         
%         progressbar(n,nSteps);
%         
%     end
%     save('growNormals.mat');
% else
%     load('growNormals.mat');
% end
% 
% corticalThickness = zeros( size( storeDist, 1 ), 1 );
% counter = 1;
% flagLoop = 1;
% storeEndPoint = zeros( size( storeDist, 1 ), 3 );
% for k = 1:size(storeDist,1)
%     
%     distanceArray = storeDist(k,:);
%     [ sorted iS ] = sort( distanceArray );
%     indexProfile = iS(1);
%         
%     corticalThickness(k) = ( indexProfile(1).*dDisp ) + startSampling;
%         
%     storeEndPoint(k,:) = storeCoord( k, :, indexProfile(1) );
%         
%     progressbar( k, size(storeDist,1) );
% end
% hist(corticalThickness)
% dlmwrite( 'thickness.1D.dset', corticalThickness, 'delimiter', ' ' );
% 
% %% surface at middle cortical thickness:
% middleCorticalThickness = corticalThickness./2;
% v_middleCorticalThickness = zeros( size(v_wm) );
% for k = 1:size(v_wm,1)
%     vertexMiddle = v_wm(k,:) + middleCorticalThickness(k).*normals_white(k,:);  
%     v_middleCorticalThickness(k,:) =  vertexMiddle;
% end
% [v_middleCorticalThickness, faces] = smoothMesh( v_middleCorticalThickness, f_wm, smoothIterations);
% dlmwrite( 'middle_cortical_thickness.1D.coord', v_middleCorticalThickness, 'delimiter', ' ' );
% 
% 
% %% distance map;
% switch methodMeshSmooth
%     case {1}
%         [v2, faces] = smoothMesh( storeEndPoint, f_wm, smoothIterations);
%     case {2}
%         % parameters for the operator
%         laplacian_type = 'distance';
%         options.symmetrize = 0;
%         options.normalize = 1; % it must be normalized for filtering
%         options.verb = 0;
%         W = compute_mesh_weight(storeEndPoint,f_wm,laplacian_type,options);
%         % This is the corresponding laplacian
%         L = compute_mesh_laplacian(storeEndPoint,f_wm,laplacian_type,options);
%         
%         v2 = storeEndPoint';
%         clf;
%         options.face_vertex_color = [];
%         for i=1:smoothIterations
%             v2 = (W*(W*v2'))';
%         end
%         v2 = v2';  
% end
% 
% dlmwrite( 'csf_distance.1D.coord', v2, 'delimiter', ' ' );
% 
% system(' . ../surfaceCommands_2 ')
% 
% cd ..


