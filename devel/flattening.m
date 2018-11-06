function flattening

coordFile = getenv('COORDNAME');
topoFile = getenv('TOPONAME');
pathDir = getenv('NUMTOURSDIR');
currDir = getenv('CURRDIR');

coordPath = [ currDir, '/', coordFile];
topoPath = [ currDir, '/', topoFile];
disp( coordPath );
disp( topoPath );
disp( pathDir )
disp( currDir )

addpath( genpath( [ pathDir , '/toolbox_general' ] ) )
addpath( genpath( [ pathDir , '/toolbox_graph' ] ) )

roiCoords = dlmread( coordPath );
roiTopo = dlmread( topoPath );
roiTopo = roiTopo+1;
opt1.face_vertex_color = repmat( linspace( 0, 1, size(roiCoords,1) ), 3, 1 )';
opt1.name = 'iii';
%figure(1)
%plot_mesh( roiCoords', roiTopo', opt1 )
%shading faceted
vertex = roiCoords';
faces = roiTopo';
n = size(vertex,2);

% Cotan weights
options.symmetrize = 1;
options.normalize = 0;
L = compute_mesh_laplacian(vertex,faces,'conformal',options);
% boundary
options.verb = 0;
boundary = compute_boundary(faces, options);
% fixed positions
p = length(boundary);
t = linspace(0,2*pi,p+1)'; t(p) = [];
x0 = cos(t); y0 = sin(t);
% system
L1 = L;
L1(boundary,:) = 0;
L1(boundary + (boundary-1)*n) = 1;
% Set up the right hand sizes with the fixed position.
Rx = zeros(n,1); Rx(boundary) = x0;
Ry = zeros(n,1); Ry(boundary) = y0;
% solve
x = L1 \ Rx;
y = L1 \ Ry;
vertexF = [x';y'];

vertexOut = [ vertexF; ones(1,n) ];
topoOut = roiTopo';

cd( currDir )

outVertexFile = split( coordFile, '_' );
outVertexFile = sprintf( '%s_coords_flat.1D', outVertexFile{1} );
outFacesFile = split( topoFile, '_' );
outFacesFile = sprintf( '%s_topo_flat.1D', outFacesFile{1} );
dlmwrite(outVertexFile, vertexOut, 'delimiter','\t','precision',3 );
dlmwrite(outFacesFile, topoOut, 'delimiter','\t','precision',3 );

% align
%vertexF = vertexF - repmat(vertexF(:,icenter), [1 n]);
%theta = -pi/2+atan2(vertexF(2,irotate),vertexF(1,irotate));
%vertexF = [vertexF(1,:)*cos(theta)+vertexF(2,:)*sin(theta); ...
%           -vertexF(1,:)*sin(theta)+vertexF(2,:)*cos(theta)];
% display
% clf;
% figure(2)
% plot_mesh(vertexF,faces,opt1);

% figure(3)
% plot_mesh( [ vertexF; ones(1,n) ], roiTopo', opt1 )
% shading faceted


