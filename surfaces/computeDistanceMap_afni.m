function [dist distCsf v f] = computeDistanceMap_afni( params )


% regular distance map

% Build an isodensity surface at the gray-white interface, and get
% its vertices:
volume = params.volume;
whiteVolume = (volume == 3 | volume == 4); 
csfVolume = volume == 1;

whiteVolume = smooth3(whiteVolume,'gaussian',params.gaussianSize,params.gaussianStd);
csfVolume = smooth3(csfVolume,'gaussian',params.gaussianSize,params.gaussianStd);

whiteVolume = smooth3(whiteVolume);
csfVolume = smooth3(csfVolume);

[f, v] = isosurface(whiteVolume, params.whiteThr);
[fCsf, vCsf] = isosurface(csfVolume, params.csfThr);
%[distPIndex distP] = nearpoints(v',vCsf');
%[distPIndex distP] = dsearchn( v, vCsf );


%[errFlag v] = AFNI_XYZcontinuous2Index( params.v_vm, params.info );
%f = params.f_vm;
%[errFlag vCsf] = AFNI_XYZcontinuous2Index( params.v_csf, params.info );
%fCsf = params.f_csf;
%n = isonormals(whiteVolume,v);

% white distance map
v1 = v(:, [2 1 3])';
%v1 = v';
whiteIso = v1;
thr = 0.5;

%[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(v,f);

% Calculate distances in the gray matter (from white matter, white matter
% surface with gray set to 0, increasing in gray and outside the brain
%whiteInds = find(~whiteVolume);
whiteInds = find( whiteVolume >= thr );
whiteVerts = indices2Coords(whiteInds,  size(whiteVolume) );
[inds, whiteDist] = nearpoints(whiteVerts, whiteIso);
oDist = repmat(NaN, size(whiteVolume));
%oDist(whiteInds) = sqrt(whiteDist).*anatomy.pixdim(1); % this volume is in the original space (244 58 244)
oDist(whiteInds) = sqrt(whiteDist).*1; % this volume is in the original space (244 58 244)

% Calculate distances in the gray matter (from white matter, white matter
% surface with gray set to 0, increasing in gray and outside the brain
%whiteInds = find(whiteVolume);
whiteInds = find(whiteVolume < thr );
whiteVerts = indices2Coords(whiteInds,  size(whiteVolume) );
[inds, whiteDist] = nearpoints(whiteVerts, whiteIso);
iDist = repmat(NaN, size(whiteVolume));
%iDist(whiteInds) = sqrt(whiteDist).*anatomy.pixdim(1); % this volume is in the original space (244 58 244)
iDist(whiteInds) = sqrt(whiteDist).*1; % this volume is in the original space (244 58 244)

dist = repmat(NaN, size(whiteVolume));
dist( whiteVolume >= thr ) = -1.*oDist( whiteVolume >= thr );
dist( whiteVolume < thr ) = iDist( whiteVolume < thr );

% csf distance map
v1 = vCsf(:, [2 1 3])';
%v1 = vCsf';
whiteIso = v1;
thr = 0.5;

% Calculate distances in the gray matter (from white matter, white matter
% surface with gray set to 0, increasing in gray and outside the brain
%whiteInds = find(~whiteVolume);
whiteInds = find( csfVolume >= thr );
whiteVerts = indices2Coords(whiteInds,  size(whiteVolume) );
[inds, whiteDist] = nearpoints(whiteVerts, whiteIso);
oDist = repmat(NaN, size(whiteVolume));
%oDist(whiteInds) = sqrt(whiteDist).*anatomy.pixdim(1); % this volume is in the original space (244 58 244)
oDist(whiteInds) = sqrt(whiteDist).*1; % this volume is in the original space (244 58 244)

% Calculate distances in the gray matter (from white matter, white matter
% surface with gray set to 0, increasing in gray and outside the brain
%whiteInds = find(whiteVolume);
whiteInds = find( csfVolume < thr );
whiteVerts = indices2Coords(whiteInds,  size(whiteVolume) );
[inds, whiteDist] = nearpoints(whiteVerts, whiteIso);
iDist = repmat(NaN, size(whiteVolume));
%iDist(whiteInds) = sqrt(whiteDist).*anatomy.pixdim(1); % this volume is in the original space (244 58 244)
iDist(whiteInds) = sqrt(whiteDist).*1; % this volume is in the original space (244 58 244)

distCsf = repmat(NaN, size(whiteVolume));
distCsf( csfVolume >= thr ) = -1.*oDist( csfVolume >= thr );
distCsf( csfVolume < thr ) = iDist( csfVolume < thr );

%save('distanceMap.mat','dist','distCsf','f','v','n','distP','distPIndex','fCsf','vCsf','Cmean');
