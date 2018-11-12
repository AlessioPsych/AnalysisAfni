maskFile = BrikLoad('motionCorrect.results/brainMask+orig.BRIK'); % brain mask, 

phaseNiiFiles = dir( strcat( phaseDir, '/*.nii' ) );

for nFile = 1:length( phaseNiiFiles )
    
    [ phaseSeries infoPhase ] = BrikLoad( [ phaseDir, '/', phaseNiiFiles(nFile).name ] );
    
    mask = logical( single( maskFile ) );
    
    phaseVolume = scaleif( single( phaseSeries(:,:,:,1) ), -pi, pi );
    plot_axialSagittalCoronal( phaseVolume, 1, [-pi pi], 'Masked, wrapped phase' )
    plot_axialSagittalCoronal( mask, 2, [0 1], 'Mask' )
    
    unwrapFlag = 1; 
    % 1 martinos code, 2 utrecht code; there are two versions    
    % of the phase unwrapping algorithm
    % one is developed in Utrecht, the second one at Martinos
    % center, the latter works betten in my opinion,
    % you can find the details of Martinos code here:
    %
    % https://www.martinos.org/~berkin/publications.html
    %
    % B. Bilgic, A.P. Fan, J.R. Polimeni, S.F. Cauley,
    % M. Bianciardi, E. Adalsteinsson, L.L. Wald, K. Setsompop;
    % Fast Quantitative Susceptibility Mapping with
    % L1-Regularization and Automatic Parameter Selection;
    % Magnetic Resonance in Medicine
    
    performRecursiveFilt = 1; % performs SHARP filtering step (described in the publication above)
    
    dimSelected = size( phaseSeries );
    volSeries = single( NaN( dimSelected ) );
    for k = 1:size(volSeries,4) % loop across phase volumes
        phaseVolume = single( phaseSeries(:,:,:,k) );
        phaseScaled = scaleif( phaseVolume, -pi, pi ) ; % scale between -pi and pi
        [unwrapVolume laplacianFilter residual] = unwrap_filter( phaseScaled, mask, unwrapFlag, 0 ,performRecursiveFilt);  % perform phase unwrapping
        volSeries(:,:,:,k) = unwrapVolume; % store
    end
    
    Opt.Prefix = [ phaseDir,  '/', phaseNiiFiles(nFile).name(1:end-4), '_unwrap' ];
    WriteBrik( volSeries, infoPhase, Opt )
    

    
end

eval( sprintf(' cd %s/', phaseDir ) )
phaseUnwrappedFiles = dir('*+orig.BRIK');
for k = 1:length( phaseUnwrappedFiles )
    convertToNiftiCommand = sprintf( '3dAFNItoNIFTI %s', phaseUnwrappedFiles(k).name );
    system( convertToNiftiCommand )
end

cd ..

% load('motionTransformationMatrices.mat') % load motion transformation matrices derived from motion correction
%                                          % applied to amplitude volumes. I derived these matrices from
%                                          % a separate mrVISTA session loading amplitude images                                        
%                                          
% %% motion correct phase unwrapped images
% warpedVolSeries = single( NaN( dimSelected ) );
% scan = 1;
% wMotion = MMwith; % within scan motion
% bMotion = MMbw; % between scans motion
% nFrames = size(volSeries,4);
% interpMethod = '*linear';
% 
% waitHandle = waitbar(0,'Motion correction...');
% for volume = 1:nFrames % loop across phase volumes and apply motion correction
%     waitbar(volume/nFrames);
%     if (scan == 1) && (volume==1)        
%         warpedVolSeries(:,:,:,volume) = volSeries(:,:,:,volume);
%     elseif (scan == 1) && (volume~=1)
%         tMat = wMotion{scan}(:,:,volume);
%         warpedVolSeries(:,:,:,volume) = warpAffine3(volSeries(:,:,:,volume),tMat,NaN,1,interpMethod);
%     elseif (scan~=1) && (volume==1)
%         tMat = bMotion(:,:,scan);
%         warpedVolSeries(:,:,:,volume) = warpAffine3(volSeries(:,:,:,volume),tMat,NaN,1,interpMethod);
%     else
%         tMat = bMotion(:,:,scan) * wMotion{scan}(:,:,volume);
%         warpedVolSeries(:,:,:,volume) = warpAffine3(volSeries(:,:,:,volume),tMat,NaN,1,interpMethod);
%     end
% end
% close(waitHandle)
% 
% %% save motion corrected phase unwrapped images
% maskFile.data = mean( warpedVolSeries, 4 ); % compute average between phase unwrapped volumes
% maskFile.fname = 'phaseUnwrapped.nii.gz';
% writeFileNifti(maskFile) % save average phase unwrapped volume
% 
% %% motion correct amplitude images
% warpedVolSeries = single( NaN( dimSelected ) );
% scan = 1;
% wMotion = MMwith;
% bMotion = MMbw;
% nFrames = size(volSeries,4);
% interpMethod = '*linear';
% volSeries = amplitudeSeries.data;
% 
% waitHandle = waitbar(0,'Motion correction...');
% for volume = 1:nFrames
%     waitbar(volume/nFrames); % loop across amplitude volumes and apply motion correction
%                              % this is kind of redundant.
%     if (scan == 1) && (volume==1)        
%         warpedVolSeries(:,:,:,volume) = volSeries(:,:,:,volume);
%     elseif (scan == 1) && (volume~=1)
%         tMat = wMotion{scan}(:,:,volume);
%         warpedVolSeries(:,:,:,volume) = warpAffine3(volSeries(:,:,:,volume),tMat,NaN,1,interpMethod);
%     elseif (scan~=1) && (volume==1)
%         tMat = bMotion(:,:,scan);
%         warpedVolSeries(:,:,:,volume) = warpAffine3(volSeries(:,:,:,volume),tMat,NaN,1,interpMethod);
%     else
%         tMat = bMotion(:,:,scan) * wMotion{scan}(:,:,volume);
%         warpedVolSeries(:,:,:,volume) = warpAffine3(volSeries(:,:,:,volume),tMat,NaN,1,interpMethod);
%     end
% end
% close(waitHandle)
% 
% %% save motion corrected phase unwrapped images
% maskFile.data = mean( warpedVolSeries, 4 ); % compute average between amplitude volumes
% maskFile.fname = 'amplitude.nii.gz';
% writeFileNifti(maskFile) % save average amplitude volume
% 
% 
% % to visualize the results of the procedure I use 3dSlicer
% % because it allows to superimpose 2 or more volumes and 
% % adjust the transparency between the two as well as
% % contrast and brightness. Examples of the visualization
% % are in the files 'amplitude.png', 'phase.png' and
% % 'amplitude_phase.png' (the combined volume)
% 
% 
