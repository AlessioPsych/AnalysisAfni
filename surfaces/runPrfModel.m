clear all
close all

dbstop if error

myStartup('+','matFileAddOn')
myStartup('+','afni_matlab')
myStartup('+','vistasoft')

PATH = getenv('PATH');
setenv( 'PATH', [PATH ':/usr/lib/afni/bin'] );

cd prfModel

%% load data and brain mask
[ts infoTs] = BrikLoad( 'meanTs.nii' );
[brainMask infoBrainMask] = BrikLoad( 'brainMask+orig' );
%[brainMask infoBrainMask] = BrikLoad( 'COPY_brainMask+orig.BRIK' );
ts = single( ts );

ts = ts(:,:,:,2:73);

brainMask = single( brainMask );

inplaneRoiIndex = find( brainMask == 1 );
inplaneRoiCoords = indices2Coords( inplaneRoiIndex,  size( brainMask ) );

tsR = single( zeros( size(ts,4), length(inplaneRoiIndex) ) );
for k = 1 : length(inplaneRoiIndex)
    tTemp = squeeze( ts( inplaneRoiCoords(1,k), ...
        inplaneRoiCoords(2,k), ...
        inplaneRoiCoords(3,k), : ) );
    tTemp = single( tTemp );
    tTemp = ( tTemp./mean(tTemp) ).*100;
    tTempIndex = find(tTemp>200);
    tTemp( tTempIndex ) = 200;
    tsR( :, k ) = tTemp;
    progressbar( k, length( inplaneRoiIndex ) );
end
tsSize = size( ts );
%clear ts;
%clear infoTs;

%% build stimuli
params.upFactor = 2;
params.physRes = [1024 768].*params.upFactor;
params.res = [1024 768].*params.upFactor;
params.screenSize = [10 10]; %cm
params.vDist = 30; %cm
params.degperpix=2*((atan(params.screenSize./(2*params.vDist))).*(180/pi))./params.res; %visual angle for 1 pixel
params.pixperdeg = 1./params.degperpix; %n of pixels for 1 visual angle;
params.pixperdeg = round(median(params.pixperdeg));
params.nSteps = 12;
params.stepSize = 0.5; % degree visual angle
params.stepSizeVoxels = params.stepSize*params.pixperdeg; % degree visual angle
params.ringsSuperimp = 0.25;
params.ringsSuperimpVoxels = params.ringsSuperimp*params.pixperdeg;
params.startCycle = 0.3;
params.startCycleVoxels = params.startCycle*params.pixperdeg;
params.sequenceType = 1;
params.blurStimuli = 0;
params.showStimuli = 0;
params.tr = 4;

stimuli = generateStimuli(params);

%% generate hrfs parameters

p1Prf = 6;%linspace( 2, 8, 2 );
p2Prf = 16;%linspace( 12, 20, 2 );
p3Prf = 1;%1;
p4Prf = 1;%linspace( 0.5, 1.5, 2 );
p5Prf = 6;%linspace( 4, 8, 3 );
p6Prf = 0;%linspace( 0, 1, 2 );

hrfMatrixAll = zeros( 7, length(p1Prf)*length(p2Prf)*length(p4Prf)*length(p5Prf)*length(p6Prf) );
count = 1;
for k = 1 : length( p1Prf )
    for a = 1 : length( p2Prf )
        for b = 1:length(p4Prf)
            for f = 1 : length( p5Prf )
                for l = 1 : length( p6Prf )
                    hrfPar = [ p1Prf(k) p2Prf(a) 1 p4Prf(b)  p5Prf(f) p6Prf(l) 32];
                    hrfMatrixAll( :, count ) = hrfPar;
                    count = count + 1;
                end
            end
        end
    end
end

tic

%% fit
multNegative = linspace( 1.1, 3, 40);
nSamples = 50;
nElementsLoop = 2;

countFitLoop = 1;
sigmaPrfArray = linspace(0.015,7,nSamples);
eccPrfArray = linspace(0.015,8,nSamples);
indexLoop = [ 1:nElementsLoop:nSamples+1 ];
sOut = [];
eOut = [];
polort = [ linspace(-1,1,size(tsR,1)); linspace(-1,1,size(tsR,1)).^2 ]';

for multIndexEcc = 1 : length(indexLoop)-1
    for multIndexSize = 1 : length(indexLoop)-1
        for hIndex = 1 : size(hrfMatrixAll,2)
            
            hrfMatrix = hrfMatrixAll( :, hIndex );
            %% build predictors
            
            sigmaPrf = sigmaPrfArray( indexLoop(multIndexSize) : indexLoop(multIndexSize+1)-1  );
            
            eccPrf = eccPrfArray( indexLoop(multIndexEcc) : indexLoop(multIndexEcc+1)-1  );
                        
            sOut = [ sOut sigmaPrf ];
            eOut = [ eOut eccPrf ];
            
            params.surroundGauss = 0;
            params.vectorVsLoop = 2;
            params.hrfPar = hrfMatrix;
            params.normalizePredictor = 1;
            tssP = single( zeros( size(tsR,1), length(sigmaPrf) * length(eccPrf) * size( hrfMatrix, 2 ) * length( multNegative ) ) );
            tssN = single( zeros( size(tsR,1), length(sigmaPrf) * length(eccPrf) * size( hrfMatrix, 2 ) * length( multNegative ) ) );
            dataFrame = single( zeros( length(sigmaPrf) * length(eccPrf) * size( hrfMatrix, 2 ) * length(multNegative) , 4 ) );
            
            disp( sprintf('build predictors...') )
            warning off
            
            count = 1;
            progress = 0.05;
            for k = 1 : length( sigmaPrf )
                for f = 1 : length( eccPrf )
                    for kk = 1 : length( multNegative )
                        for l = 1:size( hrfMatrix, 2 )
                            
                            params.surroundGauss = 0;
                            params.sigmaPrf = sigmaPrf(k);
                            params.eccPrf = eccPrf(f);
                            params.hrfPar = hrfMatrix(:,l);
                            singlePredictorPos = single( generatePredictor1D(params,stimuli)' );
                            
                            
                            params.surroundGauss = 1;
                            params.sigmaNegativePrf = params.sigmaPrf*multNegative(kk);
                            singlePredictorNeg = single( generatePredictor1D(params,stimuli)' );
                            
                            %singlePredictor = singlePredictorPos + singlePredictorNeg;
                            %singlePredictor = singlePredictor./max(singlePredictor);
                            
                            if isnan( sum( singlePredictorPos  ) )
                                singlePredictorPos = single( zeros( size(tsR,1), 1 ) );
                            end
                            if isnan( sum( singlePredictorNeg  ) )
                                singlePredictorNeg = single( zeros( size(tsR,1), 1 ) );
                            end
                            
                            tssP(:,count) = singlePredictorPos;
                            tssN(:,count) = singlePredictorNeg;
                            dataFrame(count,:) = [ sigmaPrf(k) eccPrf(f) l multNegative(kk) ];
                           
                            count = count + 1;
                            if count > size( tssP, 2 )*progress
                                fprintf(1,'.');drawnow;
                                progress = progress + 0.05;
                            end
                            
                        end
                    end
                end
                
            end
            
            warning on          
            
            %% fit model
            params.runSurroundFit = 1;
            params.saveOuput = 1;
            
            if countFitLoop == 1
                output = fitPrfModelAcrossLayers_v5(params, tsR, tssP, tssN, dataFrame, stimuli, sigmaPrf, eccPrf, hrfMatrix, polort);
                output.hrfIndex = ones(1,length(output.sigmaArray));
                output.loop = ones(1,length(output.sigmaArray));
            else
                outputTemp = fitPrfModelAcrossLayers_v5(params, tsR, tssP, tssN, dataFrame, stimuli, sigmaPrf, eccPrf, hrfMatrix, polort);
                indexUpdateValues = find( outputTemp.r2Array >= output.r2Array & ... 
                                    outputTemp.betaSlopeArray1 > outputTemp.betaSlopeArray2 );            
                %% update the output structure with the new values
                output.loop( indexUpdateValues ) = countFitLoop;
                output.r2Array( indexUpdateValues ) = outputTemp.r2Array( indexUpdateValues );
                output.sigmaArray( indexUpdateValues ) = outputTemp.sigmaArray( indexUpdateValues );
                output.eccArray( indexUpdateValues ) = outputTemp.eccArray( indexUpdateValues );
                output.betaInterceptArray( indexUpdateValues ) = outputTemp.betaInterceptArray( indexUpdateValues );
                output.betaSlopeArray1( indexUpdateValues ) = outputTemp.betaSlopeArray1( indexUpdateValues );
                output.betaSlopeArray2( indexUpdateValues ) = outputTemp.betaSlopeArray2( indexUpdateValues );
                output.betaPolortArray( indexUpdateValues, : ) = outputTemp.betaPolortArray( indexUpdateValues, : );
                output.hrfIndex( indexUpdateValues ) = hIndex;
                output.sigmaArrayNegative( indexUpdateValues ) = outputTemp.sigmaArrayNegative( indexUpdateValues );
            end
            countFitLoop = countFitLoop + 1;
        end
    end
end
output.hrfMatrix = hrfMatrixAll;

%% fit surround
if params.runSurroundFit
    output.sigmaNegative = output.sigmaArray.*output.sigmaArrayNegative;
    output.sigmaPositive = output.sigmaArray;
    output = prfSizeLayersEstimateParametersMexHat_v5( output, params, stimuli );
end

%% save output
switch params.saveOuput
    case {1}
        if params.runSurroundFit
            dataTable = [ output.sigmaArray; output.eccArray; output.r2Array; output.betaInterceptArray; output.betaSlopeArray1; output.betaSlopeArray2; ...
                output.sigmaNegative; output.sigmaPositiveUpdated; output.sigmaNegativeUpdated; output.betaPolortArray' ];
            filenameMat = [ datestr(now, 'mmddyy_HHMMSS'), '.mat' ];
            filenameTxt = [ filenameMat(1:end-4), '.txt' ];
            filenameBrik = [ 'prfModelOutput_', filenameMat(1:end-4) ];
            save( filenameMat, 'output');
            dlmwrite( filenameTxt, dataTable, 'delimiter', '\t', 'precision', 6 );
            
            emptyResults = single( zeros( [ size( brainMask ), 8 ] ) );
            for k=1:8
                tempResults = zeros( [ size( brainMask ) ] );
                if k==1; tempResults( inplaneRoiIndex ) = output.r2Array; end
                if k==2; tempResults( inplaneRoiIndex ) = output.sigmaPositiveUpdated; end
                if k==3; tempResults( inplaneRoiIndex ) = output.sigmaNegativeUpdated; end
                if k==4; tempResults( inplaneRoiIndex ) = output.sigmaArray; end
                if k==5; tempResults( inplaneRoiIndex ) = output.eccArray; end
                if k==6; tempResults( inplaneRoiIndex ) = output.betaInterceptArray; end
                if k==7; tempResults( inplaneRoiIndex ) = output.betaSlopeArray1; end
                if k==8; tempResults( inplaneRoiIndex ) = output.betaSlopeArray2; end
                emptyResults(:,:,:,k) = tempResults;
            end
            Opt.Prefix = filenameBrik;
            WriteBrik( emptyResults, infoTs, Opt )
            system( sprintf('mv %s+orig.BRIK ../anatomy', filenameBrik ) );
            system( sprintf('mv %s+orig.HEAD ../anatomy', filenameBrik ) );
            
            
            params = output.params;
            stimuli = output.stimuli;
            tsExpected = single( zeros( tsSize ) ) ;
            tsObserved = single( zeros( tsSize ) ) ;
           
            for k = 1:length( inplaneRoiIndex )
                
                voxelIndex = k;
                Y = output.tsR(:,voxelIndex);
                bSlope1 = output.betaSlopeArray1(voxelIndex);
                bSlope2 = output.betaSlopeArray2(voxelIndex);
                bIntercept = output.betaInterceptArray(voxelIndex);
                polort1 = output.betaPolortArray(voxelIndex,1);
                polort2 = output.betaPolortArray(voxelIndex,2);
                
                params.surroundGauss = 0;
                params.sigmaPrf = output.sigmaPositive(voxelIndex);
                params.eccPrf = output.eccArray(voxelIndex);
                params.hrfPar = output.hrfMatrix(:,1);
                singlePredictorPos = single( generatePredictor1D( params, stimuli )' );
                
                params.surroundGauss = 1;
                params.sigmaNegativePrf = output.sigmaNegative(voxelIndex);
                singlePredictorNeg = single( generatePredictor1D(params,stimuli)' );
                
                pTs = [bIntercept bSlope1 bSlope2 polort1 polort2] * [ ones( 1, size(output.tsR,1) ); singlePredictorPos; singlePredictorNeg; polort(:,1)'; polort(:,2)'  ];
                
                tsExpected( inplaneRoiCoords(1,k), ...
                    inplaneRoiCoords(2,k), ...
                    inplaneRoiCoords(3,k), : ) = pTs ;
                
                
                tTemp = squeeze( ts( inplaneRoiCoords(1,k), ...
                        inplaneRoiCoords(2,k), ...
                        inplaneRoiCoords(3,k), : ) );
                tTemp = single( tTemp );
                tTemp = ( tTemp./mean(tTemp) ).*100;
                tTempIndex = find(tTemp>200);
                tTemp( tTempIndex ) = 200;
                 
                tsObserved(inplaneRoiCoords(1,k), ...
                    inplaneRoiCoords(2,k), ...
                    inplaneRoiCoords(3,k), :) = tTemp;
                
                progressbar( k, length( inplaneRoiIndex ) );
                
                
                %figure(102), plot(Y,'--k','linewidth',1), hold on, plot( pTs,  '--b', 'linewidth',3 )
                
            end
            
            filenameBrikFit = [ 'prfModelOutput_Fit', filenameMat(1:end-4) ];
            Opt.Prefix = filenameBrikFit;
            WriteBrik( double( tsExpected ), infoTs, Opt )
            system( sprintf('mv %s+orig.BRIK ../anatomy', filenameBrikFit ) );
            system( sprintf('mv %s+orig.HEAD ../anatomy', filenameBrikFit ) );
            
            filenameBrikFit = [ 'prfModelOutput_Ts', filenameMat(1:end-4) ];
            Opt.Prefix = filenameBrikFit;
            WriteBrik( double( tsObserved ), infoTs, Opt )
            system( sprintf('mv %s+orig.BRIK ../anatomy', filenameBrikFit ) );
            system( sprintf('mv %s+orig.HEAD ../anatomy', filenameBrikFit ) );
            
            
        else
            dataTable = [ output.sigmaArray; output.eccArray; output.r2Array; output.distanceArray'; output.betaInterceptArray; output.betaSlopeArray; ...
                output.hrfMatrixParameterFWHM; output.hrfMatrixParameterResponseDelay; output.sigmaPositiveUpdatedFWHM ];
            filename = [ params.saveDirectory, datestr(now, 'mmddyy_HHMMSS'), params.ROIname, '.mat' ];
            save( filename, 'output');
            filenameTxt = [ params.saveDirectory, datestr(now, 'mmddyy_HHMMSS'), params.ROIname, '.txt' ];
            dlmwrite( filenameTxt, dataTable, 'delimiter', '\t', 'precision', 6 );
        end
end

overallTime = toc;

cd ..

% sum( output.sigmaPositiveUpdated==1000 )
% 
% %% plot single ts
% [sorted indx] = sort(  output.r2Array, 'descend' );
% params = output.params;
% stimuli = output.stimuli;
% 
% voxelIndex = indx(1);
% Y = output.tsR(:,voxelIndex);
% bSlope1 = output.betaSlopeArray1(voxelIndex);
% bSlope2 = output.betaSlopeArray2(voxelIndex);
% bIntercept = output.betaInterceptArray(voxelIndex);
% polort1 = output.betaPolortArray(voxelIndex,1);
% polort2 = output.betaPolortArray(voxelIndex,2);
% 
% params.surroundGauss = 0;
% params.sigmaPrf = output.sigmaPositive(voxelIndex);
% params.eccPrf = output.eccArray(voxelIndex);
% params.hrfPar = output.hrfMatrix(:,1);
% singlePredictorPos = single( generatePredictor1D( params, stimuli )' );
% 
% params.surroundGauss = 1;
% params.sigmaNegativePrf = output.sigmaNegative(voxelIndex);
% singlePredictorNeg = single( generatePredictor1D(params,stimuli)' );
% 
% pTs = [bIntercept bSlope1 bSlope2 polort1 polort2] * [ ones( 1, size(output.tsR,1) ); singlePredictorPos; singlePredictorNeg; polort(:,1)'; polort(:,2)'  ];
% 
% figure(102), plot(Y,'--k','linewidth',1), hold on, plot( pTs,  '--b', 'linewidth',3 )
% 
% modelResiduals = Y - pTs';
% varResiduals = var( modelResiduals, 0, 1);
% varVolume = var( Y );
% r2 = single( 1 - ( ( varResiduals ) ./ varVolume ) );
% output.r2Array( voxelIndex )
% 
% 
% 



% output.sigmaNegative = output.sigmaArray.*output.sigmaArrayNegative;
% output.sigmaPositive = output.sigmaArray;
% output = prfSizeLayersEstimateParametersMexHat_v5( output, params, stimuli );
% 
% indExcl = find( output.sigmaPositiveUpdated==1000 );
% r2Excl = output.r2Array( output.sigmaPositiveUpdated==1000 );
% 
% [sorted indx] = sort(  r2Excl, 'descend' );
% output.sigmaPositiveUpdated( indExcl( indx( 4 ) ) )
% output.sigmaNegativeUpdated( indExcl( indx( 4 ) ) )
% output.loop( indExcl( indx( 4 ) ) )
% 
% for k = 1:25
%     voxelIndex = indExcl( indx( k ) );
%     Y = output.tsR(:,voxelIndex);
%     bSlope1 = output.betaSlopeArray1(voxelIndex);
%     bSlope2 = output.betaSlopeArray2(voxelIndex);
%     bIntercept = output.betaInterceptArray(voxelIndex);
%     polort1 = output.betaPolortArray(voxelIndex,1);
%     polort2 = output.betaPolortArray(voxelIndex,2);
%     
%     params.surroundGauss = 0;
%     params.sigmaPrf = output.sigmaPositive(voxelIndex);
%     params.eccPrf = output.eccArray(voxelIndex);
%     params.hrfPar = output.hrfMatrix(:,1);
%     singlePredictorPos = single( generatePredictor1D( params, stimuli )' );
%     
%     params.surroundGauss = 1;
%     params.sigmaNegativePrf = output.sigmaNegative(voxelIndex);
%     singlePredictorNeg = single( generatePredictor1D(params,stimuli)' );
%     
%     pTs = [bIntercept bSlope1 bSlope2 polort1 polort2] * [ ones( 1, size(output.tsR,1) ); singlePredictorPos; singlePredictorNeg; polort(:,1)'; polort(:,2)'  ];
%     
%     figure(102),
%     subplot(5,5,k), plot(Y,'--k','linewidth',1), hold on, plot( pTs,  '--b', 'linewidth',3 )
%     
%     modelResiduals = Y - pTs';
%     varResiduals = var( modelResiduals, 0, 1);
%     varVolume = var( Y );
%     r2 = single( 1 - ( ( varResiduals ) ./ varVolume ) )
%     output.r2Array( voxelIndex )
% end
% 
% %hist( output.r2Array, 1000 )
% %hist( output.betaSlopeArray1(output.r2Array>0.25), 10000 )
