cd anatomy

[ volLeft ] = BrikLoad( roiFilenameLeft );
[ volRight ] = BrikLoad( roiFilenameRight );
[ distanceMap ] = BrikLoad( distanceMapFilename );
[ prfOutcome ] = BrikLoad( prfOutcomeFilename );

areaIndex = ( volLeft == 1 | volRight == 1 );

distanceArray = distanceMap( areaIndex ); 
storePrfMatrix = zeros( length( distanceArray ), size( prfOutcome, 4 ) );
for k = 1:size( prfOutcome, 4 )
    tempVolume = prfOutcome(:,:,:,k);
    storePrfMatrix(:,k) = tempVolume( areaIndex );
end
storePrfMatrix = [ storePrfMatrix distanceArray ];

filterIndex = storePrfMatrix(:,9) > -0.1 & storePrfMatrix(:,9) <= 2.2 & ...
    storePrfMatrix(:,1) > quantile( storePrfMatrix(:,1), 0.5 ) & ...
    storePrfMatrix(:,2) ~= 1000 & ...
    storePrfMatrix(:,5) < 4.5 & storePrfMatrix(:,5) > 0.5 & ...    
    storePrfMatrix(:,7) > 0 & storePrfMatrix(:,7) < 50 & ...
    storePrfMatrix(:,8) > 0 & storePrfMatrix(:,8) < 50;
%storePrfMatrix(:,7)*(2/3) > storePrfMatrix(:,8) & ...

storePrfMatrixFilter = storePrfMatrix( filterIndex, : );

figure(18)
names{1} = 'coherence';
names{2} = 'sigmaUpdatedCenter';
names{3} = 'sigmaUpdatedSurround';
names{4} = 'sigma';
names{5} = 'eccentricity';
names{6} = 'betaInterceptArray';
names{7} = 'betaSlopeArray1';
names{8} = 'betaSlopeArray2';
names{9} = 'distance';
for k = 1:size( storePrfMatrixFilter, 2 )
    subplot(3,3,k), hist( storePrfMatrixFilter(:,k), 20 )
    title( names{k} );
end

figure(19)
distanceCut = quantile( storePrfMatrixFilter(:,9), linspace(0,1,11) );
for k = 1 : ( length( distanceCut ) - 1 )
        
        filterDistance = storePrfMatrixFilter( :, 9 ) >= distanceCut( k ) & ...
            storePrfMatrixFilter(:,9) < distanceCut( k + 1 );
        
        sigmaPositiveTemp = storePrfMatrixFilter( filterDistance , 2 );
        sigmaNegativeTemp = storePrfMatrixFilter( filterDistance , 3 );
        eccentricityTemp = storePrfMatrixFilter( filterDistance , 5 );
        eccentricityCut = quantile( eccentricityTemp, linspace(0,1,11) );
        
        for j = 1 : ( length( eccentricityCut ) - 1 )
            filterEccentricity = eccentricityTemp >= eccentricityCut( j ) & ...
                eccentricityTemp < eccentricityCut( j + 1 );
            meanSigmaPositive(j) = mean( sigmaPositiveTemp( filterEccentricity ) );
            meanSigmaNegative(j) = mean( sigmaNegativeTemp( filterEccentricity ) );
            meanEccentricity(j) = mean( eccentricityTemp( filterEccentricity ) );
        end
        subplot(3,4,k), plot( meanEccentricity, meanSigmaPositive, 'o' ); ylim([0 2])
end

nBootRep = 1000;
sigmaPositiveProfile = zeros( length( distanceCut ) - 1, nBootRep );
sigmaNegativeProfile = zeros( length( distanceCut ) - 1, nBootRep );
meanDistance = zeros( length( distanceCut ) - 1, nBootRep );
meanVar = zeros( length( distanceCut ) - 1, nBootRep );
meanBeta = zeros( length( distanceCut ) - 1, nBootRep );
meanIntercept = zeros( length( distanceCut ) - 1, nBootRep );
for k = 1 : ( length( distanceCut ) - 1 )
    
    for n = 1:nBootRep
        
        keepIndex = ceil( length( storePrfMatrixFilter )*rand( length( storePrfMatrixFilter ), 1 ) );
        storePrfMatrixFilterLoop = storePrfMatrixFilter( keepIndex, : );
        %storePrfMatrixFilterLoop = storePrfMatrixFilter;
        
        filterDistance = storePrfMatrixFilterLoop( :, 9 ) >= distanceCut( k ) & ...
            storePrfMatrixFilterLoop(:,9) < distanceCut( k + 1 );
        
        sigmaPositiveTemp = storePrfMatrixFilterLoop( filterDistance , 2 );
        sigmaNegativeTemp = storePrfMatrixFilterLoop( filterDistance , 3 );
        eccentricityTemp = storePrfMatrixFilterLoop( filterDistance , 5 );
        eccentricityCut = quantile( eccentricityTemp, linspace(0,1,11) );
        
        for j = 1 : ( length( eccentricityCut ) - 1 )
            filterEccentricity = eccentricityTemp >= eccentricityCut( j ) & ...
                eccentricityTemp < eccentricityCut( j + 1 );
            meanSigmaPositive(j) = mean( sigmaPositiveTemp( filterEccentricity ) );
            meanSigmaNegative(j) = mean( sigmaNegativeTemp( filterEccentricity ) );
            meanEccentricity(j) = mean( eccentricityTemp( filterEccentricity ) );
        end
        
        
        
        %subplot(3,4,k), plot( meanEccentricity, meanSigmaPositive, 'o' ); ylim([0 2])
        bPositive = glmfit( meanEccentricity, meanSigmaPositive );
        bNegative = glmfit( meanEccentricity, meanSigmaNegative );
        sigmaPositiveProfile(k,n) = bPositive(1) + 2*bPositive(2);
        sigmaNegativeProfile(k,n) = bNegative(1) + 2*bNegative(2);
        meanDistance(k,n) = mean( storePrfMatrixFilterLoop(filterDistance,9) );
        meanVar(k,n) = mean( storePrfMatrixFilterLoop(filterDistance,1) );
        meanBeta(k,n) = mean( storePrfMatrixFilterLoop(filterDistance,7) );
        meanIntercept(k,n) = mean( storePrfMatrixFilterLoop(filterDistance,6) );
        
    end
    
    progressbar( k, ( length( distanceCut ) - 1 ) )
    
end

figure(21)
for k = 1:size( sigmaPositiveProfile, 1 )
    subplot(1,10,k), hist( sigmaPositiveProfile(k,:) );
end
figure(22)
for k = 1:size( sigmaPositiveProfile, 1 )
    subplot(1,10,k), hist( sigmaNegativeProfile(k,:) );
end


sigmaPositiveProfile = mean( sigmaPositiveProfile, 2 );
sigmaNegativeProfile = mean( sigmaNegativeProfile, 2 );
meanDistance = mean( meanDistance, 2 );
meanVar = mean( meanVar, 2 );
meanBeta = mean( meanBeta, 2 );
meanIntercept = mean( meanIntercept, 2 );

positiveProfileFit = smooth( sigmaPositiveProfile, 0.4, 'lowess' );
negativeProfileFit = smooth( sigmaNegativeProfile, 0.4, 'lowess' );
distanceFit = linspace( meanDistance( 1 ), meanDistance( end ), 20 );
posInterp = interp1( meanDistance, positiveProfileFit, distanceFit );
negInterp = interp1( meanDistance, negativeProfileFit, distanceFit );

distanceFitAnalysis = distanceFit( 3 : end-2 );
posInterpAnalysis = posInterp( 3 : end-2 );
negInterpAnalysis = negInterp( 3 : end-2 );

cd ..




figure(20)
subplot(2,2,1), plot( meanDistance, sigmaPositiveProfile, 'o' ); axis square;
    line( meanDistance, positiveProfileFit );
subplot(2,2,2), plot( meanDistance, sigmaNegativeProfile, 'o' ); axis square;
        line( meanDistance, negativeProfileFit );
subplot(2,2,3), plot( meanDistance, meanVar, 'o' ); axis square;
subplot(2,2,4), plot( meanDistance, meanBeta, 'o' ); axis square;

        
linear = linspace(-1,1,length(distanceFitAnalysis));
quadratic = linspace(-1,1,length(distanceFitAnalysis)).^2;

[ bSigmaPositive devSigmaPositive statsSigmaPositive ] = glmfit( [ linear; quadratic ]', posInterpAnalysis );
bSigmaPositive
statsSigmaPositive.p

[ bSigmaNegative devSigmaNegative statsSigmaNegative ] = glmfit( [ linear; quadratic ]', negInterpAnalysis );
bSigmaNegative
statsSigmaNegative.p
