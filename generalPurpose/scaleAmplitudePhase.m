cd anatomy

[ amplitude infoAmplitude ] = BrikLoad('amplitudeAnatomy.nii');
[ phase infoPhase ] = BrikLoad('phaseAnatomy.nii');
[ anatomy infoAnatomy ] = BrikLoad('t1Myelin.nii');

amplitude = ( amplitude - mean( amplitude(:) ) ) / std( amplitude(:) );
Opt.Prefix = 'amplitudeAnatomy_scaled';
WriteBrik( amplitude, infoAmplitude, Opt )
system('3dAFNItoNIFTI amplitudeAnatomy_scaled+orig')

phase = ( phase - mean( phase(:) ) ) / std( phase(:) );
Opt.Prefix = 'phaseAnatomy_scaled';
WriteBrik( phase, infoPhase, Opt )
system('3dAFNItoNIFTI phaseAnatomy_scaled+orig')

anatomy = ( anatomy - mean( anatomy(:) ) ) / std( anatomy(:) );
Opt.Prefix = 't1Myelin_scaled';
WriteBrik( anatomy, infoAnatomy, Opt )
system('3dAFNItoNIFTI t1Myelin_scaled+orig')


system('rm *.BRIK')
system('rm *.HEAD')

cd ..