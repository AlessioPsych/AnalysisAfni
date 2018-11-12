#!/bin/bash

PARAMFILE=$1

if [ -z "$1" ]
then
echo 'parameters file with this structure'
echo ''
echo 'baselineLengthStep 20, (TRs)'
echo 'baselineLengthBegin 8, (TRs)'
echo 'baselineLengthEnd 8, (TRs)'
echo 'bWidth 1, (dva)'
echo 'bStep 0.5, (dva)'
echo 'nSteps 20 note that bStep*bWidth=screenWidth'
echo 'screenWidth 10 (dva)'
echo ''
echo 'from here on there is the order of the stimuli baselines + bars'
echo 'the names define the direction of the bar, for example'
echo ''
echo 'baselineStepBegin'
echo 'stimDown'
echo 'baselineStep'
echo 'stimDownAngleSx'
echo 'stimRight'
echo 'baselineStep'
echo 'stimDownAngleDx'
echo 'stimUp'
echo 'baselineStep'
echo 'stimUpAngleDx'
echo 'stimLeft'
echo 'baselineStep'
echo 'stimUpAngleSx'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/generatePrfStimuli.R $PARAMFILE
