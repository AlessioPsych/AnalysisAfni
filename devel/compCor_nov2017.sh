#!/bin/bash

EPIDIR=$1
CLFRAC=$2
EPIDIRWRITE=$3
FILETYPEINPUT=$4
POLORTPAR=$5
STIMULUSDIR=$6
JUSTDETRENDFLAG=$7

if [ -z "$1" ]
then
echo 'Bash script to correct EPI data for compcor components.' 
echo 'using ANTsR, 4 components are selected, based on 99 percentile of time series stds'
echo 'EPIDIR=$1, epi files directory (put the backslash, as usual)' 
echo 'CLFRAC=$2, clfrac parameter for 3dAutomask, use 0.5 to begin, small values -> larger volume' 
echo 'EPIDIRWRITE=$3, newly created directory where the output is stored'
echo 'FILETYPEINPUT=$4, input file type (.nii for example)'
echo 'POLORPAR=$5, polort parameter for detrending'
echo 'STIMULUSDIR=$6, directory where the stimuli file (.1D is stored)'
echo 'JUSTDETRENDFLAG=$7, do you want just to detrend the data but not correct for compcor?'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/compCor_nov2017.R \
$EPIDIR $CLFRAC $EPIDIRWRITE $FILETYPEINPUT $POLORTPAR $STIMULUSDIR $JUSTDETRENDFLAG $AFNI_TOOLBOXDIR

