#!/bin/bash

INPUTDIR=$1
FILEFORMAT=$2
POLORT=$3


if [ -z "$1" ]
then
echo 'scales and detrends time series'
echo 'Inputs:'
echo 'INPUTDIR=$1, input directory, put the backslash!'
echo 'FILEFORMAT=$2, file format in regex (like: *.BRIK)' 
echo 'POLORT=$3, -polort parameter for detrending'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/generalPurpose/detrend_whole_dir.R $INPUTDIR $FILEFORMAT $POLORT
