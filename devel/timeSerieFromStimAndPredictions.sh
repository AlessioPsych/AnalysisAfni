#!/bin/bash

MASKFILE=$1
OUTPUTNAME=$2
FINEFIT=$3
SAMPLINGTIME=$4
FITTEDPARAMETERS=$5
TSTEST=$6
STIMTYPE=$7

if [ -z "$1" ]
then
echo 'based on a model built from twoDimensionalFit_Glasgow_eyeFix_gain_parallel_plane.R'
echo 'it takes the prf results from a first stimulus'
echo 'then it takes a second stimuli and build the predicted time series based on this new stimuli, given the prf estimated on the first stimuli' 
echo 'given the prf model fitted on the first stimuli'
echo 'Inputs:'
echo 'MASKFILE=$1, mask file'
echo 'OUTPUTNAME=$2, output name'
echo 'FINEFIT=$3, which model was fitted (first model, first stimuli)'
echo 'SAMPLINGTIME=$4, set it to 0.166'
echo 'FITTEDPARAMETERS=$5: already fitted parameters, with the model defined in FINEFIT'
echo 'TSTEST=$6: time series where you want to test the new time series based on stimuli 2' 
echo 'STIMTYPE=$7 which is the new stimuli (second stimuli) over which you want test the generated time series from stimuli 1'
echo 'it needs the following R libraries:'
echo 'library(pracma) library(abind)'
echo 'library( parallel ) library(neuRosim)'
exit 1
fi

Rscript $AFNI_TOOLBOXDIR/devel/timeSerieFromStimAndPredictions.R $MASKFILE $OUTPUTNAME $FINEFIT $SAMPLINGTIME $FITTEDPARAMETERS $TSTEST $STIMTYPE
