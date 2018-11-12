#!/bin/bash

export COORDNAME=$1
export TOPONAME=$2
export CURRDIR=$PWD

if [ -z "$1" ]
then
echo
echo 'flattens a patch from the a surface, to extract a patch see file extractPatch.sh'
echo
echo 'Inputs:'
echo 'COORDNAME=$1, coordinates file (vertices)'
echo 'TOPONAME=$2, topology file (faces)'
echo
exit 1
fi

cd $AFNI_TOOLBOXDIRSURFACES

matlab -nodisplay -nojvm -r 'flattening; exit;'

cd $CURRDIR