#!/bin/bash

INPUTVOL=$1
SAVEDATA=$2
OUTDIR=$3
FILENAME=$4

if [ -z "$1" ]
then
echo
echo 'Creates a levelset surface representations from a probabilistic or deterministic tissue' 
echo 'classification. The levelset indicates each voxel’s distance to the closest boundary. It takes' echo 'negative values inside and positive values outside of the brain.'
echo 'from the toolbox nighres, check it here: http://nighres.readthedocs.io/en/latest/index.html'
echo
echo 'INPUTVOL=$1, volume filename' 
echo 'SAVEDATA=$2, save output? [True,False]' 
echo 'OUTDIR=$3, output directory'
echo 'FILENAME=$4, output filename'
echo
exit 1
fi

CURRDIR=$PWD
INDIR=$CURRDIR/$INPUTVOL
echo $INDIR

cp $AFNI_TOOLBOXDIR/nighres_addon/probability_to_levelset_wrapper.py $NIGHRES_TOOLBOXDIR/_tttt.py

cd $NIGHRES_TOOLBOXDIR

python _tttt.py $INDIR $SAVEDATA $OUTDIR $FILENAME
rm _tttt.py

cd $CURRDIR





