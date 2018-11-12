#!/bin/bash

VOLUME=$1

if [ -z "$1" ]
then

echo
echo
echo 'projects an volume over all the surfaces'
echo 'useful for PSF estimation over the surface'
echo
echo
echo 'Inputs:'
echo 'VOLUME=$1, e.g. volume that you want to'

exit 1
fi


Rscript $AFNI_TOOLBOXDIR/surfaces/vol2Surf.R \
$VOLUME \
$AFNI_TOOLBOXDIR
