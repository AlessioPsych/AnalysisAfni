#!/bin/bash

CURRDIR=$(pwd)

VOLFILE=$1

cd $AFNI_TOOLBOXDIR
matlab -nodesktop -nosplash -r "matCall('$AFNI_MATLAB', '$VOLFILE', '$CURRDIR')"


