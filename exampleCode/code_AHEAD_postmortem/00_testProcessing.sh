#!/bin/bash

# Instructions to run:
# sudo docker run -v /media/alessiofracasso/DATADRIVE1/AHEAD_exvivo:/AHEAD_exvivo --rm -it nighres
# cd /AHEAD_exvivo
# cd code/
# sh 01_startProcessing.sh

maindir="/AHEAD_exvivo"
targetdir="122017/"

echo "$maindir"
echo "$targetdir"

cd "$maindir"
echo "Current folder: $PWD"


