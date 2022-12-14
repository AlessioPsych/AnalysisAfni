#!/bin/bash

if [ -z "$1" ]
 then
  echo 'DATADIR=$1, directory where individual subjects anatomies are stored'
  echo 'FILEANATOMY=$2, anatomy filename, MUST be common across all subjects' 
  echo 'NCORES=$3, how many cores could be used by Freesurfer' 
  echo 'RUNFLAG=$4, run the code (1) or just show the code (0)' 
  echo ''
  echo 'note:'
  echo 'data structure MUST be like the following:'
  echo 'subjDir01/'
  echo ' ---- FILEANATOMY'
  echo 'subjDir02/'
  echo ' ---- FILEANATOMY'
  echo 'subjDir03/'
  echo ' ---- FILEANATOMY'
  echo 'subjDir04/'
  echo ' ---- FILEANATOMY'
  echo 'etc'
  echo 'etc'
  exit 1
fi

DATADIR=$1
FILEANATOMY=$2
NCORES=$3
RUNFLAG=$4

Rscript runFreesurferSuma.R $DATADIR $FILEANATOMY $NCORES $RUNFLAG