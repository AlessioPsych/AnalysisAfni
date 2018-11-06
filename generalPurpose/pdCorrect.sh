#!/bin/bash

MPRAGE=$1
PD=$2

if [ -z "$1" ]
 then
  echo 'pd correct anatomy, resamples PD and correct T1 for it. Inputs:'
  echo 'MPRAGE=$1, T1 filename' 
  echo 'PD=$2, PD filename'
  exit 1
fi


3dresample -master ${MPRAGE} -inset ${PD} -prefix PD+orig

3dcalc -a $MPRAGE -b PD+orig -expr 'a / ( b*ispositive(a) )' -prefix MPRAGE_PD+orig

3dAFNItoNIFTI MPRAGE_PD+orig

rm *.BRIK

rm *.HEAD
