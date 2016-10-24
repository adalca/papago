#!/bin/bash
#
# ./loopBuckner.sh sgeSubvolReconGlobalGMM

sgefcn=$1

for f in /data/vision/polina/scratch/adalca/patchSynthesis/data/buckner/proc/*
do
  #echo $f
  cmd="./${sgefcn}.sh /data/vision/polina/scratch/adalca/patchSynthesis/data/adni ${f} dec2015"
  echo $cmd
  
  $cmd
done
