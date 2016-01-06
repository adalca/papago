#!/bin/bash

# SGE
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}
mkdir -p $SGE_LOG_PATH
echo $SGE_O_HOME

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/latest/subspace/"
CLUST_PATH="${PROJECT_PATH}/clust/"
SYNTHESIS_DATA_PATH="/data/vision/polina/scratch/adalca/patchSynthesis/data/"
ADNI_SUBVOLS_PATH="${SYNTHESIS_DATA_PATH}/adni/subvols_v2/"
ADNI_GMMS_PATH="${SYNTHESIS_DATA_PATH}/adni/gmms/"
ADNI_TMP_PATH="${SYNTHESIS_DATA_PATH}/adni/tmp/"

# GMMLL shell file
mccSh="${PROJECT_PATH}/MCC/MCC_sgeGMMLL/run_sgeGMMLL.sh"

# training parameters
K=(2 3 5 10 15 25)
volnamesvec=("brainIso2Ds5Us5sizeReg" "brainDs5Us5Reg,brainDs5Us5RegMask")
modelvec=("model0" "model3") # iso, ds
subtractpatchmean="true"
nSamples="5000"
nRunTypes="${#volnamesvec[@]}"
nGrid=`ls ${ADNI_SUBVOLS_PATH}/*.mat | wc -l`

# run different models
for r in `seq 0 ${nRunTypes}`;
do
  model=${modelvec[$r]}
  volnames=${volnamesvec[$r]}
  echo $model $r

  for K in ${K[@]}
  do

    for i in `seq 1 $nGrid`
    do
      # prepare input & output mat files
      gmmfile="${ADNI_GMMS_PATH}/gmm_${i}_${model}_K${K}.mat"
      subvolfile="${ADNI_SUBVOLS_PATH}/subvol_${i}.mat"
      llfile="${ADNI_TMP_PATH}/ll_${i}_${model}_K${K}.mat"

      # run training.
      cmd="${CLUST_PATH}/qsub-run ${mccSh} $mcr $gmmfile $subvolfile $volnames $subtractpatchmean $nSamples $llfile"
      echo $cmd
      $cmd

      echo "done gmmll r:${r} k:${K} i:${i}"
      sleep 100
    done
  done
done
