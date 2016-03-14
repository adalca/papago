#!/bin/bash
# run subspace training
#
# >$ sgeTrain datapath gmmver
#
# examples:
# >$ /data/vision/polina/scratch/adalca/patchSynthesis/data/adni dec2015
# >$ /data/vision/polina/scratch/adalca/patchSynthesis/data/buckner dec2015

###############################################################################
# Settings
###############################################################################

DATA_PATH=$1
gmmver=$2

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/";
SUBVOLS_PATH="${DATA_PATH}/subvols/"
GMMS_PATH="${DATA_PATH}/gmms/"

# training shell file
mccSh="${PROJECT_PATH}/MCC/MCC_sgeTrainCluster/run_sgeTrainCluster.sh"

# training parameters
K=(2 3 5 10 15 25)
idxvec=("1" "[3, 2]") # iso, ds
modelvec=("model0" "model3") # iso, ds
nSamples="30000"
nGrid=`ls ${SUBVOLS_PATH}/*.mat | wc -l`
nRunTypes="${#idxvec[@]}"

###############################################################################
# Running Code
###############################################################################

# run different settings
for r in `seq 0 ${nRunTypes}`; # different models
do
  model=${modelvec[$r]}

  for k in ${K[@]}
  do

    # prepare output folder.
    gmmfolder="${GMMS_PATH}/${gmmver}/${model}/K${k}/"
    sgeopath="${gmmfolder}/sge/"
    mkdir -p ${sgeopath}

    for i in `seq 1 $nGrid`
    do

      # prepare input & output mat files
      subvolfile="${SUBVOLS_PATH}/subvol_${i}.mat"
      gmmfile="${gmmfolder}/gmm_${i}_${model}_K${k}.mat"

      # sge file which we will execute
      sgerunfile="${sgeopath}/gmm_${i}_${model}_K${k}.sh"

      # prepare matlab command
      lcmd="${mccSh} $mcr $subvolfile $k ${model} $gmmfile ${idxvec[$r]} $nSamples"

      # create sge files
      sge_par_o="--sge \"-o ${sgeopath}\""
      sge_par_e="--sge \"-e ${sgeopath}\""

      cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e ${lcmd} > ${sgerunfile}"
      echo $cmd
      eval $cmd

      # run training
      sgecmd="qsub ${sgerunfile}"
      echo -e "$sgecmd\n"
      $sgecmd

      echo "done training queue r:${r} k:${K} i:${i}"
      # sleep 1
    done
  done
done
