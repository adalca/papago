#!/bin/bash
# run subvolume reconstruction
#
# >$ sgeSubvolReconGlobalGMM traindatapath testsubjpath ver <inireconfile>
#
# examples:
# ./sgeSubvolReconGlobalGMM.sh /data/vision/polina/scratch/adalca/patchSynthesis/data/adni /data/vision/polina/scratch/adalca/patchSynthesis/data/buckner/proc/buckner01 dec2015

###############################################################################
# Settings
###############################################################################

PATH_TRAIN=$1
PATH_TEST_SUBJ=$2
gmmver=$3

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/";
SUBVOLS_PATH="${PATH_TRAIN}/subvols/${gmmver}"
GMMS_PATH="${PATH_TRAIN}/gmms/${gmmver}"
iniReconFile="${PROJECT_PATH}/ini/recon.ini" # default recon ini file

# subject files
subjid=`basename ${PATH_TEST_SUBJ}`;
dsSubjInAtlFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated5_reg.nii.gz"
dsSubjInAtlMaskFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated5_dsmask_reg.nii.gz"
dsSubjFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated5.nii.gz"
dsSubjWeightFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated5_dsmask.nii.gz"
subjCorrFile="${PATH_TEST_SUBJ}/${subjid}_brain_cor_2_ds5_us2_size.mat"
PATH_RECONSUBVOLS="${PATH_TEST_SUBJ}/subvolRecons"

# training shell file
mccSh="${PROJECT_PATH}/../MCC/MCC_mccSubvolRecon/run_mccSubvolRecon.sh"

# recon parameters -- these are the parameters that should have been trained already
K=(1 2)
modelvec=("global") #

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
    subjoutFolder="${PATH_RECONSUBVOLS}/${gmmver}/${model}/gGMM${k}/"
    sgeopath="${subjoutFolder}/sge/"
    mkdir -p ${sgeopath}

    gmmfolder="${GMMS_PATH}/${model}/gGMM${k}/"

    for i in `seq 1 $nGrid`
    do

      # prepare input files
      gmmfile="${gmmfolder}/ggmm.mat"
      subvolfile="${SUBVOLS_PATH}/subvol_${i}.mat"
      subjoutFile="${subjoutFolder}/subvolRecon_${i}_${model}_gGMM${k}.mat"

      # sge file which we will execute
      sgerunfile="${sgeopath}/rec_${i}_${model}_gGMM${k}.sh"

      # prepare matlab command
      lcmd="${mccSh} $mcr $gmmfile $subvolfile $iniReconFile $dsSubjInAtlFile $dsSubjInAtlMaskFile $dsSubjFile $dsSubjWeightFile $subjCorrFile $subjoutFile"

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

      echo "done subvolRecon queue r:${r} gGMM:${K} i:${i}"
      #sleep 1
    done
  done
done
