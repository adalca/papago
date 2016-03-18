#!/bin/bash
# prepare useful volumes for a subject in the data set
#
# >$ sgePaffineWarpVol.sh datapath subjid
#
# examples:
# ./sgePaffineWarpVol.sh /data/vision/polina/scratch/adalca/patchSynthesis/data/buckner/proc buckner01

###############################################################################
# Settings
###############################################################################

OPATH=$1

if [ "$#" -ne 2 ]; then
  subjids=`ls ${OPATH}`
else
  subjids=$2
fi

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/";

# training shell file
mccSh="${PROJECT_PATH}/../MCC/MCC_mccPaffineWarpVol/run_mccPaffineWarpVol.sh"


###############################################################################
# Running Code
###############################################################################

for subjid in ${subjids}
do
  # subject files
  PATH_TEST_SUBJ="${OPATH}/${subjid}"
  dsSubjFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated2.nii.gz"
  dsusSubjMaskFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated2_dsmask.nii.gz"
  subjCorrFile="${PATH_TEST_SUBJ}/${subjid}_brain_cor_2_ds5_us2_size.mat"
  atlFile="/data/vision/polina/scratch/adalca/patchSynthesis/data/buckner/atlases/buckner61_brain_proc_ds5_us2.nii.gz"
  regOutFile="${PATH_TEST_SUBJ}/${subjid}_brain_downsampled5_reinterpolated2_regwcor.nii.gz"

  # prepare command
  lcmd="${mccSh} $mcr $dsSubjFile $dsusSubjMaskFile $subjCorrFile $atlFile $regOutFile"

  # create sge files
  sgeopath="${PATH_TEST_SUBJ}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sgerunfile="${sgeopath}/paffineWarpVol_${subjid}_brain_cor_2_ds5_us2_size.sh"

  # prepare sge run file
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e ${lcmd} > ${sgerunfile}"
  echo $cmd
  eval $cmd

  # run training
  sgecmd="qsub ${sgerunfile}"
  echo -e "$sgecmd\n"
  $sgecmd

done
