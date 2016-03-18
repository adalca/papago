
#!/bin/bash
# prepare useful volumes for a subject in the data set
#
# >$ sgePreparePaffineSubjVolumes.sh datapath subjid
#
# examples:
# ./sgePreparePaffineSubjVolumes.sh /data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/proc/wholevol buckner01

###############################################################################
# Settings
###############################################################################

if [ "$#" -lt 1 ] ; then
  echo "Usage: $0 dataPath <subjid>" >&2
  exit 1
fi

OPATH=$1

if [ "$#" -ne 2 ]; then
  subjids=`ls ${OPATH}`
else
  subjids=$2
fi

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/";

# training shell file
mccSh="${PROJECT_PATH}/../MCC/MCC_mccPreparePaffineSubjVolumes/run_mccPreparePaffineSubjVolumes.sh"

###############################################################################
# Running Code
###############################################################################

for subjid in ${subjids}
do
  # subject files
  PATH_TEST_SUBJ="${OPATH}/${subjid}"
  dsSubjFile="${PATH_TEST_SUBJ}/${subjid}_ds5_us5.nii.gz"
  dsSubjInAtlFile="${PATH_TEST_SUBJ}/${subjid}_ds5_us5_reg.nii.gz"
  dsSubjInAtlFileTform="${PATH_TEST_SUBJ}/${subjid}_ds5_us5_reg.mat"
  subjCorrFile="${PATH_TEST_SUBJ}/${subjid}_cor_2_ds5_us5_size.mat"

  if [ -f ${subjCorrFile} ] ;
  then
    echo skipping $subjid since $subjCorrFile exists
    continue;
  fi

  # prepare command
  lcmd="${mccSh} $mcr $dsSubjInAtlFileTform $dsSubjFile $dsSubjInAtlFile $subjCorrFile"

  # create sge files
  sgeopath="${PATH_TEST_SUBJ}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=100G \""
  sgerunfile="${sgeopath}/preparePaffineSubjVolumes_${subjid}_cor_2_ds5_us5_size.sh"

  # prepare sge run file
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l ${lcmd} > ${sgerunfile}"
  echo $cmd
  eval $cmd

  # run training
  sgecmd="qsub ${sgerunfile}"
  echo -e "$sgecmd\n"
  $sgecmd

done
