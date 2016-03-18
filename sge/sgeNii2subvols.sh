#!/bin/bash
# nii2subvols
#
# sgeNii2subvols.sh dataName subvolVer <dsRate>
#
# examples
# ./sgeNii2subvols.sh buckner mar10_2016 5
# ./sgeNii2subvols.sh ADNI_T1_baselines mar10_2016 5

if [ "$#" -lt 2 ] ; then
  echo "Usage: $0 dataName subvolVer dsFact" >&2
  exit 1
fi

###############################################################################
# Parameters
###############################################################################

# this version's info
dataName=$1 # stroke or buckner
subvolVer=$2
if [ "$#" -lt 3 ] ; then dsRate="7"; else dsRate="$3"; fi

procType="wholevol" # wholevol brain_pad10 brain_pad30

# prepare data
PROC_PATH="/data/vision/polina/projects/stroke/work/patchSynthesis/data/${dataName}/proc/${procType}/";
ATLAS_PATH="/data/vision/polina/projects/stroke/work/patchSynthesis/data/${dataName}/atlases/${procType}/";
ATLAS_FILE_SUFFIX="${dataName}61";

###############################################################################
# Program Settings
###############################################################################

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
OUTPUT_PATH="/data/vision/polina/projects/stroke/work/patchSynthesis/data/${dataName}/subvols/${subvolVer}/${procType}/";
OUTPUT_PATH="/data/vision/polina/scratch/adalca/patchSynthesis/tmp-data/${dataName}/subvols/${subvolVer}/${procType}/";
mkdir -p ${OUTPUT_PATH}
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/"
CLUST_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/MCC/";

# command shell file
mccSh="${CLUST_PATH}MCC_mccNii2subvols/run_mccNii2subvols.sh"

###############################################################################
# Running Code
###############################################################################

# execute
for subjfolder in `ls ${PROC_PATH}`
do
  subjid=`echo $subjfolder | cut -d _ -f 1`

  niifile="${PROC_PATH}${subjid}/${subjid}_ds${dsRate}_us${dsRate}_reg.nii.gz"
  volName="ds${dsRate}us${dsRate}reg"
  patchSize="[9,9,9]"
  gridSpacing="[5,5,5]"
  atlVolSize="[256,256,256]"
  savefile="${OUTPUT_PATH}${subjid}/${subjid}_ds${dsRate}_us${dsRate}_reg_loc%d.mat"

  checkfile="${OUTPUT_PATH}${subjid}/${subjid}_ds${dsRate}_us${dsRate}_reg_loc10000.mat"
  if [ -f $checkfile ]; then continue; fi

  # nii2subvols(niifile, volName, patchSize, gridSpacing, atlVolSize, savefile)
  lcmd="${mccSh} $mcr $niifile $volName $patchSize $gridSpacing $atlVolSize $savefile"

  # create sge file
  sgeopath="${OUTPUT_PATH}${subjid}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=100G \""
  sge_par_q="" #--sge \"-q qOnePerHost \""
  sgerunfile="${sgeopath}/mccNii2subvols.sh"
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${lcmd} > ${sgerunfile}"
  echo $cmd
  eval $cmd
  chmod a+x ${sgerunfile}

  # run sge
  sgecmd="qsub ${sgerunfile}"
  echo $sgecmd
  $sgecmd

  # sleep for a bit to give sge time to deal with the new job (?)
  # sleep 1
done
