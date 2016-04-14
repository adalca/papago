#!/bin/bash
# sgeSubvol2ecmwgmm
#
# examples
# ./sgeSubvol2ecmwgmm.sh ADNI_T1_baselines mar12_2016 Ds5Us5Reg wholevol 5

if [ "$#" -lt 4 ] ; then
  echo "Usage: $0 dataName subvolVer mod proctype <dsFact>" >&2
  exit 1
fi

# this version's info
dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
mod=$3 # e.g. Ds5Us5Reg. Mask will be added for the weight.
procType=$4 # wholevol brain_pad10 brain_pad30
if [ "$#" -lt 4 ] ; then dsRate="5"; else dsRate="$4"; fi

###############################################################################
# Paths
###############################################################################

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
OUTPUT_PATH="/data/vision/polina/projects/stroke/work/patchSynthesis/data/${dataName}/subvols/${procType}/${subvolVer}/";
SUBVOLFILES_PATH="${OUTPUT_PATH}subvols/"
RECONFILES_PATH="${OUTPUT_PATH}recons/"
mkdir -p ${RECONFILES_PATH}
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/"
CLUST_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/MCC/";

# command shell file
mccSh="${CLUST_PATH}MCC_subvol2ecmwgmm/run_subvol2ecmwgmm.sh"

# files
inifilename="${PROJECT_PATH}/ini/subvol2ecmwgmm.ini"
locfile="${OUTPUT_PATH}/selidx2loc_rest343.txt"

###############################################################################
# Running Code
###############################################################################

# execute
while read line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  dsSubvolMat="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wtSubvolMat="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}Mask_subvol${subvolInd}.mat"
  clusterIdxMat="${RECONFILES_PATH}clusterIdx_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wgmmMat="${RECONFILES_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"

  # subvol2ecmwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
  lcmd="${mccSh} $mcr $dsSubvolMat $wtSubvolMat $clusterIdxMat $wgmmMat $inifilename"

  # create sge file
  sgeopath="${OUTPUT_PATH}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=100G \""
  sge_par_q="--sge \"-q qSparse \""
  sgerunfile="${sgeopath}/subvol2ecmwgmm_${subvolInd}.sh"
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${lcmd} > ${sgerunfile}"
  echo $cmd
  eval $cmd
  chmod a+x ${sgerunfile}

  # run sge
  sgecmd="qsub ${sgerunfile}"
  echo $sgecmd
  $sgecmd

  # sleep for a bit to give sge time to deal with the new job (?)
  # sleep 100
done < ${locfile}
