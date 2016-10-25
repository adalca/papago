#!/bin/bash
# sgeHugeMatfile2subvol
#
# examples
# ./sgeHugeMatfile2subvol.sh ADNI_T1_baselines mar12_2016 mod wholevol 

if [ "$#" -lt 4 ] ; then
  echo "Usage: $0 dataName subvolVer mod proctype <dsFact>" >&2
  exit 1
fi

# this version's info
dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
mod=$3 # e.g. Ds5Us5RegMask
procType=$4 # wholevol brain_pad10 brain_pad30

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
mkdir -p ${SUBVOLFILES_PATH}
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/"
CLUST_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/MCC/";

# command shell file
mccSh="${CLUST_PATH}MCC_hugeMatfile2subvol/run_hugeMatfile2subvol.sh"

# files
matfilefile="${OUTPUT_PATH}/${dataName}_${procType}_${mod}_volumes.mat"
locfile="${OUTPUT_PATH}/selidx2loc_ds7us5.txt"
subvolSize="[18,18,18]" # patchSize + gridSpacing - 2 = 13 + 7 - 2

###############################################################################
# Running Code
###############################################################################

# execute
while read line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  subvolLoc=`echo ${line} | cut -d " " -f 2`
  subvolfile="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"

  # if [ -f $subvolfile ] ; then
  #   printf "skipping $subvolInd since $subvolfile is present \n\n"
  #   continue;
  # fi

  # hugeMatfile2subvol(niifile, volName, patchSize, gridSpacing, atlVolSize, savefile)
  lcmd="${mccSh} $mcr $matfilefile $subvolLoc $subvolSize $subvolfile"

  # create sge file
  sgeopath="${OUTPUT_PATH}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=100G \""
  sge_par_q="" #--sge \"-q qOnePerHost \""
  sgerunfile="${sgeopath}/mccHugeMatfile2subvol_${subvolInd}.sh"
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${lcmd} > ${sgerunfile}"
  # echo $cmd
  eval $cmd
  chmod a+x ${sgerunfile}

  # run sge
  sgecmd="qsub ${sgerunfile}"
  echo $sgecmd
  $sgecmd

  # sleep for a bit to give sge time to deal with the new job (?)
  # sleep 100
done < ${locfile}
