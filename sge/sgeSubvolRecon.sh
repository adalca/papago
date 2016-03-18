#!/bin/bash
# run subvolume reconstruction
#
# >$ sgeSubvolRecon traindatapath testsubjpath ver <inireconfile>
#
# examples:
# ./sgeSubvolRecon.sh ADNI_T1_baselines mar12_2016 101411 wholevol <dsFact>

###############################################################################
# Arguments
###############################################################################

if [ "$#" -lt 4 ] ; then
  echo "Usage  : $0 dataName subvolVer subjid procType <dsFact>" >&2
  echo "Example: $0 ADNI_T1_baselines mar12_2016 101411 wholevol 5" >&2
  exit 1
fi

dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
subjid=$3 # e.g. 101411
procType=$4 # wholevol brain_pad10 brain_pad30
if [ "$#" -lt 5 ] ; then dsRate="5"; else dsRate="$5"; fi
echo $dsRate

###############################################################################
# Main Paths
###############################################################################

DATA_PATH="/data/vision/polina/projects/stroke/work/patchSynthesis/data/"
OUTPUT_PATH="${DATA_PATH}${dataName}/subvols/${procType}/${subvolVer}/"
SUBVOLS_PATH="${OUTPUT_PATH}/subvols/"
RECONS_PATH="${OUTPUT_PATH}/recons/"
SUBJ_PATH="${DATA_PATH}${dataName}/proc/${procType}/${subjid}/"

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/";
# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/
# training shell file
mccSh="${PROJECT_PATH}/../MCC/MCC_mccSubvolRecon/run_mccSubvolRecon.sh"

###############################################################################
# Shared Inputs
###############################################################################

iniReconFile="${PROJECT_PATH}/ini/recon.ini" # default recon ini file

# subject files
subjid=`basename ${SUBJ_PATH}`;
dsSubjInAtlFile="${SUBJ_PATH}/${subjid}_ds5_us5_reg.nii.gz"
dsSubjInAtlMaskFile="${SUBJ_PATH}/${subjid}_ds5_us5_dsmask_reg.nii.gz"
dsSubjFile="${SUBJ_PATH}/${subjid}_ds5_us5.nii.gz"
dsSubjWeightFile="${SUBJ_PATH}/${subjid}_ds5_us5_dsmask.nii.gz"
subjCorrFile="${SUBJ_PATH}/${subjid}_cor_2_ds5_us5_size.mat"

if [ ! -f ${subjCorrFile} ] ; then
  echo "Warning Subject corr file doesn't exist, we'll use tform"
  subjCorrFile="${SUBJ_PATH}/${subjid}_ds5_us5_reg.mat"
fi

locfile="${OUTPUT_PATH}/selidx2loc_rest343.txt"
mod="Ds${dsRate}Us${dsRate}Reg"

###############################################################################
# Running Code
###############################################################################

# run different settings
while read line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  subvolfile="${SUBVOLS_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wgmmfile="${RECONS_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  subjfolder="${RECONS_PATH}/${subjid}/"
  mkdir -p ${subjfolder}
  subjoutFile="${subjfolder}/subvolRecon_${dataName}_${procType}_${mod}_subvol${subvolInd}_${subjid}.mat"

  if [ ! -f $wgmmfile ] ; then
    printf "skipping $subvolInd since I can't find $wgmmfile \n\n"
    continue;
  fi

  if [ -f $subjoutFile ] ; then
    printf "skipping $subvolInd since $subjoutFile is present \n\n"
    continue;
  fi

  # prepare matlab command
  lcmd="${mccSh} $mcr $wgmmfile $subvolfile $iniReconFile $dsSubjInAtlFile $dsSubjInAtlMaskFile $dsSubjFile $dsSubjWeightFile $subjCorrFile $subjoutFile"

  # create sge files
  sgeopath="${OUTPUT_PATH}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=100G \""
  sge_par_q="" #--sge \"-q qSparse \""
  sgerunfile="${sgeopath}/subvolRecon_${subjid}_${subvolInd}.sh"
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${lcmd} > ${sgerunfile}"
  echo $cmd
  eval $cmd

  # run training
  sgecmd="qsub ${sgerunfile}"
  echo -e "$sgecmd\n"
  $sgecmd

  # sleep for a bit to give sge time to deal with the new job (?)
  # sleep 100
done < ${locfile}
