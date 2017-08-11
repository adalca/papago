#!/bin/bash
# run subvolume reconstruction
#
# >$ sgeSubvolRecon traindatapath testsubjpath ver <inireconfile>
#
# examples:
# ./sgeSubvolRecon.sh ADNI_T1_baselines mar12_2016 101411 wholevol iterppca_K5_ppcaK10_apr21 <dsFact>

###############################################################################
# Arguments
###############################################################################

if [ "$#" -lt 5 ] ; then
  echo "Usage  : $0 dataName subvolVer subjid procType <dsFact>" >&2
  echo "Example: $0 ADNI_T1_baselines mar12_2016 101528 wholevol 5" >&2
  exit 1
fi

dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
subjid=$3 # e.g. 101411
procType=$4 # wholevol brain_pad10 brain_pad30
suffix="$5" # e.g. iterppca_K5_ppcaK10_apr21
if [ "$#" -lt 6 ] ; then dsRate="5"; else dsRate="$6"; fi
echo $dsRate

usRate=${dsRate}
# dsRate="7" 
# usRate="5"
mod="Ds${dsRate}Us${usRate}Reg"
# mod="Iso2Ds${dsRate}Us${dsRate}sizeReg"

###############################################################################
# Main Paths
###############################################################################

DATA_PATH="/data/vision/polina/projects/patchSynthesis/data/"
OUTPUT_PATH="${DATA_PATH}${dataName}/subvols/${procType}/${subvolVer}/"
SUBVOLS_PATH="${OUTPUT_PATH}/subvols/"
RECONS_PATH="${OUTPUT_PATH}/recons_${suffix}/"
TRAIN_PATH="${OUTPUT_PATH}/wgmm_${suffix}/"
SUBJ_PATH="${DATA_PATH}${dataName}/proc/${procType}/${subjid}/"

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git-papago/";
# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/
mcr=/data/vision/polina/shared_software/MCR/v91/
# training shell file
mccSh="${PROJECT_PATH}/../MCC/MCC_mccSubvolRecon/run_mccSubvolRecon.sh"

###############################################################################
# Shared Inputs
###############################################################################

oiniReconFile="${PROJECT_PATH}/ini/recon.ini" # default recon ini file
iniReconFile="${RECONS_PATH}/recon.ini" # default recon ini file
mkdir ${RECONS_PATH}
cp ${oiniReconFile} ${iniReconFile}

# subject files
subjid=`basename ${SUBJ_PATH}`;
dsSubjInAtlFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${usRate}_reg.nii.gz"
dsSubjInAtlMaskFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${usRate}_dsmask_reg.nii.gz"
dsSubjFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${usRate}.nii.gz"
dsSubjWeightFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${usRate}_dsmask.nii.gz"
subjCorrFile="${SUBJ_PATH}/${subjid}_cor_2_ds${dsRate}_us${usRate}_size.mat"

if [ ! -f ${subjCorrFile} ] ; then
  echo "Warning Subject corr file doesn't exist, we'll use tform"
  subjCorrFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${dsRate}_reg.mat"
fi

locfile="${OUTPUT_PATH}/selidx2loc_top343.txt"
locfile="${OUTPUT_PATH}/selidx2loc_rest343.txt"
# locfile="${OUTPUT_PATH}/selidx2loc_ds7us5_top343.txt"
# locfile="${OUTPUT_PATH}/selidx2loc_ds7us5_rest343.txt"

###############################################################################
# Running Code
###############################################################################

# run different settings
while read line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  subvolfile="${SUBVOLS_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wgmmfile="${TRAIN_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  subjfolder="${RECONS_PATH}/${subjid}/"
  mkdir -p ${subjfolder}
  subjoutFile="${subjfolder}/subvolRecon_${dataName}_${procType}_${mod}_subvol${subvolInd}_${subjid}.mat"

  if [ ! -f $wgmmfile ] ; then
    # printf "skipping $subvolInd since I can't find $wgmmfile \n\n"
    continue;
  fi

  if [ -f $subjoutFile ] ; then
    # printf "skipping $subvolInd since $subjoutFile is present \n\n"
    continue;
  fi

  # prepare matlab command
  lcmd="${mccSh} $mcr $wgmmfile $subvolfile $iniReconFile $dsSubjInAtlFile $dsSubjInAtlMaskFile $dsSubjFile $dsSubjWeightFile $subjCorrFile $subjoutFile"

  # create sge files
  sgeopath="${RECONS_PATH}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=10G \""
  # sge_par_q="--sge \"-q qSparse \""
  sgerunfile="${sgeopath}/subvolRecon_${subjid}_${subvolInd}.sh"
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${lcmd} > ${sgerunfile}"
  # echo $cmd
  eval $cmd
  chmod a+x ${sgerunfile}
    # ${sgerunfile} % if want to run dirrectly

  # run training
  sgecmd="qsub ${sgerunfile}"
  echo -e "$sgecmd\n"
  $sgecmd

  # sleep for a bit to give sge time to deal with the new job (?)
  # sleep 100
done < ${locfile}
