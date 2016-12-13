#!/bin/bash
# run subvolume reconstruction
#
# >$ sgesubvolReconFlexible traindatapath testsubjpath ver <inireconfile>
#
# examples:
# ./sgeSubvolReconFlexible.sh ADNI_T1_baselines mar12_2016 101411 wholevol iterppca_K5_ppcaK10_apr21 latentSubspaceR

###############################################################################
# Arguments
###############################################################################

if [ "$#" -lt 6 ] ; then
  echo "Usage  : $0 dataName subvolVer subjid procType <dsFact>" >&2
  echo "Example: $0 ADNI_T1_baselines mar12_2016 101528 wholevol 5" >&2
  exit 1
fi

dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
subjid=$3 # e.g. 101411
procType=$4 # wholevol brain_pad10 brain_pad30
suffix="$5" # e.g. iterppca_K5_ppcaK10_apr21
reconModel="$6" # 'latentMissingR' or 'latentSubspaceR'

echo $dataName

# set some settings based on datatype
if [ "$dataName" = "stroke" ] ; then
  dsRate="7"
  usRate="5"
  locfilePost="selidx2loc_ds7us5_top343.txt"
  locfilePost="selidx2loc_ds7us5_rest343.txt"
else 
  dsRate="5"
  usRate="5"
  # locfilePost="grididx2loc.txt"
  locfilePost="selidx2loc_top343.txt" 
  locfilePost="selidx2loc_rest343.txt"
fi

# cluster parameters
clustertype='sge'
clustertype='vision'
nVisionCmdQueues=20

# mod name
mod="Ds${dsRate}Us${usRate}Reg"

# possible vision machines
possMachines=("asia" "africa" "america" "europe" "antarctica" "australia" "monday" "tuesday" "wednesday" "thursday" "friday" "saturday" "quickstep" "tango" "waltz" "swing" "rumba" "vision01" "vision02" "vision03" "vision04" "vision05" "vision06" "vision07" "vision08" "vision09" "vision10" "vision11" "vision12" "vision13" "vision14" "vision15" "vision16" "vision17" "vision18" "vision19" "vision20" "vision21" "vision22" "vision23" "vision24" "vision25" "vision26" "vision27" "vision28" "vision29" "vision30" "vision31" "vision32" "vision33" "vision34" "vision35" "vision36" "vision37" "vision38")

###############################################################################
# Main Paths
###############################################################################

DATA_PATH="/data/vision/polina/projects/patchSynthesis/data/"
OUTPUT_PATH="${DATA_PATH}${dataName}/subvols/${procType}/${subvolVer}/"
SUBVOLS_PATH="${OUTPUT_PATH}/subvols/"
RECONS_PATH="${OUTPUT_PATH}/recons_${suffix}_${reconModel}/"
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
mccSh="${PROJECT_PATH}/../MCC/MCC_mccSubvolReconFlexible/run_mccSubvolReconFlexible.sh"

###############################################################################
# Shared Inputs
###############################################################################

oiniReconFile="${PROJECT_PATH}/ini/reconFlexible_${reconModel}.ini" # default recon ini file
iniReconFile="${RECONS_PATH}/reconFlexible_${reconModel}.ini" # default recon ini file
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
  echo "Warning Subject corr file ${subjCorrFile} doesn't exist, we'll use tform"
  subjCorrFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${dsRate}_reg.mat"
fi

# set up locations files
locfile="${OUTPUT_PATH}/${locfilePost}" 

# makse sge path
sgeopath="${RECONS_PATH}/sge/"
mkdir -p ${sgeopath}

###############################################################################
# Running Code
###############################################################################

# prepare batch files
stacki="0"
batchCoreName="`date +%s`_${RANDOM}_${subjid}"
batchfile="${sgeopath}/subvolReconFlexible_batch_${batchCoreName}.sh"

# prepare machine execution
tot="1"
cpui="1"
tot="1"

# run different settings
while read -u10 line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  subvolfile="${SUBVOLS_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wgmmfile="${TRAIN_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  subjfolder="${RECONS_PATH}/${subjid}/"
  mkdir -p ${subjfolder}
  subjoutFile="${subjfolder}/subvolReconFlexible_${dataName}_${procType}_${mod}_subvol${subvolInd}_${subjid}.mat"

  cmdfile="${sgeopath}/subvolReconFlexible_${subjid}_${subvolInd}_cmd.sh"
  echo "# command:" > ${cmdfile}
  chmod a+x $cmdfile

  if [ ! -f $wgmmfile ] ; then
    printf "skipping $subvolInd since I can't find $wgmmfile \n\n"
    continue;
  fi

  if [ -f $subjoutFile ] ; then
    # printf "skipping $subvolInd since $subjoutFile is present \n\n"
    continue;
  fi

  # prepare matlab command
  lcmd="${mccSh} $mcr $wgmmfile $subvolfile $iniReconFile $dsSubjFile $dsSubjWeightFile $subjCorrFile $subjoutFile;"
  echo "$lcmd" > $cmdfile
  echo "$cmdfile;" >> ${batchfile}

  # increment the queue
  stacki=`expr $stacki + 1`
  m=`expr $stacki % $nVisionCmdQueues`
  if [ $m -eq "0" ] ; then
    chmod a+x ${batchfile}

    if [ "$clustertype" = 'vision' ] ; then
      # run on vision machine
      # cpui=`expr $cpui % 35` # if it's 34 then it shoudl stay 34 and call vision35. If it's 35 it should call vision01
      cpui=`expr $cpui + 1`
      cpui=`expr $cpui % ${#possMachines[@]}`
      machine="${possMachines[$cpui]}"
      # machine="vision`echo $cpui | awk '{printf "%02d", $0}'`"
      cmd="ssh ${machine} \"nohup ${batchfile} > ${batchfile}.out 2> ${batchfile}.err &\" " # took out adalca@${machine}
      echo $cmd
      eval $cmd

      tot=`expr $tot + 1`
      if [ $tot -ge 300 ] ; then
        exit;
      fi

    else 
      # create sge files
      sge_par_o="--sge \"-o ${sgeopath}\""
      sge_par_e="--sge \"-e ${sgeopath}\""
      sge_par_l="--sge \"-l mem_free=10G \""
      sge_par_q="--sge \"-q half \""
      sgerunfile="${sgeopath}/subvolReconFlexible_${batchCoreName}_sge.sh"

      cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${batchfile} > ${sgerunfile}"
      
      # echo $cmd
      eval $cmd
      chmod a+x ${sgerunfile}
      # ${sgerunfile} % if want to run dirrectly

      # run training
      sgecmd="qsub ${sgerunfile}"
      echo -e "$sgecmd\n"
      $sgecmd
    fi

    # prepare next batch file
    batchCoreName="`date +%s`_${RANDOM}_${subjid}"
    batchfile="${sgeopath}/subvolReconFlexible_batch_${batchCoreName}.sh" 

    # sleep 1
  fi

  # sleep for a bit to give sge time to deal with the new job (?)
  # sleep 100
done 10< ${locfile}

echo "Done $0"