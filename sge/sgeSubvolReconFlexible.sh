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
  echo "Usage  : $0 dataName subvolVer subjid procType <modelType>" >&2
  echo "Example: $0 ADNI_T1_baselines mar12_2016 101528 wholevol latentSubspaceR" >&2
  exit 1
fi

dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
subjid=$3 # e.g. 101411
procType=$4 # wholevol brain_pad10 brain_pad30
suffix="$5" # e.g. iterppca_K5_ppcaK10_apr21
reconModel="$6" # 'latentMissingR' or 'latentSubspaceR'

hackyPostWgmm="_reducedSize"; # hack for reduced size wgmms... use "_reducedSize" or ""
hackyPostWgmm=""

# set some settings based on datatype
if [ "$dataName" = "stroke" ] ; then
  dsRate="7"
  usRate="5"
  locfilePost="selidx2loc_ds7us5_top343.txt"
  locfilePost="selidx2loc_ds7us5_v2.txt"
  locfilePost="selidx2loc_ds7us5_v2_rest343.txt"
  logfilePost="selidx2loc_ds7us5_special343.txt"
else 
  dsRate="5"
  usRate="5"
  # locfilePost="grididx2loc.txt"
  # 
  # locfilePost="selidx2loc_top343.txt" 
  # locfilePost="selidx2loc.txt"
  locfilePost="selidx2loc_rest343.txt"
  # locfilePost="selidx2loc_special343.txt"
fi

# cluster parameters
clustertype='vision'
clustertype='sge'
nBatchStack=100

# mod name
mod="Ds${dsRate}Us${usRate}Reg"

# possible vision machines
possMachines=("vision14" "vision15" "vision16" "vision17" "vision18" "vision19" "vision20" "vision21" "vision22" "vision23" "vision24" "vision25" "vision26" "vision27" "vision28" "vision29" "vision30" "vision31" "vision32" "vision33" "vision34" "vision35" "vision36" "vision37" "vision38" "asia" "africa" "america" "europe" "antarctica" "australia" "vision01" "vision02" "vision03" "vision04" "vision05" "vision06" "vision07" "vision08" "vision09" "vision10" "vision11" "vision12" "vision13") 
# "monday" "tuesday" "wednesday" "thursday" "friday" "saturday" "quickstep" "tango" "waltz" "swing" "rumba"
# possMachines=("vision02" "vision04" "vision05" "vision06" "vision07" "vision08" "vision09" "vision10" "vision11" "vision12" "vision13" "vision14" "vision15" "vision16" "vision17" "vision18" "vision19" "vision20" "vision21" "vision22" "vision23" "vision24" "vision25" "vision26" "vision27" "vision28" "vision29" "vision30" "vision32" "vision33" "vision34" "vision35" "vision37" "vision38" "asia" "africa" "america" "europe" "antarctica" "australia" )
# possMachines=("vision10" "vision25" "vision27" "vision38" "vision07" "vision11" "vision01")
# possMachines=("vision01" "vision02" "vision03" "vision04" "vision05" "vision06" "vision07" "vision10" "vision11" "vision12" "vision13" "vision14" "vision16" "vision17" "vision18" "vision19" "vision21" "vision22" "vision25" "vision27" "vision28" "vision29" "vision32" "vision35" "vision37" "vision38")
totlimit=250

# possMachines=("sorbet" "redbean")
# possMachines=("redbean")
# totlimit=20

###############################################################################
# Main Paths
###############################################################################

DATA_PATH="/data/vision/polina/projects/patchSynthesis/data/"
OUTPUT_PATH="${DATA_PATH}${dataName}/subvols/${procType}/${subvolVer}/"
SUBVOLS_PATH="${OUTPUT_PATH}/subvols/"
RECONS_PATH="${OUTPUT_PATH}/recons_${suffix}_${reconModel}/"
TRAIN_PATH="${OUTPUT_PATH}/wgmm_${suffix}${hackyPostWgmm}/"
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
mccSh="${PROJECT_PATH}/../MCC/MCC_mccSubvolReconFlexibleWrap/run_mccSubvolReconFlexibleWrap.sh"

###############################################################################
# Shared Inputs
###############################################################################

oiniReconFile="${PROJECT_PATH}/ini/reconFlexible_${reconModel}.ini" # default recon ini file
iniReconFile="${RECONS_PATH}/reconFlexible_${reconModel}.ini" # default recon ini file
mkdir -p ${RECONS_PATH}
cp ${oiniReconFile} ${iniReconFile}

# subject files
subjid=`basename ${SUBJ_PATH}`;
dsSubjFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${usRate}.nii.gz"
dsSubjWeightFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${usRate}_dsmask.nii.gz"
subjCorrFile="${SUBJ_PATH}/${subjid}_cor_2_ds${dsRate}_us${usRate}_size.mat"

if [ ! -f ${subjCorrFile} ] ; then
  printf "\nWarning Subject corr file ${subjCorrFile} doesn't exist, we'll use tform\n"
  subjCorrFile="${SUBJ_PATH}/${subjid}_ds${dsRate}_us${dsRate}_reg.mat"
fi

# set up locations files
locfile="${OUTPUT_PATH}/${locfilePost}" 
printf "locfile: ${locfile} \n\n"

# makse sge path
sgeopath="${RECONS_PATH}/sge/"
mkdir -p ${sgeopath}

# subject files
subjfolder="${RECONS_PATH}/${subjid}/"
mkdir -p ${subjfolder}

###############################################################################
# Running Code
###############################################################################

# prepare batch files
stacki="0"
batchCoreName="`date +%s`_${RANDOM}_${subjid}"
batchfile="${sgeopath}/subvolReconFlexible_batch_${batchCoreName}.sh"
batchVolFile="${sgeopath}/subvolReconFlexible_batch_${batchCoreName}_subvols.txt"

# prepare machine execution
tot="1"
cpui="1"
tot="1"
m="1"

totLines=`wc -l ${locfile} | cut -f 1 -d " "`
lineCount=0

echo "total lines: $totLines"

# run different settings
while read -u10 line
do
  lineCount=`expr $lineCount + 1`
  subvolInd=`echo ${line} | cut -d " " -f 1`
  
  skip=0
  wgmmfile="${TRAIN_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  if [ ! -f $wgmmfile ] ; then
    # printf "skipping $subvolInd since I can't find $wgmmfile \n\n"
    skip=1
  fi

  subjoutFile="${subjfolder}/subvolReconFlexible_${dataName}_${procType}_${mod}_subvol${subvolInd}_${subjid}.mat"
  if [ -f $subjoutFile ] ; then
    # printf "skipping $subvolInd since $subjoutFile is present \n\n"
    skip=1
  fi

  # continue if skippping and we're NOT at the last line (otherwise we need to execute what we have)
  if [ $skip -eq 1 ] && [ $lineCount -ne $totLines ] ; then continue; fi

  # as long as we'te not skipping, add this subvol to the list
  if [ $skip -eq 0 ]  ; then 
    
    # add subvolume to batch 
    printf "${subvolInd}\n" >> ${batchVolFile} 

    # increment the queue
    stacki=`expr $stacki + 1`
    m=`expr $stacki % $nBatchStack`
  fi
  
  
  if [ $m -eq "0" ] || [ $lineCount -eq $totLines ] ; then # ready to run
    
    # prepare matlab command
    wgmmfile="${TRAIN_PATH}wgmm_${dataName}_${procType}_${mod}_subvol%%d.mat"
    subjoutFile="${subjfolder}/subvolReconFlexible_${dataName}_${procType}_${mod}_subvol%%d_${subjid}.mat"
    subvolfile="${SUBVOLS_PATH}${dataName}_${procType}_${mod}_subvol%%d.mat"
    lcmd="${mccSh} $mcr $wgmmfile $subvolfile $iniReconFile $dsSubjFile $dsSubjWeightFile $subjCorrFile $subjoutFile ${batchVolFile};"
    printf "$lcmd" > ${batchfile} 
    chmod a+x ${batchfile}


    if [ "$clustertype" = 'vision' ] ; then

      # run on vision machine
      # cpui=`expr $cpui % 35` # if it's 34 then it shoudl stay 34 and call vision35. If it's 35 it should call vision01
      cpui=`expr $cpui + 1`
      cpui=`expr $cpui % ${#possMachines[@]}`
      machine="${possMachines[$cpui]}"

      rand=`expr $RANDOM % ${#possMachines[@]}`
      machine="${possMachines[$rand]}"
      echo "Warning: randomly assigned machine $machine"

      cmd="ssh ${machine} \"nohup ${batchfile} > ${batchfile}.out 2> ${batchfile}.err &\" " # took out adalca@${machine}
      echo $cmd
      eval $cmd &

      tot=`expr $tot + 1`
      if [ $tot -ge $totlimit ] ; then
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
    batchVolFile="${sgeopath}/subvolReconFlexible_batch_${batchCoreName}_subvols.txt"

    # echo "TRYME"
    # sleep 100
  fi

  # echo "$lineCount/$totLines (m: $m)"

done 10< ${locfile}

echo "Done $0"

