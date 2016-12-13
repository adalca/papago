#!/bin/bash
# sgeSubvol2itersLSwgmm
#
# examples
# ./sgeSubvol2itersLSwgmm.sh ADNI_T1_baselines mar12_2016 Ds5Us5Reg wholevol LS_K5_dec02

if [ "$#" -lt 5 ] ; then
  echo "Usage: $0 dataName subvolVer mod proctype <dsFact>" >&2
  exit 1
fi

# this version's info
dataName=$1 # ADNI_T1_baselines or stroke or buckner
subvolVer=$2 # e.g. mar12_2016
mod=$3 # e.g. Ds5Us5Reg. Mask will be added for the weight.
procType=$4 # wholevol brain_pad10 brain_pad30
suffix="$5" # e.g. LM_K5_nov29

# set some settings based on datatype
if [ "$dataName" = "stroke" ] ; then
  ininame="subvol2itersLSwgmm_stroke.ini"
  logfilePost="selidx2loc_ds7us5_top343.txt"
  logfilePost="selidx2loc_ds7us5_rest343.txt"
else 
  ininame="subvol2itersLSwgmm.ini"
  logfilePost="selidx2loc_top343.txt"
  logfilePost="selidx2loc_rest343.txt"
fi

# cluster parameters
clustertype='sge'
clustertype='vision'
nVisionCmdQueues=10

# possible vision machines
possMachines=("asia" "africa" "america" "europe" "antarctica" "australia" "monday" "tuesday" "wednesday" "thursday" "friday" "saturday" "quickstep" "tango" "waltz" "swing" "rumba" "vision01" "vision02" "vision03" "vision04" "vision05" "vision06" "vision07" "vision08" "vision09" "vision10" "vision11" "vision12" "vision13" "vision14" "vision15" "vision16" "vision17" "vision18" "vision19" "vision20" "vision21" "vision22" "vision23" "vision24" "vision25" "vision26" "vision27" "vision28" "vision29" "vision30" "vision31" "vision32" "vision33" "vision34" "vision35" "vision36" "vision37" "vision38")


###############################################################################
# Paths
###############################################################################

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/
mcr=/data/vision/polina/shared_software/MCR/v91/

# project paths
OUTPUT_PATH="/data/vision/polina/projects/patchSynthesis/data/${dataName}/subvols/${procType}/${subvolVer}/";
SUBVOLFILES_PATH="${OUTPUT_PATH}subvols/"
RECONFILES_PATH="${OUTPUT_PATH}wgmm_${suffix}/"
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git-papago/"
CLUST_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/MCC/";

# command shell file
mccSh="${CLUST_PATH}MCC_mccSubvol2itersLSwgmm/run_mccSubvol2itersLSwgmm.sh"

# setting files
locfile="${OUTPUT_PATH}/${logfilePost}"
oinifilename="${PROJECT_PATH}/ini/${ininame}"
inifilename="${RECONFILES_PATH}/${ininame}"
if [ ! -f ${inifilename} ] ; then
    mkdir -p ${RECONFILES_PATH}
    cp ${oinifilename} ${inifilename}
else
    echo "skipping copying of ini file. Already found ${inifilename}"
fi  

###############################################################################
# Running Code
###############################################################################

# execute
tot="1"
cpui="1"
stackcmd=""
stacki="0"
while read -u10 line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  dsSubvolMat="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wtSubvolMat="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}Mask_subvol${subvolInd}.mat"
  clusterIdxMat="${RECONFILES_PATH}clusterIdx_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wgmmMat="${RECONFILES_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"

  if [ -f ${wgmmMat} ] ; then
    # echo "Skipping $subvolInd since `basename $wgmmMat` exists."
    continue;
  fi

  # subvol2itersLSwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
  lcmd="${mccSh} $mcr $dsSubvolMat $wtSubvolMat $clusterIdxMat $wgmmMat $inifilename"

  # create sge file
  sgeopath="${RECONFILES_PATH}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=10G \""
  sge_par_q="--sge \"-q half \""
  cmdrunfile="${sgeopath}/subvol2itersLSwgmm_${subvolInd}.sh"
  sgerunfile="${sgeopath}/subvol2itersLSwgmm_${subvolInd}_sge.sh"
  echo ${lcmd} > ${cmdrunfile};
  chmod a+x ${cmdrunfile}
  
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l $sge_par_q ${cmdrunfile} > ${sgerunfile}"
  # echo $cmd
  eval $cmd
  chmod a+x ${sgerunfile}

  # run sge
  if [ "$clustertype" = 'sge' ] ; then
    sgecmd="qsub ${sgerunfile}"
    echo $sgecmd
    $sgecmd
  fi

  if [ "$clustertype" = 'vision' ] ; then
    stackcmd="${stackcmd} ${sgerunfile}; "
    
    stacki=`expr $stacki + 1`
    m=`expr $stacki % $nVisionCmdQueues`
    if [ $m -eq "0" ] ; then # if we've reached our limit of jobs in queue for this machine, launch them'
      stackcmdrunfile="${sgeopath}/subvol2itersLSwgmm_stack_end${subvolInd}.sh"
      echo ${stackcmd} > ${stackcmdrunfile}
      chmod a+x ${stackcmdrunfile}

      # run on vision machine
      # cpui=`expr $cpui % 35` # if it's 34 then it shoudl stay 34 and call vision35. If it's 35 it should call vision01
      cpui=`expr $cpui + 1`
      cpui=`expr $cpui % ${#possMachines[@]}`
      machine="${possMachines[$cpui]}"
      # machine="vision`echo $cpui | awk '{printf "%02d", $0}'`"
      cmd="ssh adalca@${machine} \"nohup ${stackcmdrunfile} > ${stackcmdrunfile}.out 2> ${stackcmdrunfile}.err &\" "
      echo $cmd
      eval $cmd

      # clean up stack commands
      stackcmd=""
      stacki="0"
      tot=`expr $tot + 1`
    fi
  fi

  if [ $tot -ge 300 ] ; then
    exit;
  fi

  # sleep for a bit to give sge time to deal with the new job (?)
  # exit
  # sleep 1
done 10< ${locfile}

echo "done"