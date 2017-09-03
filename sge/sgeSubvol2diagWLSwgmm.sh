#!/bin/bash
# sgeSubvol2diagLSwgmm
#
# examples
# ./sgeSubvol2diagWLSwgmm.sh ADNI_T1_baselines aug_2017 Ds5Us5Reg wholevol LS_K5_dec02

if [ "$#" -lt 5 ] ; then
  echo "Usage: $0 dataName subvolVer mod proctype" >&2
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
  ininame="subvol2diagWLSwgmm_stroke.ini"
  logfilePost="selidx2loc_ds7us5_top343.txt"
  logfilePost="selidx2loc_ds7us5_rest343.txt"
  logfilePost="selidx2loc_ds7us5_v2.txt"
else 
  ininame="subvol2diagWLSwgmm.ini" 
  logfilePost="selidx2loc_ds5us5_gs11x11x11.txt"
fi

# cluster parameters
clustertype='vision'
clustertype='sge' 
nVisionCmdQueues=50

# possible vision machines
possMachines=("vision11" "vision12" "vision13" "vision14" "vision15" "vision16" "vision17" "vision18" "vision19" "vision20" "vision21" "vision22" "vision23" "vision24" "vision25" "vision26" "vision27" "vision28" "vision29" "vision30" "vision31" "vision32" "vision33" "vision34" "vision35" "vision36" "vision37" "vision38" "asia" "africa" "america" "europe" "antarctica" "australia" "vision01" "vision02" "vision03" "vision04" "vision05" "vision06" "vision07" "vision08" "vision09" "vision10" )
# possMachines=("asia" "africa" "america" "antarctica" "australia" "monday" "tuesday" "wednesday" "thursday" "friday" "saturday")
#  "quickstep" "tango" "waltz" "swing" "rumba"
possMachines=("monday" "tuesday" "wednesday" "thursday" "friday" "saturday")
possMachines=("vision01" "vision02" "vision03" "vision04" "vision05" "vision06" "vision07" "vision10" "vision11" "vision12" "vision13" "vision14" "vision16" "vision17" "vision18" "vision19" "vision21" "vision22" "vision25" "vision27" "vision28" "vision29" "vision32" "vision35" "vision37" "vision38")
# possMachines=("vision08" "vision09" "vision15" "vision18" "vision20" "vision22" "vision30" "vision34")
possMachines=("perilla"  "lemongrass" "quassia" "orris" "lovage" "boldo" "olida" "cassia" "jimbu" "salep" "hickory" "peppercorn" "cardamom" "fenugreek" "rosemary" "paprika")
# "athelas" "mugwort" "advieh" "sassafras"
totlimit=500


# possMachines=("sorbet" "redbean")
# possMachines=("redbean")
# totlimit=20

###############################################################################
# Paths
###############################################################################

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/
mcr=/data/vision/polina/shared_software/MCR/v91/

# project paths
OUTPUT_PATH="/data/vision/polina/projects/patchSynthesis/data/${dataName}/subvols/${procType}/${subvolVer}/";
SUBVOLFILES_PATH="${OUTPUT_PATH}/subvols/"
RECONFILES_PATH="${OUTPUT_PATH}/wgmm/${suffix}"
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git-papago/"
CLUST_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/MCC/";

# command shell file
mccSh="${CLUST_PATH}MCC_mccSubvol2diagLSwgmm/run_mccSubvol2diagWLSwgmm.sh"

# setting files
locfile="${OUTPUT_PATH}/subvol_loc_files/${logfilePost}"
printf  "$locfile \n"
oinifilename="${PROJECT_PATH}/ini/${ininame}"
inifilename="${RECONFILES_PATH}/${ininame}"
if [ ! -f ${inifilename} ] ; then
    mkdir -p ${RECONFILES_PATH}
    cp ${oinifilename} ${inifilename}
else
    echo "skipping copying of ini file. Already found ${inifilename}"
fi  


# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH="${OUTPUT_PATH}/sge/"
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}


###############################################################################
# Running Code
###############################################################################

# execute
tot="1"
cpui="1"
stackcmd=""
stacki="0"
submitted="0"
while read -u10 line
do
  subvolInd=`echo ${line} | cut -d " " -f 1`
  dsSubvolMat="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wtSubvolMat="${SUBVOLFILES_PATH}${dataName}_${procType}_${mod}Mask_subvol${subvolInd}.mat"
  clusterIdxMat="${RECONFILES_PATH}clusterIdx_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"
  wgmmMat="${RECONFILES_PATH}wgmm_${dataName}_${procType}_${mod}_subvol${subvolInd}.mat"

  if [ -f ${wgmmMat} ] ; then
    echo "Skipping $subvolInd since `basename $wgmmMat` exists."
    continue;
  fi

  # subvol2diagLSwgmm(dsSubvolMat, wtSubvolMat, clusterIdxMat, wgmmMat, iniFilename)
  lcmd="${mccSh} $mcr $dsSubvolMat $wtSubvolMat $clusterIdxMat $wgmmMat $inifilename"

  # create sge file
  sgeopath="${RECONFILES_PATH}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=10G \""
  sge_par_q="--sge \"-q half \""
  cmdrunfile="${sgeopath}/subvol2diagLSwgmm_${subvolInd}.sh"
  sgerunfile="${sgeopath}/subvol2diagLSwgmm_${subvolInd}_sge.sh"
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
    submitted=`expr $submitted + 1`
    echo $submitted
  fi

  if [ "$clustertype" = 'vision' ] ; then
    stackcmd="${stackcmd} ${sgerunfile}; "
    
    stacki=`expr $stacki + 1`
    m=`expr $stacki % $nVisionCmdQueues`
    if [ $m -eq "0" ] ; then # if we've reached our limit of jobs in queue for this machine, launch them'
      stackcmdrunfile="${sgeopath}/subvol2diagLSwgmm_stack_end${subvolInd}.sh"
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
      submitted=`expr $submitted + 1`
    fi
  fi

  if [ $tot -ge $totlimit ] ; then
    exit;
  fi

  # sleep for a bit to give sge time to deal with the new job (?)
  # exit
  # sleep 1
done 10< ${locfile}

echo "done $0; submitted $submitted jobs"