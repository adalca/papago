#!/bin/bash
# prepare useful volumes for a subject in the data set
#
# >$ sgeProcessSubjectStroke.sh datapath mdpath intensityNorm atlModsPath
#
# examples:
# ./sgeProcessSubjectStroke.sh \
#   /data/vision/polina/scratch/adalca/patchSynthesis/data/stroke/proc \
#   /data/vision/polina/scratch/adalca/patchSynthesis/data/stroke/md/adalca_restor_md_2016_02_27.mat \
#   400 \
#   /data/vision/polina/scratch/adalca/patchSynthesis/data/stroke/atlases/atlMods.mat

if [ "$#" -lt 4 ] ; then
  echo "Usage: $0 procpatch mdpatch normfact atlMods" >&2
  exit 1
fi

###############################################################################
# Settings
###############################################################################

subjpath=$1
mdpath=$2
intensityNorm=$3
atlModsPath=$4
regtype=$5
padval=$6
nSubjects=`ls $subjpath | wc -l`

# prepare SGE variables necessary to move SGE environment away from AFS.
export SGE_LOG_PATH=/data/vision/polina/scratch/adalca/patchSynthesis/sge/
export SGE_O_PATH=${SGE_LOG_PATH}
export SGE_O_HOME=${SGE_LOG_PATH}

# MCR file. This has to match the MCC version used in mcc.sh
mcr=/data/vision/polina/shared_software/MCR/v82/

# project paths
PROJECT_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace/git/";

# training shell file
mccSh="${PROJECT_PATH}/../MCC/MCC_processSubjectStroke/run_processSubjectStroke.sh"


###############################################################################
# Running Code
###############################################################################



for subjid in `seq 1 $nSubjects`
do
  # input files. mdpath, daRate and subjid are passed in

  # prepare command
  lcmd="${mccSh} $mcr $mdpath $subjid $intensityNorm $atlModsPath $regtype $padval"

  # create sge files
  sgeopath="${subjpath}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=100G \""
  sgerunfile="${sgeopath}/processSubjectStroke_${subjid}.sh"

  # prepare sge run file
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l ${lcmd} > ${sgerunfile}"
  echo $cmd
  eval $cmd

  # run training
  sgecmd="qsub ${sgerunfile}"
  chmod a+x ${sgerunfile}
  echo -e "$sgecmd\n"
  $sgecmd

done
