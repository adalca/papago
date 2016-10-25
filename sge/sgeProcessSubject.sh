#!/bin/bash
# prepare useful volumes for a subject in the data set
#
# >$ sgeProcessSubject.sh datapath mdpath dsRate intensityNorm atlModsPath
#
# examples:
# ./sgeProcessSubject.sh \
#   /data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/proc/wholevol \
#   /data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_03_16.mat \
#   7 \
#   255 \
#   /data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/atlases/wholevol/atlMods.mat \
#   atlas \
#   0

###############################################################################
# Settings
###############################################################################

subjpath=$1
mdpath=$2
dsRate=$3
intensityNorm=$4
atlModsPath=$5
regtype=$6
padval=$7
echo $padval
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
mccSh="${PROJECT_PATH}/../MCC/MCC_processSubject/run_processSubject.sh"


###############################################################################
# Running Code
###############################################################################



for subjid in `seq 1 $nSubjects`
do
  # input files. mdpath, daRate and subjid are passed in

  # prepare command
  # processSubject(md, 1, 5, intensityNorm, BUCKNER_ATLAS_MODS);
  lcmd="${mccSh} $mcr $mdpath $subjid $dsRate $intensityNorm $atlModsPath $regtype $padval"

  # create sge files
  sgeopath="${subjpath}/sge/"
  mkdir -p ${sgeopath}
  sge_par_o="--sge \"-o ${sgeopath}\""
  sge_par_e="--sge \"-e ${sgeopath}\""
  sge_par_l="--sge \"-l mem_free=10G \""
  sgerunfile="${sgeopath}/processSubject_${subjid}.sh"

  # prepare sge run file
  cmd="${PROJECT_PATH}sge/qsub-run -c $sge_par_o $sge_par_e $sge_par_l ${lcmd} > ${sgerunfile}"
  # echo $cmd
  eval $cmd

  # run training
  sgecmd="qsub ${sgerunfile}"
  chmod a+x ${sgerunfile}
  echo -e "$sgecmd\n"
  $sgecmd
  # sleep 1

done
