#!/usr/bin/env bash
# mcc files relevant to subspace restoration project.
#
# usage:
# mcc relpath1 relpath2
#
# where relpath are paths relative to ${PROJECT_PATH}
# Example:
# >$ ./mcc mcc/mccTrain mcc/mccPrepareSubvols

# prepare project and toolbox paths
# if want to make current: "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MAIN_PATH="/data/vision/polina/users/adalca/patchSynthesis/subspace"
PROJECT_PATH="${MAIN_PATH}/git-papago"
TOOLBOX_PATH="/data/vision/polina/users/adalca/MATLAB/toolboxes"
EXTTOOLBOX_PATH="/data/vision/polina/users/adalca/MATLAB/external_toolboxes"

# MCC-related paths
MCCBUILD_PATH="${TOOLBOX_PATH}/mgt/src/mcc/" # mccBuild script path
MCC_RUN_DIR="/afs/csail.mit.edu/system/common/matlab/2013b/bin/mcc"
MCC_RUN_DIR="/afs/csail.mit.edu/system/common/matlab/2016b/bin/mcc"

# need to add main path to system path (?).
export PATH="${MAIN_PATH}:$PATH"

## run mcc on desired (*.m) files.
for pfilename in "$@"
do
  filename=`basename ${pfilename}`

  # run via mccBuild.sh
  ${MCCBUILD_PATH}/mccBuild.sh \
    ${MCC_RUN_DIR} \
    ${PROJECT_PATH}/${pfilename}.m \
    ${MAIN_PATH}/MCC/MCC_${filename} \
    ${PROJECT_PATH}/src/ \
    ${PROJECT_PATH}/ext/ \
    ${TOOLBOX_PATH} \
    ${EXTTOOLBOX_PATH}
done
