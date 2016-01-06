#!/usr/bin/env bash
# mcc files relevant to subspace restoration project.

# prepare project and toolbox paths
# if want to make current: "$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MAIN_PATH="/data/vision/polina/users/adalca/patchSynthesis"
PROJECT_PATH="${MAIN_PATH}/subspace/"
TOOLBOX_PATH="${MAIN_PATH}/toolboxes"

# MCC-related paths
MCCBUILD_PATH="${TOOLBOX_PATH}/mgt/src/mcc/" # mccBuild script path
MCC_RUN_DIR="/afs/csail.mit.edu/system/common/matlab/2013b/bin/mcc"

# need to add main path to system path (?).
export PATH="${MAIN_PATH}:$PATH"

## run mcc on desired (*.m) files.
for filename in sgeTrain sgePrepareSubvols
do
  # run via mccBuild.sh
  ${MCCBUILD_PATH}/mccBuild.sh \
    ${MCC_RUN_DIR} \
    ${PROJECT_PATH}/clust/${filename}.m \
    ${MAIN_PATH}/MCC/MCC_${filename} \
    ${PROJECT_PATH}/src/ \
    ${TOOLBOX_PATH}
done
