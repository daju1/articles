#!/bin/sh

set -e
set -x

DISPLAY=$(env | grep DISPLAY= | sed 's/DISPLAY=//')
export DISPLAY=$DISPLAY
PROJECT_ROOT=$(dirname ${PWD})
DOCKER_IMAGE=ubuntu_sage_build:18.04
PROJECTS_DIR=$(dirname $(dirname ${PWD}))
USR3_DIR=$(dirname $(dirname $(dirname ${PWD})))

# docker network create sagemath

docker run -it --rm --name sage_build_container \
    --cap-add=SYS_PTRACE \
    --security-opt seccomp=unconfined \
    --device /dev/net/tun \
    --workdir=${PROJECT_ROOT} \
    -p 8880:8888 \
    -e HOME=/home/${USER} \
    -v ${HOME}/.local:/home/${USER}/.local \
    -v ${HOME}/.ssh:/home/${USER}/.ssh \
    -v ${USR3_DIR}/winlibghemical:${PROJECT_ROOT}/winlibghemical \
    -v ${USR3_DIR}/moldyn:${PROJECT_ROOT}/moldyn \
    -v ${USR3_DIR}/science:${PROJECT_ROOT}/science \
    -v ${USR3_DIR}/study:${PROJECT_ROOT}/study \
    -v ${USR3_DIR}/tss_bot:${PROJECT_ROOT}/tss_bot \
    -v ${USR3_DIR}/jenkins:${PROJECT_ROOT}/jenkins \
    -v ${USR3_DIR}/test:${PROJECT_ROOT}/test \
    -v ${USR3_DIR}/LaTex2Docx:${PROJECT_ROOT}/LaTex2Docx \
    -v ${USR3_DIR}/Prompt-Engineering-Guide:${PROJECT_ROOT}/Prompt-Engineering-Guide \
    -v ${PROJECTS_DIR}:${PROJECT_ROOT}/work \
    -v ${PROJECT_ROOT}:${PROJECT_ROOT} \
    -e DISPLAY=$DISPLAY --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" \
    --user ${USER}:${USER} ${DOCKER_IMAGE} /bin/bash


