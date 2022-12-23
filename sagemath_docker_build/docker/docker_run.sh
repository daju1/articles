#!/bin/sh

set -e
set -x

DISPLAY=$(env | grep DISPLAY= | sed 's/DISPLAY=//')
export DISPLAY=$DISPLAY
PROJECT_ROOT=$(dirname ${PWD})
DOCKER_IMAGE=ubuntu_sage_build:18.04
PROJECTS_DIR=$(dirname $(dirname ${PWD}))
USR3_DIR=$(dirname $(dirname $(dirname ${PWD})))

docker run -it --rm --name sage_build_container --cap-add=NET_ADMIN --device /dev/net/tun \
    --workdir=${PWD} \
    -p 8888:8888 \
    -v /home/${USER}/.local:/home/${USER}/.local \
    -v ${USR3_DIR}/winlibghemical:${PROJECT_ROOT}/winlibghemical \
    -v ${USR3_DIR}/moldyn:${PROJECT_ROOT}/moldyn \
    -v ${USR3_DIR}/science:${PROJECT_ROOT}/science \
    -v ${PROJECTS_DIR}:${PROJECT_ROOT}/work \
    -v ${PROJECT_ROOT}:${PROJECT_ROOT} \
    -v /opt/output/target/root:/opt/output/target/root \
    -v /home/${USER}/.ssh_daju/id_rsa:/home/${USER}/.ssh/id_rsa \
    -v /home/${USER}/.ssh_daju/id_rsa.pub:/home/${USER}/.ssh/id_rsa.pub \
    -v /home/${USER}/.ssh/known_hosts:/home/${USER}/.ssh/known_hosts \
    -e DISPLAY=$DISPLAY --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" \
    --user ${USER}:${USER} ${DOCKER_IMAGE} /bin/bash


