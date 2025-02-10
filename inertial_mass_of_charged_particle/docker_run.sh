#!/bin/bash

DOCKER_IMAGE=sagemath/sagemath:latest
PROJECT_ROOT=${PWD}

DISPLAY=$(env | grep DISPLAY= | sed 's/DISPLAY=//')
export DISPLAY=$DISPLAY

docker run -it --rm -v ${PROJECT_ROOT}:/tmp/build \
    -v /home/${USER}/.ssh/id_rsa:/home/${USER}/.ssh/id_rsa \
    -v /home/${USER}/.ssh/id_rsa.pub:/home/${USER}/.ssh/id_rsa.pub \
    -v /home/${USER}/.ssh/known_hosts:/home/${USER}/.ssh/known_hosts \
    -e DISPLAY=$DISPLAY --volume="/tmp/.X11-unix:/tmp/.X11-unix:rw" \
    --user sage ${DOCKER_IMAGE} /bin/bash && cd /tmp/build && sage

RUN cd /tmp/build && sage
