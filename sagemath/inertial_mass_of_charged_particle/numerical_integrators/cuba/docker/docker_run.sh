#!/bin/bash

DOCKER_IMAGE=ubuntu_cuba:16.04
PROJECT_ROOT=${PWD}/..

docker run -it --rm -v ${PROJECT_ROOT}:/tmp/build \
    --workdir=${PWD} \
    -v ${PROJECT_ROOT}:${PROJECT_ROOT} \
    -v /home/${USER}/.ssh/id_rsa:/home/${USER}/.ssh/id_rsa \
    -v /home/${USER}/.ssh/id_rsa.pub:/home/${USER}/.ssh/id_rsa.pub \
    -v /home/${USER}/.ssh/known_hosts:/home/${USER}/.ssh/known_hosts \
    --user ${USER}:${USER} ${DOCKER_IMAGE} /bin/bash


