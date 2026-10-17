#!/bin/bash

# Container engine: podman (default, free) or docker (requires license via manager)
CONTAINER_ENGINE=${CONTAINER_ENGINE:-podman}

if [ -z "${1}" ] || [ -z "${2}" ]
then
    PROJECT_PATH=$(dirname ${0})
    DOCKERFILE=${PROJECT_PATH}/cuba.dockerfile
    echo "------------------------------------------------------"
    echo "no arguments are given! Using the default settings ..."
    echo "$(basename ${0}) ${PROJECT_PATH} ${DOCKERFILE}"
    echo "------------------------------------------------------"
else
    PROJECT_PATH=${1}
    DOCKERFILE=${1}/${2}
fi

DOCKER_IMAGE=ubuntu_cuba:16.04

${CONTAINER_ENGINE} image build \
    --build-arg USER_NAME=${USER} \
    --build-arg GROUP_NAME=${USER} \
    --build-arg USER_ID=$(id -u ${USER}) \
    --build-arg GROUP_ID=$(id -g ${USER}) \
    ${PROJECT_PATH} -f ${DOCKERFILE} -t ${DOCKER_IMAGE}
