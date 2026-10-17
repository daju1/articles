#!/bin/sh

set -e
set -x

# Container engine: podman (default, free) or docker (requires license via manager)
CONTAINER_ENGINE=${CONTAINER_ENGINE:-podman}

# Mount the parent of this folder (the "cuba" project root) so build_cuba.sh,
# cuba.mk and all the calc_* sources are visible inside the container.
PROJECT_ROOT=$(dirname ${PWD})
DOCKER_IMAGE=ubuntu_cuba:16.04

${CONTAINER_ENGINE} run -it --rm --name cuba_container \
    --workdir=${PROJECT_ROOT} \
    -e HOME=/home/${USER} \
    -v ${HOME}/.local:/home/${USER}/.local \
    -v ${PROJECT_ROOT}:${PROJECT_ROOT} \
    --user $(id -u ${USER}):$(id -g ${USER}) ${DOCKER_IMAGE} /bin/bash
