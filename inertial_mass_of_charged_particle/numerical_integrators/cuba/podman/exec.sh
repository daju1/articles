#!/bin/sh

# Open an extra shell inside the already-running cuba_container.
CONTAINER_ENGINE=${CONTAINER_ENGINE:-podman}

${CONTAINER_ENGINE} exec -it --workdir=$(dirname ${PWD}) cuba_container /bin/bash
