#!/bin/sh

# Open an extra shell inside the already-running articles_container.
CONTAINER_ENGINE=${CONTAINER_ENGINE:-podman}

${CONTAINER_ENGINE} exec -it --workdir=$(dirname ${PWD}) articles_container /bin/bash
