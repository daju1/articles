#!/bin/bash

set -e
set -x

PROJECT_ROOT=${PWD}
SAGE_VER=9.1
export SAGE_ROOT=${PROJECT_ROOT}/sage-${SAGE_VER} && ${PROJECT_ROOT}/sage-${SAGE_VER}/local/bin/sage -n jupyter  --ip=0.0.0.0 --allow-root

