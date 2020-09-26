#!/bin/bash

set -e
set -x

export SAGE_ROOT=${PROJECT_ROOT}/sage-9.1 && ${PROJECT_ROOT}/sage-9.1/local/bin/sage -n jupyter  --ip=0.0.0.0 --allow-root

