#!/bin/bash

set -e
set -x

PROJECT_ROOT=${PWD}
SAGE_VER=
SAGE_ROOT="${PROJECT_ROOT}/sage${SAGE_VER}"
SAGE_LOCAL="$SAGE_ROOT/local"
PYTHON_EXEC="$SAGE_LOCAL/bin/python3"

if [ -f "$SAGE_LOCAL/bin/jupyter-notebook" ]; then
    JUPYTER_EXEC="$SAGE_LOCAL/bin/jupyter-notebook"
    JUPYTER_ARGS="--ip=0.0.0.0 --allow-root"
elif [ -f "$SAGE_LOCAL/bin/jupyter" ]; then
    JUPYTER_EXEC="$SAGE_LOCAL/bin/jupyter"
    JUPYTER_ARGS="notebook --ip=0.0.0.0 --allow-root"
else
    echo "ERROR: jupyter or jupyter-notebook not found in $SAGE_LOCAL/bin"
    exit 1
fi

if [ ! -x "$PYTHON_EXEC" ]; then
    echo "Error: Python Snot found $PYTHON_EXEC"
    exit 1
fi

if [ ! -f "$JUPYTER_EXEC" ]; then
    echo "Jupiter run script not found in $JUPYTER_EXEC"
    ls -l "$SAGE_LOCAL/bin/" | grep sage || true
    exit 1
fi


export SAGE_ROOT="$SAGE_ROOT"
export SAGE_LOCAL="$SAGE_LOCAL"
export PATH="$SAGE_LOCAL/bin:$PATH"

SAGE_DIR=${SAGE_ROOT}/local/bin

if [ ! -x "./sage" ]; then
    echo "./sage not found and not executable"
    ls -l ./sage || true
    exit 1
fi

gdb --args "$PYTHON_EXEC" "$JUPYTER_EXEC" $JUPYTER_ARGS

# sage/local/bin/jupyter kernelspec list
# Available kernels:
#   root        /home/.../.local/share/jupyter/kernels/root
#   python3     /.../articles/sagemath_docker_build/sage/local/share/jupyter/kernels/python3
#   sagemath    /.../articles/sagemath_docker_build/sage/local/share/jupyter/kernels/sagemath


