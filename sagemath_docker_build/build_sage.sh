#!/bin/sh

set -e
set -x

SAGE_VER=9.1
SAGE_DIR=sage-${SAGE_VER}
SAGE_TAR=$SAGE_DIR.tar.gz
WORK_DIR=`pwd`

SAGE_INSTALL=/opt/sagemath-${SAGE_VER}

export SAGE_ROOT=$WORK_DIR/$SAGE_DIR
export SAGE_LOCAL=$SAGE_ROOT/local


if [ ! -f $SAGE_TAR ]; then
	wget http://mirrors.mit.edu/sage/src/$SAGE_TAR
fi

if [ ! -d $SAGE_DIR ]; then
	tar xvf $SAGE_TAR
fi

cd $SAGE_DIR

./configure  --prefix=$SAGE_LOCAL
# MAKE='make -jNUM' make
make -j$(nproc)
sudo make install
