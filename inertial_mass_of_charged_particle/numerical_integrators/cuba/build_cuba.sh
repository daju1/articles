#!/bin/bash
set -e

# ---------------------------------------------------------------------------
# Build the Cuba 4.2 library into ./local
#
# Idempotent: the tarball is downloaded only if it is not already present, so
# if Cuba-4.2.tar.gz already sits next to this script (e.g. copied in by hand
# to avoid a corporate TLS proxy), no network access happens at all.
#
# Corporate TLS proxy (Zscaler): set ZSCALER_CA to the path of the Zscaler
# root certificate and wget will verify against it, e.g.
#   ZSCALER_CA=$HOME/Downloads/ZscalerRootCertificate-2048-SHA256.crt ./build_cuba.sh
# As a last-resort override you can skip verification with
#   WGET_OPTS=--no-check-certificate ./build_cuba.sh
# ---------------------------------------------------------------------------

CUBA_VERSION=${CUBA_VERSION:-4.2}
CUBA_TARBALL=Cuba-${CUBA_VERSION}.tar.gz
CUBA_URL=${CUBA_URL:-https://feynarts.de/cuba/${CUBA_TARBALL}}

WGET_OPTS=${WGET_OPTS:-}
if [ -n "${ZSCALER_CA}" ]; then
    WGET_OPTS="--ca-certificate=${ZSCALER_CA} ${WGET_OPTS}"
fi

# 1. fetch the sources only if we don't already have them
if [ -f "${CUBA_TARBALL}" ]; then
    echo "==> ${CUBA_TARBALL} already present, skipping download"
else
    echo "==> downloading ${CUBA_URL}"
    wget ${WGET_OPTS} "${CUBA_URL}"
fi

# 2. unpack (only if not already unpacked)
if [ -d "Cuba-${CUBA_VERSION}" ]; then
    echo "==> Cuba-${CUBA_VERSION}/ already unpacked, skipping tar"
else
    tar zxvf "${CUBA_TARBALL}"
fi

# 3. configure, build and install into ./local
cd "Cuba-${CUBA_VERSION}"
./configure CFLAGS="-fPIC -O3" CPPFLAGS="-fPIC -O3" --prefix="$PWD/../local"
make
make install
