#!/usr/bin/env bash

# The make step requires something like:
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$PREFIX/lib"
# further steps additionally require:
# export PATH="$PATH:$PREFIX/bin"

# fail on non-zero return code from a subprocess
set -e

if [ -z "$3" ]
then
    >&2 echo "Usage: $0 <workdir> <prefix> <branch>"
    >&2 echo "<workdir>  Working directory (must exists)"
    >&2 echo "<prefix>   Install prefix"
    >&2 echo "<branch>   Branch name (passed to git clone --branch)"
    exit 1
fi

WORKDIR="$1"
export INSTALL_PREFIX="$2"
BRANCH="$3"

cd $WORKDIR

# GRASS GIS

git clone https://github.com/OSGeo/grass.git --branch "$BRANCH" --depth=1

cd grass

./configure \
    --prefix="$INSTALL_PREFIX/" \
    --enable-largefile \
    --with-cxx \
    --with-zstd \
    --with-bzlib \
    --with-readline \
    --with-openmp \
    --with-pthread \
    --with-tiff \
    --with-freetype \
    --with-freetype-includes="/usr/include/freetype2/" \
    --with-proj-share=/usr/share/proj \
    --with-geos \
    --with-sqlite \
    --with-fftw \
    --with-netcdf \
    --without-pdal

make
make install

# Delete the source code.
# The source code should not be needed anymore and
# it may cause problems when it is in the same directory used further
# for the compilation of the module itself
# (e.g., g.extension will recurse to it).
cd ..
rm -rf grass
