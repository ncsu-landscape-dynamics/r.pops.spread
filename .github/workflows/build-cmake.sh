#!/usr/bin/env bash

# Build and install GRASS GIS using its CMake build system.
# This mirrors build.sh (which uses the Autotools build) so that the CMake
# code path of g.extension (and thus the addon's CMakeLists.txt) is exercised.

# The build step requires something like:
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
INSTALL_PREFIX="$2"
BRANCH="$3"

cd "$WORKDIR"

# GRASS source

git clone https://github.com/OSGeo/grass.git --branch "$BRANCH" --depth=1

cd grass

# GRASS build with CMake. The set of enabled features mirrors build.sh.
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
    -DWITH_LARGEFILES=ON \
    -DWITH_ZSTD=ON \
    -DWITH_BZLIB=ON \
    -DWITH_READLINE=ON \
    -DWITH_OPENMP=ON \
    -DWITH_TIFF=ON \
    -DWITH_FREETYPE=ON \
    -DWITH_GEOS=ON \
    -DWITH_SQLITE=ON \
    -DWITH_FFTW=ON \
    -DWITH_NETCDF=ON \
    -DWITH_OPENGL=OFF \
    -DWITH_PDAL=OFF \
    -DWITH_GUI=OFF \
    -DWITH_DOCS=OFF \
    -DWITH_NLS=OFF

cmake --build build
cmake --install build

# Delete the source code.
# The source code should not be needed anymore and
# it may cause problems when it is in the same directory used further
# for the compilation of the tool itself
# (e.g., g.extension will recurse to it).
cd ..
rm -rf grass
