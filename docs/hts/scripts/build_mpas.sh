#!/usr/bin/env bash
# Build unmodified MPAS-Atmosphere (and init_atmosphere) the way Phase 1 did.
#
# I/O is the bundled SMIOL library over PnetCDF. PIO is not required.
# Default real kind is single precision (MPAS default; do not pass PRECISION=single).
#
# Ubuntu 24.04 packages:
#   gfortran g++ make cmake pkg-config
#   libopenmpi-dev openmpi-bin
#   libnetcdf-dev libnetcdff-dev netcdf-bin
#   libhdf5-dev zlib1g-dev
#   libpnetcdf-dev pnetcdf-bin
#
# Distro libraries do not live in a single prefix with include/ and lib/,
# which is what the MPAS Makefile looks for. This script builds a symlink
# prefix. The PnetCDF Fortran module shipped by Ubuntu is under
# gfortran-mod-15; the default gfortran (13.3) accepts that module on this image.
#
# Usage:
#   ./docs/hts/scripts/build_mpas.sh [atmosphere|init_atmosphere|both]
#
# Environment:
#   MPAS_DEPS   symlink prefix (default: $HOME/mpas-deps)
#   MPAS_SRC    model root (default: repo root, two levels above this script)
#   MPAS_JOBS   make parallelism (default: 2)

set -euo pipefail

CORE_ARG="${1:-both}"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MPAS_SRC=${MPAS_SRC:-$(cd "$SCRIPT_DIR/../../.." && pwd)}
MPAS_DEPS=${MPAS_DEPS:-"$HOME/mpas-deps"}
MPAS_JOBS=${MPAS_JOBS:-2}

link_prefix() {
    local prefix=$1
    mkdir -p "$prefix/netcdf/include" "$prefix/netcdf/lib" \
             "$prefix/pnetcdf/include" "$prefix/pnetcdf/lib"

    local f
    for f in /usr/include/netcdf.h /usr/include/netcdf.inc /usr/include/*.mod; do
        [ -e "$f" ] || continue
        ln -sfn "$f" "$prefix/netcdf/include/$(basename "$f")"
    done
    if [ -d /usr/include/hdf5 ]; then
        ln -sfn /usr/include/hdf5 "$prefix/netcdf/include/hdf5"
    fi
    ln -sfn /usr/lib/x86_64-linux-gnu/libnetcdf.so "$prefix/netcdf/lib/libnetcdf.so"
    ln -sfn /usr/lib/x86_64-linux-gnu/libnetcdff.so "$prefix/netcdf/lib/libnetcdff.so"

    ln -sfn /usr/include/pnetcdf.h "$prefix/pnetcdf/include/pnetcdf.h"
    if [ -f /usr/include/x86_64-linux-gnu/pnetcdf.inc ]; then
        ln -sfn /usr/include/x86_64-linux-gnu/pnetcdf.inc "$prefix/pnetcdf/include/pnetcdf.inc"
    fi
    local pmod
    pmod=$(find /usr/lib/x86_64-linux-gnu/fortran -name pnetcdf.mod 2>/dev/null | head -n 1 || true)
    if [ -z "$pmod" ]; then
        echo "ERROR: pnetcdf.mod not found. Install libpnetcdf-dev." >&2
        exit 1
    fi
    ln -sfn "$pmod" "$prefix/pnetcdf/include/pnetcdf.mod"
    ln -sfn /usr/lib/x86_64-linux-gnu/libpnetcdf.so "$prefix/pnetcdf/lib/libpnetcdf.so"
}

build_core() {
    local core=$1
    echo "===== building CORE=$core ====="
    make -C "$MPAS_SRC" -j"$MPAS_JOBS" gfortran \
        "CORE=$core" \
        "NETCDF=$MPAS_DEPS/netcdf" \
        "PNETCDF=$MPAS_DEPS/pnetcdf"
}

link_prefix "$MPAS_DEPS"

case "$CORE_ARG" in
    atmosphere) build_core atmosphere ;;
    init_atmosphere) build_core init_atmosphere ;;
    both)
        build_core atmosphere
        build_core init_atmosphere
        ;;
    *)
        echo "Usage: $0 [atmosphere|init_atmosphere|both]" >&2
        exit 1
        ;;
esac

echo "Built:"
ls -l "$MPAS_SRC/atmosphere_model" "$MPAS_SRC/init_atmosphere_model" 2>/dev/null || true
