#!/bin/bash

LAPACK_DIR=submodules/lapack

pushd ${LAPACK_DIR}

# Use default config (for gfortran)
cp make.inc.example make.inc

# Bump up -O2 to -O3
sed -i "s/-O2/-O3/g" make.inc

# Only these libs are required
make blaslib
make lib

# From LAPACK dir
popd

