#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# GPU build of a checkout: $1 checkout, $2 ON|OFF (CPPWRF_DOUBLE), $3 build dir, rest: extra cmake args
# (the no-FMA FP32 realization: -DCMAKE_CUDA_FLAGS=--fmad=false). Host compiler g++-14 (CUDA 12.9 rejects g++ 15).
C=$1; P=$2; B=$3; shift 3
NVCC=${NVCC:-$HOME/cuda129/bin/nvcc}; export PATH=$(dirname $NVCC):$PATH
cd $C && rm -rf $B
cmake -S . -B $B -G Ninja -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=$P -DCPPWRF_ENABLE_CUDA=ON \
  -DCMAKE_CUDA_COMPILER=$NVCC -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++-14 -DCMAKE_CXX_COMPILER=/usr/bin/g++-14 \
  -DCMAKE_CUDA_ARCHITECTURES=${ARCH:-89} "$@" > $B.cfg.log 2>&1 || { tail -30 $B.cfg.log; exit 1; }
ninja -C $B -j ${JOBS:-8} cppwrf > $B.log 2>&1; rc=$?; grep -E 'error' $B.log | sort -u | head -20; echo "$C $B rc=$rc"
