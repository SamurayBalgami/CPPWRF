#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# CUDA 12.9 with glibc >= 2.41 (Ubuntu 26.04: 2.43): glibc declares sinpi/cospi/rsqrt(+f) noexcept, CUDA's
# crt/math_functions.h does not, so nvcc fails ("exception specification is incompatible"). Without touching the
# system install: a user-local view ~/cuda129 (symlinks, a real copy of bin/ since nvcc finds its root from its own
# path, and a patched math_functions.h). Use -DCMAKE_CUDA_COMPILER=$HOME/cuda129/bin/nvcc.
set -e
R=${CUDA_ROOT:-/usr/local/cuda-12.9}; L=${1:-$HOME/cuda129}
rm -rf $L; mkdir -p $L
for e in $R/*; do n=$(basename $e); case $n in bin|targets) ;; *) ln -s $e $L/$n ;; esac; done
cp -a $R/bin $L/bin
mkdir -p $L/targets; cp -as $R/targets/x86_64-linux $L/targets/
H=$L/targets/x86_64-linux/include/crt/math_functions.h
rm $H; cp $R/targets/x86_64-linux/include/crt/math_functions.h $H
sed -i -E 's/^(extern __DEVICE_FUNCTIONS_DECL__ __device_builtin__ (double|float) +(sinpi|sinpif|cospi|cospif|rsqrt|rsqrtf)\((double|float) x\));/\1 noexcept (true);/' $H
grep -cE '(sinpi|sinpif|cospi|cospif|rsqrt|rsqrtf)\((double|float) x\) noexcept' $H
