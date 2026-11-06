#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# The CPU side of tools/colab/colab_quickstart.sh in the Debian 13 chroot of tools/colab/remote_chroot_job.sh.
#
#   quickstart_body.sh J STAMP
#
# J: /opt/job (src/: the clean source tree; wps.tar.gz: WPS's executables and tables). Work in /opt/qs (WRF's release
# source, the low-resolution geography, the CPU build, the example). Each step appends "PASS|FAIL step (seconds) note"
# to J/quickstart.txt; progress in J/quickstart.log. The GPU side (tools/colab/quickstart_gpu.sh, outside the chroot)
# waits for /opt/qs/gpu/cpu_ready (STAMP), runs the example's model step and answers with /opt/qs/gpu/gpu_done.
J=$1 STAMP=$2
Q=/opt/qs G=/opt/qs/gpu
R=$J/quickstart.txt
mkdir -p $Q
rm -rf $G $Q/ex && mkdir -p $G
touch $J/setup.out $R
L() { echo "$(date +%H:%M:%S) $*" >> $J/quickstart.log; }
step() {  # NAME COMMAND...: runs COMMAND, records PASS / FAIL with the time
  local n=$1 t0=$SECONDS rc; shift
  L "$n: $*"
  "$@" > $J/$n.log 2>&1; rc=$?
  echo "$([ $rc = 0 ] && echo PASS || echo FAIL) $n ($((SECONDS - t0)) s)$([ $rc = 0 ] || echo " exit $rc, $n.log")" >> $R
  L "$n: exit $rc"
  return $rc
}
set -x
export DEBIAN_FRONTEND=noninteractive
apt-get -qq update > $J/apt.log 2>&1 && apt-get -qq install -y curl ca-certificates libpng16-16t64 >> $J/apt.log 2>&1 ||
  echo "QS_SETUP_FAILED apt (apt.log)" >> $J/setup.out
NP=$(nproc)
L "start: $NP vCPUs"

downloads() {  # WRF v4.8.0's run/ and test/em_quarter_ss/ (no build), WPS's low-resolution geography, WPS
  mkdir -p $Q/WRF $Q/geog $Q/WPS
  curl -sfL --retry 3 https://github.com/wrf-model/WRF/releases/download/v4.8.0/v4.8.0.tar.gz -o $Q/wrf.tar.gz &&
    tar xzf $Q/wrf.tar.gz -C $Q/WRF --strip-components=1 --wildcards '*/run/*' '*/test/em_quarter_ss/*' &&
    rm -f $Q/wrf.tar.gz && [ -f $Q/WRF/run/LANDUSE.TBL ] && [ -f $Q/WRF/test/em_quarter_ss/input_sounding ] || return 1
  curl -sfL --retry 3 https://www2.mmm.ucar.edu/wrf/src/wps_files/geog_low_res_mandatory.tar.gz | tar xzf - -C $Q/geog &&
    [ -d $Q/geog/WPS_GEOG_LOW_RES ] || return 1
  tar xzf $J/wps.tar.gz -C $Q/WPS
}
step downloads downloads &
DL=$!

# README.md, "Build": the CPU FP64 build
cd $J/src || exit 1
step build_cpu bash -c "cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_DOUBLE=ON &&
  cmake --build build-cpu -j $NP"
wait $DL
# USER_GUIDE.md, "Testing": ctest and the translator's tests
step ctest_cpu ctest --test-dir build-cpu --output-on-failure --timeout 1800
step translator_tests python3 -m unittest tests.tools.test_wrf2cppwrf
# USER_GUIDE.md, "Quick start": the analytic Straka case and WRF's em_quarter_ss
step straka_cpu env OMP_NUM_THREADS=$NP build-cpu/cppwrf configs/straka.ini --set output.dir=$Q/straka_cpu
step quarter_ss_cpu env OMP_NUM_THREADS=$NP build-cpu/cppwrf configs/wrf_quarter_ss.ini \
  --set case.sounding_file=$Q/WRF/test/em_quarter_ss/input_sounding --set output.dir=$Q/qss_cpu
# USER_GUIDE.md, "End-to-end example": GFS -> WPS -> cppwrf_real -> cppwrf -> wrfout with the CPU build
step example_cpu env OMP_NUM_THREADS=$NP WPS_DIR=$Q/WPS WRF_DIR=$Q/WRF GEOG_DIR=$Q/geog/WPS_GEOG_LOW_RES \
  CPPWRF_BUILD=build-cpu examples/quickstart/run.sh $Q/ex
for f in wps/geogrid.log wps/ungrib.log wps/metgrid.log run/real.log run/translate.log run/cppwrf.log run/run.ini; do
  [ -f $Q/ex/$f ] && cp $Q/ex/$f $J/ex_$(basename $f).txt
done
ls -lR $Q/ex/out >> $J/quickstart.log 2>&1
grep -E "cppwrf: done|stepping" $Q/ex/run/cppwrf.log >> $J/quickstart.log 2>&1

# the GPU side: the same run.ini, outputs in /opt/qs/gpu/out_*
echo "$STAMP" > $G/cpu_ready
end=$((SECONDS + 3600))
until [ -e $G/gpu_done ] || [ $SECONDS -ge $end ]; do sleep 10; done
cat $G/gpu.txt >> $R 2> /dev/null || echo "FAIL gpu (no gpu.txt)" >> $R
cp $G/*.log $J/ 2> /dev/null
compare_out() {  # NAME DIR: every wrfout file of DIR against the CPU run (GPU FP64: tolerance 0; FP32: relative L2)
  python3 - "$1" "$Q/ex/out" "$2" <<'EOF'
import glob, os, sys
import numpy as np, netCDF4
sys.path.insert(0, "tests/verify")
import exact
name, ref, run = sys.argv[1:]
files = sorted(os.path.relpath(f, ref) for f in glob.glob(ref + "/**/wrfout_d0*", recursive=True))
bad, worst = [], {}
for f in files:
    if not os.path.exists(os.path.join(run, f)):
        bad.append(f"{f}: missing"); continue
    a, b = exact.raw(netCDF4.Dataset(os.path.join(ref, f))), exact.raw(netCDF4.Dataset(os.path.join(run, f)))
    for v in a.variables:
        if v not in b.variables:
            bad.append(f"{f} {v}: missing"); continue
        x, y = a[v][:], b[v][:]
        if x.dtype.kind not in "fiu":
            continue
        if name == "fp64":
            d = exact.describe(x, y)
            if d: bad.append(f"{f} {v}: {d}")
        elif v in ("U", "V", "W", "T", "QVAPOR", "PH", "MU", "PSFC", "T2", "RAINNC"):
            x, y = np.asarray(x, np.float64), np.asarray(y, np.float64)
            n = np.linalg.norm(x)
            worst[v] = max(worst.get(v, 0.0), float(np.linalg.norm(x - y) / n) if n else 0.0)
print(f"{len(files)} files")
if name == "fp64":
    print("identical" if not bad else "DIFFER\n" + "\n".join(bad[:40]))
    sys.exit(1 if bad or not files else 0)
print("relative L2 (largest over the files): " + " ".join(f"{v} {e:.2e}" for v, e in sorted(worst.items())))
sys.exit(1 if bad or not files else 0)
EOF
}
step compare_gpu_fp64_vs_cpu compare_out fp64 $G/out_cuda_dp
step compare_gpu_fp32_vs_cpu compare_out fp32 $G/out_cuda
L "end"
echo END >> $R
exit 0
