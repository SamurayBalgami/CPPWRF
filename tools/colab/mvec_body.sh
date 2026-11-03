#!/usr/bin/env bash
# The proof of src/core/libmvec_x86.hpp (glibc libmvec's vector log / log10 with Intel's rcpps) on one machine, run by
# tools/colab/colab_mvec.sh in the Debian 13 chroot of a Colab runtime (glibc 2.41 as the GCP VMs and the references):
#
#   mvec_body.sh J N
#
# J holds src/ (git archive). Logs in J: cpu.txt; unit.log (tests/unit/test_libmvec built by CMake as in the suite:
# the Intel hashes on any CPU, on Intel also libmvec and rcpps themselves); check.log (tools/libmvec_tables.py --check:
# this CPU's rcpps against the table); regen.log (Intel only: the header made again from glibc 2.41's source tarball
# and this CPU, compared with the committed one); proof.log (tools/libmvec_proof.cpp with N random inputs per
# function: on Intel every input must equal libmvec; elsewhere it counts libmvec's differences); binary.log (the
# cppwrf binary: approximate-reciprocal instructions, the libmvec functions it calls); bug_table.log, bug_order.log
# (a deliberate fault each: one step of the rcpps table moved, log's last two sums reassociated; test_libmvec must
# fail).
J=$1 N=${2:-10000000000}
set -x
touch $J/setup.out
cd $J/src || exit 1
lscpu | grep -E '^(Vendor ID|Model name|CPU\(s\)):' > $J/cpu.txt
INTEL=0; grep -q GenuineIntel /proc/cpuinfo && INTEL=1
FL=(-std=c++20 -O3 -march=skylake -mtune=skylake -ffp-contract=off -fopenmp -DCPPWRF_DOUBLE)
cmake -S . -B $J/build -DCMAKE_BUILD_TYPE=Release -DCPPWRF_DOUBLE=ON -DCPPWRF_ENABLE_CUDA=OFF -DCPPWRF_NATIVE=OFF \
  -DCMAKE_CXX_FLAGS="-march=skylake -mtune=skylake" > $J/cmake.log 2>&1 &&
  cmake --build $J/build -j"$(nproc)" --target test_libmvec cppwrf >> $J/cmake.log 2>&1 || echo "BUILD_FAILED (cmake.log)" >> $J/setup.out
$J/build/test_libmvec > $J/unit.log 2>&1; echo "exit $?" >> $J/unit.log
python3 tools/libmvec_tables.py --check src/core/libmvec_tables.hpp > $J/check.log 2>&1; echo "exit $?" >> $J/check.log
if [ $INTEL = 1 ]; then
  python3 -c "import io, sys, tarfile, urllib.request as u
t = tarfile.open(fileobj=io.BytesIO(u.urlopen('https://ftp.gnu.org/gnu/glibc/glibc-2.41.tar.xz').read()), mode='r:xz')
t.extractall(sys.argv[1], [m for m in t.getmembers() if m.name.startswith('glibc-2.41/sysdeps/x86_64/fpu/')], filter='data')" \
    $J > $J/regen.log 2>&1
  python3 tools/libmvec_tables.py $J/glibc-2.41 > $J/regen.hpp 2>> $J/regen.log; echo "exit $?" >> $J/regen.log
  diff <(grep -v 'rcpps dumped on' src/core/libmvec_tables.hpp) <(grep -v 'rcpps dumped on' $J/regen.hpp) >> $J/regen.log
  echo "diff exit $?" >> $J/regen.log; grep 'rcpps dumped on' $J/regen.hpp >> $J/regen.log
fi
g++ "${FL[@]}" -Isrc tools/libmvec_proof.cpp -o $J/libmvec_proof -lmvec -lm > $J/proof.log 2>&1
$J/libmvec_proof "$N" >> $J/proof.log 2>&1; echo "exit $?" >> $J/proof.log
{ echo "approximate reciprocal instructions in cppwrf: $(objdump -d $J/build/cppwrf | grep -c -E '\s(v?rcp(ps|ss)|v?rsqrt(ps|ss)|vrcp14|vrsqrt14)\s')"
  echo "libmvec functions cppwrf calls:"; nm -D --undefined-only $J/build/cppwrf | grep _ZGV; } > $J/binary.log 2>&1
for bug in table order; do
  rm -rf $J/bug; mkdir -p $J/bug; cp -r src tests $J/bug/
  case $bug in
    table) sed -i '0,/^    1, 3, 5,/s//    2, 3, 5,/' $J/bug/src/core/libmvec_tables.hpp ;;
    order) sed -i 's/  const double l = d(kLogLA\[int(rcp) - 512\]) + rp;/  const double l = (d(kLogLA[int(rcp) - 512]) + r) + r2;/' \
             $J/bug/src/core/libmvec_x86.hpp ;;
  esac
  diff -r src $J/bug/src > $J/bug_$bug.log
  g++ "${FL[@]}" -I$J/bug/src $J/bug/tests/unit/test_libmvec.cpp -o $J/bug_test -lmvec -lm >> $J/bug_$bug.log 2>&1 &&
    $J/bug_test >> $J/bug_$bug.log 2>&1
  echo "exit $?" >> $J/bug_$bug.log
done
rm -rf $J/bug $J/bug_test $J/glibc-2.41
