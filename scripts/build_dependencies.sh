#!/usr/bin/env bash
set -euo pipefail
export PATH="/usr/sw/bin:$PATH"
root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
build=${1:-"$root/build/dependencies"}
mkdir -p "$build"
build=$(cd "$build" && pwd)
jobs=${JOBS:-8}
cmake -S "$root/third_party/libdeflate-1.20" -B "$build/libdeflate" \
  -DCMAKE_C_COMPILER="$(command -v swgcc)" -DCMAKE_BUILD_TYPE=Release \
  -DLIBDEFLATE_BUILD_SHARED_LIB=OFF -DLIBDEFLATE_BUILD_GZIP=OFF \
  -DLIBDEFLATE_BUILD_TESTS=OFF
cmake --build "$build/libdeflate" -j "$jobs"
mkdir -p "$build/htslib"
cd "$build/htslib"
CC=swgcc CXX=swg++ CPPFLAGS="-I$root/third_party/libdeflate-1.20" \
  LDFLAGS="-L$build/libdeflate" \
  "$root/third_party/htslib-1.20/configure" --host=x86 \
  --disable-bz2 --disable-lzma --disable-libcurl --disable-s3 --disable-gcs
make -j "$jobs" libhts.a
printf 'Dependencies built in %s\n' "$build"
