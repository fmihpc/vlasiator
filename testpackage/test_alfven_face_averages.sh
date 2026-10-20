#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
"${CXX:-c++}" -std=c++23 -O2 -DDP -include array -I"$root" \
   "$root/testpackage/alfven_face_averages.cpp" \
   "$root/backgroundfield/integratefunction.cpp" \
   "$root/backgroundfield/quadr.cpp" -o "$tmp/alfven_face_averages"
"$tmp/alfven_face_averages"
