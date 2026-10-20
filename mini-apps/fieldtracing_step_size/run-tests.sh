#!/usr/bin/env bash
set -eu
ulimit -c 0
directory=$(cd "$(dirname "$0")" && pwd)
temporary=$(mktemp -d)
trap 'rm -rf "$temporary"' EXIT
"${CXX:-c++}" -std="${STANDARD:-c++20}" -Wall -Wextra -Werror -pedantic "$directory/step_size_test.cpp" -o "$temporary/test"
"$temporary/test"
for precision in float double; do
   for mode in reversed zero-min negative-min zero-max zero-step negative-step nan-min nan-max nan-step inf-min inf-max inf-step; do
      status=0
      { "$temporary/test" "$precision" "$mode" > "$temporary/output" 2>&1; } 2>/dev/null || status=$?
      if [ "$status" -ne 134 ]; then
         printf 'FAIL: %s %s exited %s instead of SIGABRT\n' "$precision" "$mode" "$status" >&2
         exit 1
      fi
      grep -Fq '(fieldtracing) Error: step lengths must be finite and positive' "$temporary/output"
   done
done
printf 'PASS: 24 invalid step-size cases rejected (float and double)\n'
