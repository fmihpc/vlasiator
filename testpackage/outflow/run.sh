#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
for ranks in 1 2; do
   for scenario in single different reversed identical empty schemes; do
      cfg=different
      case "$scenario" in
         single|empty) cfg=empty ;;
         identical|schemes) cfg=$scenario ;;
      esac
      mpirun -np "$ranks" "$root/outflow-test-bin" "$scenario" \
         --run_config="$root/testpackage/outflow/$cfg.ini"
   done
done
if [ -n "${OUTFLOW_TEST_LOG_DIR:-}" ]; then
   scratch=$OUTFLOW_TEST_LOG_DIR
   mkdir -p "$scratch"
else
   scratch=$(mktemp -d)
   trap 'rm -rf "$scratch"' EXIT
fi
log="$scratch/invalid-scheme.log"
if mpirun -np 1 "$root/outflow-test-bin" schemes \
   --run_config="$root/testpackage/outflow/invalid-scheme.ini" > "$log" 2>&1; then
   echo "Invalid outflow scheme was accepted" >&2
   exit 1
fi
grep -q 'invalid Outflow Vlasov scheme' "$log"
echo "OUTFLOW_INVALID_SCHEME rejected"
for ranks in 1 2; do
   mkdir "$scratch/valid-$ranks" "$scratch/invalid-$ranks"
   (
      cd "$scratch/valid-$ranks"
      mpirun -np "$ranks" "$root/vlasiator" --run_config="$root/testpackage/outflow/solver-valid.ini" > stdout 2>&1
      if grep -Ei '\(BAILOUT\):|bailing out' stdout logfile.txt; then exit 1; fi
      grep -Eq 'tstep = 2([[:space:]]|$)' logfile.txt
   )
   (
      cd "$scratch/invalid-$ranks"
      if mpirun -np "$ranks" "$root/vlasiator" --run_config="$root/testpackage/outflow/solver-missing-species.ini" > stdout 2>&1; then
         echo "Incomplete species coverage was accepted" >&2
         exit 1
      fi
      grep -Fq "Missing system boundary condition for population 'proton' on face z+" stdout
   )
   echo "OUTFLOW_COVERAGE ranks=$ranks valid=step2 incomplete=rejected"
done
bash "$root/testpackage/outflow/run-mixed.sh"
