#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
if [ -n "${OUTFLOW_TEST_LOG_DIR:-}" ]; then
   scratch=$OUTFLOW_TEST_LOG_DIR/mixed
   mkdir -p "$scratch"
else
   scratch=$(mktemp -d)
   trap 'rm -rf "$scratch"' EXIT
fi
common=("--run_config=$root/testpackage/outflow/solver-mixed.ini")
for pop in proton helium; do
   for face in x+ x- z+ z-; do
      common+=("--${pop}_maxwellian.file_${face}=$root/testpackage/outflow/$pop-input.txt")
   done
done
for ranks in 1 2; do
   for scenario in outflow-higher maxwellian-higher tie-maxwellian-last tie-outflow-last explicit-none corner-outflow-complete corner-maxwellian-higher corner-tie-outflow-last; do
      extra=()
      reject=no
      case "$scenario" in
         outflow-higher) reject=yes ;;
         maxwellian-higher) extra+=(--outflow.precedence=2) ;;
         tie-maxwellian-last) extra+=(--outflow.precedence=3) ;;
         tie-outflow-last|corner-tie-outflow-last)
            extra+=(--outflow.precedence=3 --boundaries.boundary=Maxwellian --boundaries.boundary=Outflow)
            reject=yes ;;
         explicit-none)
            extra+=(--proton_outflow.face=z+ --proton_outflow.vlasovScheme_face_z+=None) ;;
         corner-outflow-complete) extra+=(--proton_outflow.face=z+) ;;
         corner-maxwellian-higher) extra+=(--outflow.precedence=2) ;;
      esac
      if [[ "$scenario" = corner-* ]]; then
         extra+=(--boundaries.periodic_x=no --maxwellian.face=x+ --maxwellian.face=x- --maxwellian.face=z+ --maxwellian.face=z-)
      fi
      mkdir "$scratch/$scenario-$ranks"
      (
         cd "$scratch/$scenario-$ranks"
         if [ "$reject" = yes ]; then
            if mpirun -np "$ranks" "$root/vlasiator" "${common[@]}" "${extra[@]}" > stdout 2>&1; then
               echo "Mixed case $scenario accepted incomplete coverage" >&2
               exit 1
            fi
            grep -Fq "Missing system boundary condition for population 'proton' on face z+ (effective boundary: Outflow)" stdout
         else
            mpirun -np "$ranks" "$root/vlasiator" "${common[@]}" "${extra[@]}" > stdout 2>&1
            if grep -Ei '\(BAILOUT\):|bailing out' stdout logfile.txt; then exit 1; fi
            grep -Eq 'tstep = 2([[:space:]]|$)' logfile.txt
         fi
      )
      echo "OUTFLOW_MIXED ranks=$ranks scenario=$scenario reject=$reject pass"
   done
done
