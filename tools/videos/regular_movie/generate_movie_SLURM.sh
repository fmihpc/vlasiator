#!/bin/bash


frameStart=0
frameEnd=2875
jobNumber=20
increment=$(( ($frameEnd - $frameStart) / $jobNumber ))

submit() {

sbatch << EOF
#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -J movie
#SBATCH -p serial
#SBATCH -n 16
#SBATCH -N 1
#SBATCH --no-requeue

~/visit_taito/bin/visit $1 $2 -lb-random -cli -nowin -l srun -nn 1 -np 16 -la "--mpi=pmi2"  -s generate_frames.py

EOF

}



for i in $( seq 0 $(( $jobNumber - 2 )) )
do
   start=$(( $frameStart + $i * $increment ))
   end=$(( $start + $increment ))
   submit $start $end &
done

submit $end $frameEnd &

