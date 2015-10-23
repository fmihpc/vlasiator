#!/bin/bash

#SBATCH -t 00:30:00
#SBATCH -J movie
#SBATCH -p parallel
#SBATCH -n 64
#SBATCH -N 4
#SBATCH --no-requeue

frameStart=2150
frameEnd=2200
jobNumber=4
increment=$(( ($frameEnd - $frameStart) / $jobNumber ))

for i in $( seq 0 $(( $jobNumber - 2 )) )
do
   start=$(( $frameStart + $i * $increment ))
   end=$(( $start + $increment ))
   ~/visit_taito/bin/visit -lb-random -cli -nowin -l srun -nn 1 -np 16 -s generate_frames.py $start $end &
   sleep 5
done

~/visit_taito/bin/visit -lb-random -cli -nowin -l srun -nn 1 -np 16 -s generate_frames.py $end $frameEnd

wait
