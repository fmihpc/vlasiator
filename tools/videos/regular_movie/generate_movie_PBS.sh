#!/bin/bash

#PBS -N movie
#PBS -l nodes=5
#PBS -l walltime=0:25:00
#PBS -W umask=007

frameStart=4238
frameEnd=4300
jobNumber=5
increment=$(( ($frameEnd - $frameStart) / $jobNumber ))

cd $PBS_O_WORKDIR

for i in $( seq 0 $(( $jobNumber - 2 )) )
do
   start=$(( $frameStart + $i * $increment ))
   end=$(( $start + $increment ))
   /zhome/academic/HLRS/pri/ipryakem/visit/bin/visit $start $end -lb-random -cli -nowin -l aprun -nn 1 -np 24 -s generate_frames.py &
   sleep 3
done

/zhome/academic/HLRS/pri/ipryakem/visit/bin/visit $end $frameEnd -lb-random -cli -nowin -l aprun -nn 1 -np 24 -s generate_frames.py

wait

