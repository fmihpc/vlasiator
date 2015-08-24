#!/bin/bash

frameStart=0
frameEnd=4000
jobNumber=10
increment=$(( ($frameEnd - $frameStart) / $jobNumber ))

for i in $( seq 0 $(( $jobNumber - 2 )) )
do
   start=$(( $frameStart + $i * $increment ))
   end=$(( $start + $increment ))
   /lustre/tmp/yann/visit/bin/visit -lb-random -cli -nowin -l aprun -nn 1 -np 20 -s generate_frames.py $start $end &
   sleep 5
done

/lustre/tmp/yann/visit/bin/visit -lb-random -cli -nowin -l aprun -nn 1 -np 20 -s generate_frames.py $end $frameEnd
