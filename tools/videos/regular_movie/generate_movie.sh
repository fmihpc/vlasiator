#!/bin/bash

frameStart=0
frameEnd=4000
jobNumber=10
increment=$(( ($frameEnd - $frameStart) / $jobNumber ))

for i in $( seq 0 $(( $jobNumber - 2 )) )
do
   start=$(( $frameStart + $i * $increment ))
   end=$(( $start + $increment ))
   /lustre/tmp/yann/visit/bin/visit -lb-random -cli -nowin -debug 1 -s generate_frames.py $start $end &
   sleep 15
done

/lustre/tmp/yann/visit/bin/visit -lb-random -cli -nowin -debug 1 -s generate_frames.py $end $frameEnd
