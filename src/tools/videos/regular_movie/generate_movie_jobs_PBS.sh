#!/bin/bash

frameStart=0
frameEnd=1371
jobNumber=50
increment=$(( ($frameEnd - $frameStart) / $jobNumber ))

submit() {

qsub << EOF
#!/bin/bash
#PBS -N fieldlines
#PBS -l nodes=1
#PBS -l walltime=1:00:00
#PBS -W umask=007

cd $1

/lustre/tmp/yann/visit/bin/visit -lb-random -cli -nowin -l aprun -nn 1 -np 2 -s generate_frames.py $2 $3

EOF

}

dir=$( pwd )

for i in $( seq 0 $(( $jobNumber - 2 )) )
do
   start=$(( $frameStart + $i * $increment ))
   end=$(( $start + $increment ))
   submit $dir $start $end
done

submit $dir $end $frameEnd
