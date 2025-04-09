# script run by srun overlap when overlap_inspection.sh is run.
# As it can be slow at scale, all are run concurrently and we wait at the end.

for proc in $( ps x | grep vlasiator | grep -v srun | grep -v grep | cut --delimiter=" " -f 1 )
do
  # add more -ex commands, or pass all commands in myfile as -x myfile
  gdb -p $proc --batch -ex bt -ex q &> ${JOBID}_${NODE}_${proc} &
done

wait
