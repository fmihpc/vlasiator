#!/bin/bash
# Modify target JOBID below, then execute to run the inspector.sh on every node for every srun vlasiator instance.

if [ ! $# -eq 1 ]
then
    cat <<EOF
overlap_inspection.sh jobid
    Script for running inspector.sh on every node of a running job usign srun --overlap..
    
    jobid    SLURM job ID.
EOF
    exit
fi

export JOBID=$1

squeue --job $JOBID -o "%N" | grep nid | cut -b 5- | cut --delimiter "]" -f 1 > .nodelist

python3 - $( cat .nodelist ) > .fullnodelist << EOS
import sys

nodes_str_list = sys.argv[1].split(",")
outnodes = []

for n in nodes_str_list:
    if "-" in n:
        n_range = n.split("-")
        for i in range(int(n_range[0]),int(n_range[1])+1):
            outnodes.append("nid"+str(i).zfill(6))
    else:
        outnodes.append("nid"+n)

for t in outnodes:
    print(str(t))

EOS


for node in $( cat .fullnodelist )
do
  export NODE=$node
  srun --jobid $JOBID --overlap -w $node bash ./inspector.sh &
done

wait
