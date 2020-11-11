#!/bin/bash

procs=4

rm -f *vlsv *txt

aprun -n $procs ../vlasiator \
  --run_config ./Flowthrough_amr.cfg \
  --io.restart_walltime_interval 10000 \
  --gridbuilder.timestep_max 5 \
  --loadBalance.rebalanceInterval 100 \
  --loadBalance.algorithm RANDOM

