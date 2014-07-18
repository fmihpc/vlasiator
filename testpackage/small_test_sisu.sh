#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -J vlasiator-test
#SBATCH -o test.%j
#SBATCH -e test.%j
#SBATCH -p test
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH --no-requeue
# here we just ask SLURM for N full nodes, and tell below (to ALPS/aprun) how to use them

t=4
run_command="aprun -n 16 -d 4 -N 4"


#folder for reference data 
reference_dir="/homeappl/home/alfthan/vlasiator-reference/sandybridge-master-c46ee5f53a46d020f1066c8f8708818518a37440"

#If 1, the reference vlsv files are generated
# if 0 then we check the validity against the reference
create_verification_files=0

# Define test
source small_test_definitions.sh
wait
# Run tests
source run_tests.sh
