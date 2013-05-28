#!/bin/bash

################
### FUNCTIONS
################

function get_next_job {
   next_job=$(cat $JOBLIST | head -n $1 | tail -n 1)
}

function test_job {
   job=$1
   echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Testing whether to run $job..."
   cd $job
   
   #do not run if there is an ongoing simulation, or if the last has aborted for whatever reason 
   if [ -e logfile.txt ]
   then
      if [ -e *.OU ]
      then
         # Job running in another slot
         is_running=1
         echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $job is already running."
      else
         is_running=0
         # Has it exited cleanly?
         has_exited_cleanly=$( tail -2 logfile.txt |grep Exiting |wc -l )
         if [ $has_exited_cleanly -eq 1 ]
         then
            echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $job has exited cleanly before."
            # Has it reached the maximum time wished for in the cfg file?
            simulated=$( grep "total simulated" logfile.txt | cut -d " " -f 19 )
            wanted=$( grep "t_max" Magnetosphere.cfg | cut -d " " -f 3 )
            has_completed=$( echo $simulated $wanted | gawk '{if($1>$2) print 1; else print 0}' )
            if [ $has_completed -eq 1 ]
            then
               echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $job has completed according to cfg specifications."
            fi
         else
            echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $job has not exited cleanly before."
            has_completed=0
         fi
      fi
   else
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $job has not run yet."
      is_running=0
      has_exited_cleanly=1
      has_completed=0
   fi
   
   if [[ $is_running -eq 0 && $has_exited_cleanly -eq 1 && $has_completed -eq 0 ]]
   then
      do_run=1
   else
      do_run=0
   fi
   
   cd ..
}

function vlasiator_setup_next {
   next_job=$1
   cd $next_job
   
   # Create a unique cfg file
   cp Magnetosphere.cfg Magnetosphere.$BATCH_JOBID.cfg
   
   # Set up restart interval and number of restarts
   # Subtract extra time for restart io
   RESTART_WALLTIME_INTERVAL=$( echo "scale=2; ($WALLTIME * 3600.0 - $RESTART_IO_EXTRA_TIME) / $NUMBER_OF_RESTARTS" | bc )
   
   sed -i'' -e 's/RESTART_WALLTIME_INTERVAL/'${RESTART_WALLTIME_INTERVAL}'/' Magnetosphere.$BATCH_JOBID.cfg
   sed -i'' -e 's/NUMBER_OF_RESTARTS/'${NUMBER_OF_RESTARTS}'/' Magnetosphere.$BATCH_JOBID.cfg
   
   # Set up the correct restart file
   last_restart=$( ls restart.*.vlsv | tail -1 )
   
   if [ -z ${last_restart-unset} ]
   then
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) No restart file! Not running."
      do_run=0
   else
      echo " " >> Magnetosphere.$BATCH_JOBID.cfg
      echo "[restart]" >> Magnetosphere.$BATCH_JOBID.cfg
      echo "filename = $last_restart" >> Magnetosphere.$BATCH_JOBID.cfg
   fi

}

function vlasiator_run {
   # Launch the OpenMP job to the allocated compute node
   aprun -n $NUM_PROCESSES -N $((32/$OMP_NUM_THREADS)) -d $OMP_NUM_THREADS ./vlasiator --run_config=Magnetosphere.$BATCH_JOBID.cfg
   cd ..
}

####################
### HERE WE START
####################

cd $PBS_O_WORKDIR

doing_something=0
job_to_check_in_list=0
while [ $doing_something -eq 0 ]
do
   # Getting next job  to do from $JOBLIST
   let job_to_check_in_list=$job_to_check_in_list+1
   get_next_job $job_to_check_in_list
   
   if [ -z ${next_job-unset} ]
   then
      continue
   fi
   
   if [ "$job_to_check_in_list" -gt "$( cat $JOBLIST | wc -l )" ]
   then
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Nothing can be run, exiting. Stop wasting queuing time!"
      exit
   fi

   test_job $next_job

   if [ $do_run -eq 0 ]
   then
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Tried $next_job, another job is running, or it has exited uncleanly, or it has completed already. Not doing this one."
   else
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Running $next_job."
      
      vlasiator_setup_next $next_job
      if [ $do_run = 1 ]
      then
         vlasiator_run $next_job
         doing_something=1
         echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Done $next_job."
      else
         doing_something=0
      fi
   fi
done
