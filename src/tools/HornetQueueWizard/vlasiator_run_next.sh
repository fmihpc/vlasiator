#!/bin/bash

################
### FUNCTIONS
################

function get_next_job {
   next_job=$(cat $JOBLIST | head -n $job_to_check_in_list | tail -n 1)
}

function test_job {
   echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Testing whether to run $next_job..."
   cd $next_job
   
   #do not run if there is an ongoing simulation, or if the last has aborted for whatever reason 
   if [ -e logfile.txt ]
   then
      has_not_run_yet=0
      if [ -e I_AM_RUNNING ]
      then
         # Job running in another slot
         is_running=1
         echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $next_job is already running."
      else
         is_running=0
         # Has it exited cleanly?
         has_exited_cleanly=$( tail -2 logfile.txt |grep Exiting |wc -l )
         if [ $has_exited_cleanly -eq 1 ]
         then
            echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $next_job has exited cleanly before."
            # Has it reached the maximum time wished for in the cfg file?
            simulated=$( grep "dt =" logfile.txt | tail -n 1 | cut -d "=" -f 3 | cut -d " " -f 2 )
            wanted=$( grep "^t_max" Magnetosphere.cfg | cut -d "=" -f 2 )
            has_completed=$( echo $simulated $wanted | gawk '{if($1>$2) print 1; else print 0}' )
            if [ $has_completed -eq 1 ]
            then
               echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $next_job has completed according to cfg specifications."
            fi
         else
            echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $next_job has not exited cleanly before."
            has_completed=0
         fi
      fi
   else
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores)   $next_job has not run yet."
      is_running=0
      has_exited_cleanly=1
      has_completed=0
      has_not_run_yet=1
   fi
   
   if [[ $is_running -eq 0 && $has_exited_cleanly -eq 1 && $has_completed -eq 0 ]]
   then
      do_run=1
   else
      do_run=0
   fi
   
   cd $PBS_O_WRKDIR
}

function vlasiator_setup_next {
   cd $next_job
   
   # Create a unique cfg file
   cp Magnetosphere.cfg Magnetosphere.$PBS_JOBID.cfg
   
   # Set up restart interval and number of restarts
   # Subtract extra time for restart io
   RESTART_WALLTIME_INTERVAL=$( echo "scale=2; ($WALLTIME * 3600.0 - $RESTART_IO_EXTRA_TIME) / $NUMBER_OF_RESTARTS" | bc )
   
   sed -i'' -e 's/RESTART_WALLTIME_INTERVAL/'${RESTART_WALLTIME_INTERVAL}'/' Magnetosphere.$PBS_JOBID.cfg
   sed -i'' -e 's/NUMBER_OF_RESTARTS/'${NUMBER_OF_RESTARTS}'/' Magnetosphere.$PBS_JOBID.cfg
   
   # Set up the correct restart file
   last_restart=$( ls restart.*.vlsv | tail -1 )
   
   if [[ -z ${last_restart-unset} && $has_not_run_yet -eq 0 ]]
   then
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) No restart file! Not running."
      do_run=0
   else
      if [ $has_not_run_yet -eq 0 ]
      then
         echo " " >> Magnetosphere.$PBS_JOBID.cfg
         echo "[restart]" >> Magnetosphere.$PBS_JOBID.cfg
         echo "filename = $last_restart" >> Magnetosphere.$PBS_JOBID.cfg
      fi
   fi
   cd $PBS_O_WRKDIR
}

function vlasiator_run {
   cd $next_job

   if [[ $DRY_RUN -eq 1 ]]
   then
      echo "Dry run: I would be running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)."
   else
      touch I_AM_RUNNING
      # Launch the OpenMP job to the allocated compute node
      echo "Running $exec on $tasks mpi tasks, with $t threads per task on $nodes nodes ($ht threads per physical core)"
      aprun -n1 $exec --version
      aprun -n $NUM_PROCESSES -N $tasks_per_node -d $OMP_NUM_THREADS -j $ht ./vlasiator --run_config=Magnetosphere.$PBS_JOBID.cfg
      aprun -n $NUM_PROCESSES -N 1 /zhome/academic/HLRS/pri/ipryakem/vlasiator/tools/HornetQueueWizard/test.sh
      rm I_AM_RUNNING
   fi

   cd $PBS_O_WRKDIR
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
   get_next_job
   
   if [ "$job_to_check_in_list" -gt "$( cat $JOBLIST | wc -l )" -a $DRY_RUN -ne 1 ]
   then
      export message="($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Nothing can be run, exiting. Stop wasting queuing time!"
      echo $message
      echo $message | mailx -s "Job ended @Hornet" vlasiator-runs@fmihpc.flowdock.com
      exit
   fi
   
   if [ -z ${next_job-unset} ]
   then
      continue
   fi
   
   test_job

   if [ $do_run -eq 0 ]
   then
      echo "($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Tried $next_job, not doing this one (see above)."
   else
      export message="($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Running $next_job."
      echo $message
      
      vlasiator_setup_next
      if [ $do_run = 1 ]
      then
         if [ $DRY_RUN -ne 1 ]
         then
            echo $message | mailx -s "New job running @Hornet" vlasiator-runs@fmihpc.flowdock.com
         fi
         vlasiator_run
         doing_something=1
         export message="($(date) $(($NUM_PROCESSES*$OMP_NUM_THREADS)) cores) Done $next_job."
         echo $message
         if [ $DRY_RUN -ne 1 ]
         then
            echo $message | mailx -s "Job ended @Hornet" vlasiator-runs@fmihpc.flowdock.com
         fi
      else
         doing_something=0
      fi
   fi
done
