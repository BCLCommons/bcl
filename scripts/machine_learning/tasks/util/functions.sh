#!/bin/bash
# Commands used by BclModelCrossValidation
# author Jeffrey Mendenhall
# date April 4, 2013

# Global Variables used by this script for returning values from functions:
# LOGGER_FILE - Updated by GetBclLoggerFile, indicates log file for BCL executable
# JOB_COUNT - Updated by UpdateJobCount, indicates number of pids (given as args to the function) still running

# Global variables used by this script as inputs:
# BCL_MAX_RERUNS - maximum number of types to retry each bcl execution. Typically 1 is appropriate, but on really finicky
#                  file systems 2 might be better
# BCL_CV_STATUS_FILENAME - 

# choose a random wait time up to the max wait time to actually wait for this lock. If the times were all exactly equal, 
# then it would still be relatively likely for processes to end up with the same wait time during the prolog
FLOCK_WAIT_TIME=$(( ( `od -A n -t dI -N 4 /dev/urandom | tr -d -` % $FLOCK_MAX_WAIT) + 1))
FLOCK_FAIL_MESSAGE="flock failed to lock $BCL_CV_STATUS_LOCK_FILENAME after $FLOCK_WAIT_TIME seconds on $HOSTNAME. Rarely, this could result in status file corruption. Request cluster administrators to restart lockd and statd on $HOSTNAME if this happens often"

# this function returns the name of the file logged to by the bcl, assuming the command is of the form ... --logger File 
LOGGER_FILE=""
GetBclLoggerFile()
{
  while [ "$1" != "-logger" ] && [ "$#" != "0" ]; do
    shift 1
  done
  if [ "$1" == "-logger" ] && [ "$2" == "File" ]; then
    LOGGER_FILE="$3"
  fi
}

# Count the number of PIDs (given as arguments) that are still running, maintain number in JOB_COUNT variable
JOB_COUNT=0
UpdateJobCount()
{
  cnt=0
  for var in "$@"
  do
    kill -0 $var 2> /dev/null
    es=$?
    if [[ "$es" -eq "0" ]]
    then
      cnt=$(( cnt + 1 ))
    fi
  done
  JOB_COUNT=$cnt
}

# this function updates the status file to indicate that another job is running
# argument (1) 
PrologueStatus()
{
  log_file=$BCL_CV_STATUS_FILENAME
  (
    flock -x -w $FLOCK_WAIT_TIME 200
    [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
    mv $log_file $log_file.back
    awk '{print ($1-1)" "($2+1)" "$3" "$4" "$5" "$6;}' $log_file.back > $log_file
    rm $log_file.back
  ) 200> $BCL_CV_STATUS_LOCK_FILENAME
}

# Increment a file containing a single number using a file lock to handle race conditions
IncrementBundleCount()
{
  filnm=$BCL_CV_SCRIPTS_PATH/bundlecount.$1.txt
  (
    flock -x -w $FLOCK_WAIT_TIME 201
    [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
    nmbr=`cat $filnm`
    nmbrp1=$(( nmbr + 1 ))
    echo $nmbrp1 > $filnm
    rm -f $filnm.wakeup
  ) 201> $filnm.lock
}

# Decrement a file containing a single number using a file lock to handle race conditions
DecrementBundleCount()
{
  filnm=$BCL_CV_SCRIPTS_PATH/bundlecount.$1.txt
  (
    flock -x -w $FLOCK_WAIT_TIME 201
    [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
    nmbr=`cat $filnm`
    nmbrm1=$(( nmbr - 1 ))
    echo $nmbrm1 > $filnm
    if [[ "$nmbrm1" -le 0 ]]; then
      touch $filnm.wakeup
    fi
  ) 201> $filnm.lock
}

# this function is used to check that all merging happened properly
# arguments: 1 - number of reruns left, 2 - status filename, 3 - actual command to run
run_merge_command()
{
  n_reruns_left=$BCL_MAX_RERUNS
  "$@"
  error_status=$?
  if [[ $error_status -ne 0 ]]
    then
    while [[ $error_status -ne 0 ]] && [[ $n_reruns_left -ne 0 ]]
    do
      echo "Command $@ failed with status $?, retrying"
      sleep 3
      "$@"
      error_status=$?
      n_reruns_left=$(($n_reruns_left - 1))
    done
    if [[ $error_status -ne 0 ]]
    then
      echo "aborting with exit status $error_status after running $*"
      GetBclLoggerFile $*
      echo "File: $LOGGER_FILE" >> $BCL_CV_STATUS_FILENAME.bclerr
      tail $LOGGER_FILE >> $BCL_CV_STATUS_FILENAME.bclerr
      exit $error_status
    fi
  fi
  return 0
}

# this function performs all epilogue script duties given the job filename
PbsEpilogue()
{
  if [[ ${10} -ne 0 ]]
  then
    echo "Job Failed! $@" >> $1
    (
      flock -x -w $FLOCK_WAIT_TIME 200
      [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
      mv $BCL_CV_STATUS_FILENAME $BCL_CV_STATUS_FILENAME.back
      awk '{print $1" "($2-1)" "($3+1)" "$4" "($5+1)" "$6;}' $BCL_CV_STATUS_FILENAME.back > $BCL_CV_STATUS_FILENAME
      rm $BCL_CV_STATUS_FILENAME.back
      pbs_err=`grep -i 'PBS: job killed' $1 | sed 's/.*job killed//'`
      if [ "$pbs_err" != "" ]; then
        echo "$5 $2 killed: $pbs_err" >> $BCL_CV_STATUS_FILENAME.pbserr
      else
        echo -ne "Job Name: $5\nResource List: $7\nResources Used: $8\nJob exit code: ${10}"  >> $BCL_CV_STATUS_FILENAME.pbserr
      fi
      touch $BCL_CV_WAKEUP_FILENAME
    ) 200> $BCL_CV_STATUS_LOCK_FILENAME
  else
    echo "Job ID: $2 succeeded" >> $1
    echo "Resource List: $7"  >> $1
    echo "Resources Used: $8"  >> $1
  fi
}

# this function performs all epilogue script duties given the job filename
# It takes a single argument: the file given to --output in the srun command
SLURMEpilogue()
{
  # wait around for various unsundry parts of SLURM to finally notice that the job is done and update their info on how
  # many resources that it consumed. 5 seconds isn't always long enough; fortunately this info isn't always critical 
  # and can be retrieved by the user later if they need to. Making it longer than 5 seconds runs a real risk that
  # SLURM decides to kill the epilogue before it's done (which happens after 10 seconds in the default configuration)
  # This is how SLURM rolls: lazy updates on info that counts,
  # but kill processes that might wait around for that info to be available
  sleep 5
  INFO=`scontrol -o show job $SLURM_JOB_ID`
  JOBSTATE=`echo "$INFO" | tr ' ' '\n' | grep '^JobState' | sed 's/.*=//'`
  EXITCODE_JOB_SIGNAL=`echo "$INFO" | tr ' ' '\n' | grep '^ExitCode' | sed 's/.*=//' | tr ':' ' '`
  EXITCODE_JOB=`echo $EXITCODE_JOB_SIGNAL | awk '{print $1}'`
  EXITCODE_SLURM=`echo $EXITCODE_JOB_SIGNAL | awk '{print $2}'`
  EXITCODE_MAX=`echo $EXITCODE_JOB_SIGNAL | awk '{print ($1 != 0 ? $1 : $2)}'`
  USEFULINFO=`echo "$INFO" | tr ' ' '\n' | egrep '^(JobId|JobName|JobState|RunTime|TimeLimit|StartTime|EndTime|NodeList|BatchHost|NumNodes|NumCPUs|NumTasks|CPUs/Task|MinMemoryCPU|MinMemoryNode|Command)='`
  # call sacct to retrieve completed job information. Reformat it using an obscenely long and obtuse awk command since
  # the default output format is presently practically unreadable
  awk_command='{if(NR==1){for( x=1; x <= NF; ++x){FIELD_NAMES[x-1]=$x;}}else if(NR == 2){for( x=1; x <= NF; ++x){FIELD_SIZES[x-1]=length($x);}}else{ y=1; for( x=1; x <= length(FIELD_SIZES); ++x){print FIELD_NAMES[x-1],substr($0,y,FIELD_SIZES[x-1]); y += FIELD_SIZES[x-1]+1;}}}'
  full_out=`sacct -j $SLURM_JOB_ID -o JobID,User,JobName,ExitCode,DerivedExitCode,AveCPU,MaxRSS,Timelimit,UserCPU,SystemCPU,MaxPages | awk "$awk_command"`
  if [[ ${EXITCODE_MAX} -ne 0 || "$JOBSTATE" -ne "COMPLETED" ]]; then
    echo "Job Failed! $INFO" >> $1
    (
      flock -x -w $FLOCK_WAIT_TIME 200
      [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
      grep "slurmstepd: error: " $1 >> $BCL_CV_STATUS_FILENAME.pbserr
      mv $BCL_CV_STATUS_FILENAME $BCL_CV_STATUS_FILENAME.back
      awk '{print $1" "($2-1)" "($3+1)" "$4" "($5+1)" "$6;}' $BCL_CV_STATUS_FILENAME.back > $BCL_CV_STATUS_FILENAME
      rm $BCL_CV_STATUS_FILENAME.back      
      if [[ $EXITCODE_SLURM -ne 0 || $EXITCODE_JOB -eq 9 || $EXITCODE_JOB -eq 137 ]]; then
        echo -ne "$SLURM_JOB_NAME $SLURM_JOB_ID killed: $JOBSTATE $USEFULINFO\n" >> $BCL_CV_STATUS_FILENAME.pbserr
      else
        echo -ne "Job died (bad flag? check log file) $USEFULINFO" >> $BCL_CV_STATUS_FILENAME.pbserr
      fi
      echo -ne "$full_out" >> $1
      echo -ne "$full_out" >> $BCL_CV_STATUS_FILENAME.pbserr
      touch $BCL_CV_WAKEUP_FILENAME
    ) 200> $BCL_CV_STATUS_LOCK_FILENAME
  else
    echo "Job ID: $SLURM_JOB_ID $JOBSTATE" >> $1
    echo "$USEFULINFO"  >> $1
    echo "$full_out"  >> $1
    echo "ECM: $EXITCODE_MAX"  >> $1
  fi
}

# this function updates the status file to indicate that 1 less job is running and update the # of errors
EpilogueStatus()
{
  error_status=$1
  log_file=$BCL_CV_STATUS_FILENAME
  (
    flock -x -w $FLOCK_WAIT_TIME 200
    [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
    mv $log_file $log_file.back
    if [[ $error_status -ne 0 ]]; then
      awk '{print $1" "($2-1)" "($3+1)" "($4+1)" "$5" "$6;}' $log_file.back > $log_file
    else
      awk '{print $1" "($2-1)" "($3+1)" "$4" "$5" "$6;}' $log_file.back > $log_file
      tail -n30 $LOGGER_FILE | grep 'independent data: ' | sed 's/.*: \([-0-9.eE+]*\) time:.*/ \1/' >> $BCL_CV_STATUS_FILENAME.ind_values.txt
      python -c 'import sys; print sum([ float(x) for x in sys.argv[1:]]) / float(len(sys.argv) - 1)' `cat $BCL_CV_STATUS_FILENAME.ind_values.txt | tr '\n' ' '` > $BCL_CV_STATUS_FILENAME.raw_ind_ave.txt
    fi
    rm $log_file.back
    n_noncomplete=`awk '{print $2+$1;}' $log_file`
    n_error=`awk '{print $4;}' $log_file`
    if [[ $n_error -ne 0 ]]; then
      touch $BCL_CV_WAKEUP_FILENAME
    elif [[ $n_noncomplete -eq 0 ]]; then
      if [[ -e "$BCL_MERGE_SCRIPT" ]]; then
        # run the merge script
        mv $log_file $log_file.back
        awk '{print $1" "$2" "$3" "$4" "$5" 2";}' $log_file.back > $log_file
        rm $log_file.back
        $BCL_MERGE_SCRIPT &> $BCL_CV_LOGFILES_PATH/log_merge.job.txt
        merge_error_status=$?
        mv $log_file $log_file.back   
        if [[ $merge_error_status -ne 0 ]]; then
          awk '{print $1" "$2" "$3" "($4+1)" "$5" -1";}' $log_file.back > $log_file
        else
          awk '{print $1" "$2" "$3" "$4" "$5" 1";}' $log_file.back > $log_file
        fi
        rm $log_file.back
      fi 
      touch $BCL_CV_WAKEUP_FILENAME
    fi
  ) 200> $BCL_CV_STATUS_LOCK_FILENAME
}

# this function will be called if there was an error running the job on slurm
EpilogueStatusSlurm()
{
  error_status=$1
  log_file=$BCL_CV_STATUS_FILENAME
  (
    flock -x -w $FLOCK_WAIT_TIME 200
    [[ $? -eq 0 ]] || echo $FLOCK_FAIL_MESSAGE
    mv $log_file $log_file.back
    if [[ $error_status -eq 9 || $error_status -eq 137 ]]; then
      awk '{print $1" "($2-1)" "($3+1)" "$4" "($5+1)" "$6;}' $log_file.back > $log_file
    elif [[ $error_status -ne 0 ]]; then
      awk '{print $1" "($2-1)" "($3+1)" "($4+1)" "$5" "$6;}' $log_file.back > $log_file
    else
      awk '{print $1" "($2-1)" "($3+1)" "$4" "$5" "$6;}' $log_file.back > $log_file
      tail $LOGGER_FILE | grep 'independent data: ' | sed 's/.*: \([-0-9.eE+]*\) time:.*/ \1/' >> $BCL_CV_STATUS_FILENAME.ind_values.txt
      python -c 'import sys; print sum([ float(x) for x in sys.argv[1:]]) / float(len(sys.argv) - 1)' `cat $BCL_CV_STATUS_FILENAME.ind_values.txt | tr '\n' ' '` > $BCL_CV_STATUS_FILENAME.raw_ind_ave.txt
    fi
    rm $log_file.back
    touch $BCL_CV_WAKEUP_FILENAME
  ) 200> $BCL_CV_STATUS_LOCK_FILENAME
}

# this function runs a command (given as arguments to the function)
# if the command fails, after up to the max # of rerun attempt times, exit with the same status as the command
run_command_check_status()
{
  "$@"
  error_status=$?
  if [[ $error_status -ne 0 ]]
    then
    n_reruns_left=$BCL_MAX_RERUNS
    while [[ $error_status -ne 0 ]] && [[ $n_reruns_left -ne 0 ]]
    do
      echo "Command $@ failed with status $?, retrying"
      sleep 3
      "$@"
      error_status=$?
      n_reruns_left=$(($n_reruns_left - 1))
    done
    if [[ $error_status -ne 0 ]]
    then
      echo "aborting with exit status $error_status after running $*"
      GetBclLoggerFile $*
      echo "File: $LOGGER_FILE" >> $BCL_CV_STATUS_FILENAME.bclerr
      tail $LOGGER_FILE >> $BCL_CV_STATUS_FILENAME.bclerr
      EpilogueStatus $error_status
      exit 0
    fi
  fi
  GetBclLoggerFile $*
  EpilogueStatus $error_status 
  return 0
}

# this function runs a command (given as arguments to the function)
# if the command fails, after up to the max # of rerun attempt times, exit with the same status as the command
run_command_check_status_slurm()
{
  "$@"
  error_status=$?
  if [[ $error_status -ne 0 ]]
    then
    n_reruns_left=$BCL_MAX_RERUNS
    while [[ $error_status -ne 0 ]] && [[ $error_status -ne 9 ]] && [[ $error_status -ne 137 ]] &&[[ $n_reruns_left -ne 0 ]]
    do
      echo "Command $@ failed with status $?, retrying"
      sleep 3
      "$@"
      error_status=$?
      n_reruns_left=$(($n_reruns_left - 1))
    done
    if [[ $error_status -ne 0 ]]
    then
      echo "aborting with exit status $error_status after running $*"
      GetBclLoggerFile $*
      if [[ $error_status -ne 9 ]] && [[ $error_status -ne 137 ]]
      then
        echo "File: $LOGGER_FILE" >> $BCL_CV_STATUS_FILENAME.bclerr
        tail $LOGGER_FILE >> $BCL_CV_STATUS_FILENAME.bclerr
        EpilogueStatus $error_status
      else
        JOB_FILE=`echo $LOGGER_FILE | sed 's/\.txt$/.job.txt/'`
        EpilogueStatusSlurm $error_status
      fi
      # set the exit status for slurm's purposes
      exit $error_status
    fi
  fi
  GetBclLoggerFile $*
  EpilogueStatus $error_status 
  return 0
}

# runs a bundle of jobs on PBS nodes
RunPBSBundle()
{
  cd $PBS_O_WORKDIR
  echo "Running jobs on "`cat $PBS_NODEFILE`
  BUNDLESTART=$(( BUNDLE * $1 ))
  echo 0 > $BCL_CV_SCRIPTS_PATH/bundlecount.$BUNDLE.txt
  pids=()
  for JB in `seq 0 $(( $1 - 1 ))`
  do
    SCRIPT=`sed -n $(( JB + BUNDLESTART + 1 ))p $BCL_CV_SCRIPTS_PATH/bundles.txt`
    if [ -n "$SCRIPT" ]; then
      pbsdsh -n $JB $SCRIPT &
      pids+=($!)
      IncrementBundleCount $BUNDLE
      JOB_COUNT=$(( JOB_COUNT + 1 ))
    fi
  done
  while [ ! -f $BCL_CV_SCRIPTS_PATH/bundlecount.$BUNDLE.txt.wakeup ]
  do 
    UpdateJobCount "${pids[@]}"
    if [ "$JOB_COUNT" -eq "0" ]; then
      exit 0
    fi
    sleep 60
  done
  kill -9 "${pids[@]}" 2> /dev/null
}

# runs a bundle of jobs on PBS nodes
RunSLURMBundle()
{
  cd $SLURM_SUBMIT_DIR
  echo "Running jobs on "`cat $SLURM_NODELIST`
  BUNDLESTART=$(( BUNDLE * $1 ))
  echo 0 > $BCL_CV_SCRIPTS_PATH/bundlecount.$BUNDLE.txt
  pids=()
  for JB in `seq 0 $(( $1 - 1 ))`
  do
    SCRIPT=`sed -n $(( JB + BUNDLESTART + 1 ))p $BCL_CV_SCRIPTS_PATH/bundles.txt`
    if [ -n "$SCRIPT" ]; then
      srun $SCRIPT &
      pids+=($!)
      IncrementBundleCount $BUNDLE
      JOB_COUNT=$(( JOB_COUNT + 1 ))
    fi
  done
  while [ ! -f $BCL_CV_SCRIPTS_PATH/bundlecount.$BUNDLE.txt.wakeup ]
  do 
    UpdateJobCount "${pids[@]}"
    if [ "$JOB_COUNT" -eq "0" ]; then
      exit 0
    fi
    sleep 60
  done
  kill -9 "${pids[@]}" 2> /dev/null
}

# remove redundant descriptor files. This is not to be used if the descriptors are not all the same; the assumption is 
# that they are
remove_duplicate_descriptor_files()
{
  cd $1
  descriptor_files=`ls model[0-9][0-9][0-9][0-9][0-9][0-9].descriptor* 2> /dev/null`
  nfiles=`echo -n "$descriptor_files" | wc -l`
  if [[ "$nfiles" -gt "0" ]]
  then
    echo "Removing duplicate descriptor files in "$1
    first_file=`echo -n "$descriptor_files" | head -n1`
    compression=`echo -n "$first_file" | sed 's/^model[0-9][0-9][0-9][0-9][0-9][0-9].descriptor//'`
    if [[ "$compression" == "" ]]
    then
      mv $first_file model.descriptor
      bzip2 -f model.descriptor
    else
      mv $first_file model.descriptor$compression
    fi
    # remove the extra files
    echo "$descriptor_files" | tail -n+2 | xargs -r -n1 rm -f
  fi
  cd - >& /dev/null
}

# Merge duplicate descriptor files; compress all model files except the miniscule result file
# This function is useful if you have model directories generated by an old version of the BCL, where the models are 
# not compressed
condense_models_directory()
{
  remove_duplicate_descriptor_files $1
  cd $1
  uncompressed_stuff=`ls model[0-9][0-9][0-9][0-9][0-9][0-9].model model[0-9][0-9][0-9][0-9][0-9][0-9].info 2> /dev/null`
  if [[ `echo -n "$uncompressed_stuff" | wc -l` -ne "0" ]]
  then
    echo "Compressing model and info files in "$1
    echo -n "$uncompressed_stuff" | xargs -r -n1 bzip2 -f
  fi
  cd - >& /dev/null
}
