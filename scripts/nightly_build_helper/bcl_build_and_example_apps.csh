#!/bin/tcsh

# make sure that at least one argument was given and that not too many arguments were given
if( $#argv > 4 ) then
  echo "Usage: $0 [architecture] [results_folder] [build_machine] [run_machine]"
  echo "Script should be run from the checkout folder"
  echo "  architecture: uses by default cmake toolchain x86_64-unknown-linux-gnu (only the architecture is needed)"
  echo "  results_folder: uses by default ./example_results to output the results of each example to"
  echo "  build_machine and run_machine: uses by default the local machine that the script is run from"
  echo "This script will attempt to run all the examples in the current environment; if that fails, it will try running them under wine"
  echo "Your command line was: $*"
  exit 1
endif

# set the arguments to named variables
set TOOLCHAIN = $1
set RESULTS_DIR = $2
set BUILD_MACHINE = $3
set RUN_MACHINE = $4
set CHECKOUT_DIR = `pwd`
set BUILD_TARGET = "bcl-example-apps-static"

# set default toolchain if none was given
if( $#argv == 0 ) then
  set TOOLCHAIN = "x86_64-unknown-linux-gnu"
endif
	
set EXECUTABLE = "./build/${TOOLCHAIN}/bin/${BUILD_TARGET}.exe"

# set default results folder if non was given
if( $#argv <= 1 ) then
  set RESULTS_DIR = "${CHECKOUT_DIR}/example_results"
endif
# create results folder 
if( ! -d ${RESULTS_DIR} ) then
  echo "STATUS: Creating results directory ${RESULTS_DIR}"
  mkdir -p ${RESULTS_DIR}
else
  echo "STATUS: Using results directory ${RESULTS_DIR}"
  rm -f ${RESULTS_DIR}/output_App*
  rm -f ${RESULTS_DIR}/example_apps_make_log.txt
endif
	
# get absolute folder path
cd ${RESULTS_DIR}
set RESULTS_DIR = `pwd`
cd -
# set log file
set REDIRECT_FILE = "${RESULTS_DIR}/example_apps_make_log.txt"

# set cmake args
set CMAKE_ARGS = '-DCMAKE_BUILD_TYPE=Release'

# build bcl
set BUILD_COMMAND = "./scripts/nightly_build_helper/bcl_build.csh ${BUILD_TARGET} ${TOOLCHAIN} ${CMAKE_ARGS} ${BUILD_MACHINE}"
echo "STATUS: Run build script with command:"
echo "STATUS: ${BUILD_COMMAND} >& ${REDIRECT_FILE}"
${BUILD_COMMAND} >& ${REDIRECT_FILE}
set EXIT_STATUS = $?

## if an error has occured
if( ${EXIT_STATUS} != 0 ) then
  set LD_ERROR = `grep "ld returned" ${REDIRECT_FILE} | wc -l`

  if ( ${LD_ERROR} != 0 ) then
    echo "ERROR: Linking failed."
    exit 2
  else
    echo "ERROR: Compilation failed."
    exit 3
  endif
endif

# run examples; no more than one at a time, since some of the example apps use/test threads
set COMMAND = "perl ./scripts/nightly_build_helper/bcl_run_individual_examples.pl ${RESULTS_DIR} ${EXECUTABLE} 1"
echo "STATUS: Parsing and running examples with command:"
if( "$RUN_MACHINE" != "" ) then
  echo "STATUS: testing whether $RUN_MACHINE is accessible"
  # try ssh-ing into the machine once and just doing an env command
  # macs sometimes sleep-through the first attempt to ssh onto them, so this command may time out, but will allow subsequent calls to ssh to succeed
  ssh -Y $RUN_MACHINE "sleep 1"
  
  # try ssh-ing again to test whether the machine is accessible
  ssh -Y $RUN_MACHINE "cd '${CHECKOUT_DIR}'"
  if( $? != 0 ) then
    echo "ERROR: ${RUN_MACHINE} is inaccessible or does not have access to ${CHECKOUT_DIR}"
    exit 5
  endif
  ssh -Y ${RUN_MACHINE} "cd ${CHECKOUT_DIR}; ${COMMAND}"
else
  echo "STATUS: ${COMMAND}"
  ${COMMAND}
endif

set EXAMPLE_ERRORS = `cat ${RESULTS_DIR}/example_errors.txt`
echo "STATUS: examples finished, #Errors = ${EXAMPLE_ERRORS}"

#if some examples failed, return 4
if( ${EXAMPLE_ERRORS} != "0") then
  exit 4
endif

exit 0
