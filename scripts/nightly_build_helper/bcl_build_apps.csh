#!/bin/tcsh

# make sure that at least one argument was given and that not too many arguments were given
if( $#argv > 3 ) then
  echo "Usage: $0 [architecture] [results_file] [build_machine]"
  echo "Script should be run from the checkout folder"
  echo "  architecture: uses by default cmake toolchain x86_64-unknown-linux-gnu (only the architecture is needed)"
  echo "  results_file: uses by default ./apps_make_log.txt to output cmake and compilation messages"
  echo "  build_machine: uses by default the local machine that the script is run from"
  echo "Your command line was: $*"
  exit 1
endif

# set the arguments to named variables
set TOOLCHAIN = $1
set REDIRECT_FILE = $2
set BUILD_MACHINE = $3
set CHECKOUT_DIR = `pwd`
set BUILD_TARGET = "bcl-apps-static"

# set default toolchain if none was given
if( $#argv == 0 ) then
  set TOOLCHAIN = "x86_64-unknown-linux-gnu"
endif
	
# set default results folder if non was given
if( $#argv <= 1 ) then
  set REDIRECT_FILE = "${CHECKOUT_DIR}/apps_make_log.txt"
endif
	
# build bcl
set BUILD_COMMAND = "./scripts/nightly_build_helper/bcl_build.csh ${BUILD_TARGET} ${TOOLCHAIN} -DCMAKE_BUILD_TYPE=Release ${BUILD_MACHINE}"
echo "STATUS: Run build script with command:"
echo "STATUS: ${BUILD_COMMAND} > ${REDIRECT_FILE}"
${BUILD_COMMAND} > ${REDIRECT_FILE}
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

exit 0
