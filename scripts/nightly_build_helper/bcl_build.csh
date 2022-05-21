#!/bin/tcsh

# make sure that at least one argument was given and that not too many arguments were given
if( $#argv > 4 ) then
  echo "Usage: $0 [target] [architecture] [cmake arguments] [build_machine]"
  echo "Script should be run from the checkout folder."
  echo "  target: uses by default bcl-example-static"
  echo "  architecture: uses by default cmake toolchain x86_64-unknown-linux-gnu (only the architecture is needed)"
  echo "  cmake arguments: by default, no extra arguments are given to cmake"
  echo "  build_machine: uses by default the local machine that the script is run from"
  exit 1
endif

# set the arguments to named variables and some constant variables
set BUILD_TARGET = $1
set TOOLCHAIN = $2
set CMAKE_ARGS = $3
set BUILD_MACHINE = $4
set CHECKOUT_DIR = `pwd`

# set default toolchain if none was given
if( $#argv == 0 ) then
  set BUILD_TARGET = "bcl-example-static"
  set TOOLCHAIN = "x86_64-unknown-linux-gnu"
else if( $#argv == 1 ) then
  set BUILD_TARGET = "bcl-example-static"
endif

# set variable for the toolchain
set TOOLCHAIN_FILE = ${CHECKOUT_DIR}/cmake/toolchain/${TOOLCHAIN}.cmake
# check that the toolchain existed
if( ! -e ${TOOLCHAIN_FILE} ) then
  if( -d ${CHECKOUT_DIR}/cmake/toolchain/ ) then
    echo "ERROR: ${TOOLCHAIN} is not a valid toolchain. Valid toolchains are: " 
    ls ${CHECKOUT_DIR}/cmake/toolchain/ | grep ".cmake" | sed 's/\.cmake//'
  else
    echo "ERROR: $0 must be started from the bcl directory."
  endif
  exit 1
endif

# print start time
echo -n "STATUS: Start building at: "
date

# create the build directories; do not clean to reuse a previous build; clean manually if needed
set BUILD_DIR = ${PWD}/build/${TOOLCHAIN}
set EXECUTABLE_DIR = "${BUILD_DIR}/bin"
if( ! -d ${BUILD_DIR} ) then
  echo "STATUS: Creating the build directory at ${BUILD_DIR}"
  mkdir -p ${BUILD_DIR}
else
  echo "STATUS: Using existing build directory at ${BUILD_DIR}"
endif

# change directory
cd ${BUILD_DIR}
echo "STATUS: Changing into ${BUILD_DIR}"

# run cmake with the desired toolchain
set CMAKE_COMMAND = "cmake -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} ${CMAKE_ARGS} ${CHECKOUT_DIR}"
echo "STATUS: Start cmake with command:"
echo "STATUS: ${CMAKE_COMMAND}"
${CMAKE_COMMAND}

# check whether cmake succeeded
set EXIT_STATUS = $?
if( ${EXIT_STATUS} != 0 ) then
  exit 1
endif

echo "STATUS: Completed cmake."

# call make on build machine or localhost
set COMPILE_COMMAND = "pump make -k -j80 ${BUILD_TARGET}"
echo "STATUS: Start compiling with command:"
if( $BUILD_MACHINE != "" )  then
  # variables have to kept out of single quotes or they are not replaced by their value
  # single quotes have to be used, because escaping double quotes in double quotes does not work in tcsh
  echo 'STATUS: ssh '${BUILD_MACHINE}' "cd '${BUILD_DIR}'; '${COMPILE_COMMAND}'"'
  ssh ${BUILD_MACHINE} "cd ${BUILD_DIR}; ${COMPILE_COMMAND}"
  # the set exit status command must come immediately after the compile command, otherwise, no error is detected
  set EXIT_STATUS = $?
else
  echo "STATUS: ${COMPILE_COMMAND}"
  ${COMPILE_COMMAND}
  # the set exit status command must come immediately after the compile command, otherwise, no error is detected
  set EXIT_STATUS = $?
endif
	
echo "STATUS: Compiling completed."

# print finish time
echo -n "STATUS: Finished building at: "
date

echo "STATUS: Changing into ${CHECKOUT_DIR}" 
cd ${CHECKOUT_DIR}
exit ${EXIT_STATUS}
