#!/bin/tcsh

# preparing make file
tcsh ${BCL_DIRECTORY}/cmake/scripts/bcl_recreate_cmake_lists_arch.csh -d ${BCL_DIRECTORY} 

# check for error
if ( $status != 0) then
  echo "ERROR: error while updating Makefiles, cannot build"
  exit -1
endif

# preparing build
set NUMBER_JOBS = 8
set PUMP_COMMAND
set TIME_CMD

# only use pump if distcc is available
which distcc > &/dev/null
if( $status == 0) then
  set PUMP_COMMAND = pump
  set NUMBER_JOBS = `distcc -j`
endif

# use /usr/bin/time (usually gnu time) if possible
/usr/bin/time -v ls > &/dev/null
if( $status == 0) then
  set TIME_COMMAND = "/usr/bin/time -v"
else
  # otherwise, use the builtin time
  set TIME_COMMAND = "time"
endif

# check if this is the correct directory
set CURRENT_DIR = `pwd`
set BUILD_DIR = ${TOOLCHAIN_ALIAS}_${BUILD_TYPE_ALIAS}

if( ! (${CURRENT_DIR} =~ ${BUILD_DIR})) then
  cd ${BCL_DIRECTORY}/build/${TOOLCHAIN_ALIAS}_${BUILD_TYPE_ALIAS}
endif

# executing make with or without pump
echo "${TIME_COMMAND} ${PUMP_COMMAND} make $* -j ${NUMBER_JOBS}"
${TIME_COMMAND} ${PUMP_COMMAND} make $* -j ${NUMBER_JOBS}
