#!/bin/bash

TOOLCHAIN_ALIAS=linux64 
TOOLCHAIN=x86_64-unknown-linux-gnu 
BUILD_TYPE_ALIAS=release
BCL_DIRECTORY=`pwd`
export TOOLCHAIN_ALIAS TOOLCHAIN BUILD_TYPE_ALIAS BCL_DIRECTORY
mkdir -p build/${TOOLCHAIN_ALIAS}_${BUILD_TYPE_ALIAS}

python2 ./cmake/scripts/CheckCmakeLists.py ./ -o  

python2 ./scripts/code/CreateNamespaceForwardHeaders.py ./ -o

if [ -z "${SHARED+x}" ] ; then
  tcsh ./scripts/build/check_pump_make.csh -k static 
else
  tcsh ./scripts/build/check_pump_make.csh -k shared
fi
