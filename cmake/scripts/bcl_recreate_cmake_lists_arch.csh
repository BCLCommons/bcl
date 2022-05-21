#!/bin/tcsh

# @author Jeff Mendenhall, Nils Woetzel
# @date 02/06/2012
# this script generates a folder and the cmake files according to a given toolchain and build type

# list of valid build types and their aliases
set BUILD_TYPES = ( Debug Release )
set BUILD_TYPES_ALIAS = ( debug release )

@ argnr = 1
while( $argnr <= $#argv)
  switch ($argv[$argnr])
    case -d:
    case --directory:
      @ argnr++
      if ( $argnr > $#argv) then
        echo "ERROR: missing parameter for -d;--directory!" 
        exit -1
      endif
      set CHECKOUT_DIR = $argv[$argnr]
      breaksw
    case -t:
    case --toolchain:
      @ argnr++
      if ( $argnr > $#argv) then
        echo "ERROR: missing parameter for -t;--toolchain!" 
        exit -1
      endif
      set TOOLCHAIN = $argv[$argnr]
      breaksw
    case -b:
    case --buildtype:
      @ argnr++
      if ($argnr > $#argv) then
        echo "ERROR: missing parameter for -b;--buildtype!" 
        exit -1
      endif
      set BUILD_TYPE_ALIAS = $argv[$argnr]
      breaksw
    case -a:
    case --alias:
      @ argnr++
      if ( $argnr > $#argv) then
        echo "ERROR: missing parameter for -a;--alias!" 
        exit -1
      endif
      set TOOLCHAIN_ALIAS = "$argv[$argnr]"
      breaksw
    case -r:
    case --recreate:
      set RECREATE = true
      breaksw
    case -h:
    case --help:
      goto usageblock
      breaksw
  endsw
  @ argnr++;
end

# checkout dir
if( ! $?CHECKOUT_DIR) then
  set CHECKOUT_DIR = `pwd`
  echo "no directory given; assuming the work directory is the checkout directory: ${CHECKOUT_DIR}"
endif

# always add the compiler wrapper
setenv LIBRARY_PATH /usr/lib/x86_64-linux-gnu

#sanity checks
if( ! $?TOOLCHAIN) then
  echo "ERROR: no toolchain given"
  goto usageblock
endif

if( ! $?BUILD_TYPE_ALIAS) then
  echo "ERROR: no buildtype given"
  goto usageblock
endif

if( ! $?TOOLCHAIN_ALIAS) then
  set TOOLCHAIN_ALIAS = $TOOLCHAIN
endif

set TOOLCHAIN_FILE = ${CHECKOUT_DIR}/cmake/toolchain/${TOOLCHAIN}.cmake

# check that the toolchain file exists
if( ! -e ${TOOLCHAIN_FILE} ) then
  if( -d ${CHECKOUT_DIR}/cmake/toolchain/ ) then
    echo "ERROR: ${TOOLCHAIN} is not a valid toolchain. Valid toolchains are: " 
    ls ${CHECKOUT_DIR}/cmake/toolchain/ | grep ".cmake" | sed 's/\.cmake//'
  else
    echo "ERROR: $0 must be started from the bcl directory."
  endif
  goto usageblock
endif

# build type from alias
@ argnr = 1
while ( $argnr <= $#BUILD_TYPES_ALIAS)
  if ( $BUILD_TYPES_ALIAS[$argnr] == ${BUILD_TYPE_ALIAS}) then
    set BUILD_TYPE = $BUILD_TYPES[$argnr]
  endif
  @ argnr++
end

# check if alias was found - only then a valid build type was supplied
if( ! $?BUILD_TYPE) then
  echo "ERROR: non valid build type supplied: ${BUILD_TYPE_ALIAS} is not in list" ${BUILD_TYPES_ALIAS}
  goto usageblock
endif

cd $CHECKOUT_DIR
set BUILD_DIR = ./build/${TOOLCHAIN_ALIAS}_$BUILD_TYPE_ALIAS
if ( $?RECREATE) then
  echo "Updating cmake lists"
  # update the cmake lists
  ${CHECKOUT_DIR}/cmake/scripts/CheckCmakeLists.py ./ -o

  echo "Removing existing files in $BUILD_DIR"
  rm -rf $BUILD_DIR
endif

# create directoy if not present
if( ! -d $BUILD_DIR) then
  mkdir -p $BUILD_DIR
endif

cd $BUILD_DIR

# calling cmake
if( ! -f Makefile || $?RECREATE) then
  echo "Calling cmake in $BUILD_DIR with "
  echo "cmake -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} -DCMAKE_BUILD_TYPE=$BUILD_TYPE ${CHECKOUT_DIR}"
  cmake -DCMAKE_TOOLCHAIN_FILE=${TOOLCHAIN_FILE} -DCMAKE_BUILD_TYPE=$BUILD_TYPE ${CHECKOUT_DIR}
  echo "Done"
endif

cd $CHECKOUT_DIR
exit( 0)

usageblock:
  echo -n "\nUsage:  [t;--toolchain] {toolchain} [-a;--alias] {toolchain alias} [-b;--buildtype] (${BUILD_TYPES_ALIAS}) [-h;--help]"
  echo "  -t the toolchain to use"
  echo "  -a alias for toolchain - directory name for build"
  echo "  -b build type from list (${BUILD_TYPES_ALIAS})"
  echo "  -h print this help"
  exit 0
