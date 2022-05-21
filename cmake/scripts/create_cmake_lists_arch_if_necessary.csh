#!/bin/tcsh

set BUILD_TYPES_ALIAS = ( debug release )
set TOOLCHAINS = ( x86_64-unknown-linux-gnu x86_64-w64-mingw32 x86_64-apple-darwin8 )
set SHORT_NAMES = ( linux64 win32 apple64 )

# assumes that this script runs within bcl root
set CHECKOUT_DIR = `pwd`

# update the cmake lists
echo "Updating cmake lists"
${CHECKOUT_DIR}/cmake/scripts/CheckCmakeLists.py ./ -o

# iterate through toolchains
@ TOOLCHAIN_INDEX = 1
while ( $TOOLCHAIN_INDEX <= $#TOOLCHAINS )
  set TOOLCHAIN = $TOOLCHAINS[$TOOLCHAIN_INDEX]
  set SHORT_NAME = $SHORT_NAMES[$TOOLCHAIN_INDEX]

  # iterate through all build types
  @ BUILD_TYPE_INDEX = 1
  while ( $BUILD_TYPE_INDEX <= $#BUILD_TYPES_ALIAS )
    set BUILD_TYPE_ALIAS = $BUILD_TYPES_ALIAS[$BUILD_TYPE_INDEX]
    
    # create Make files
    tcsh ${CHECKOUT_DIR}/cmake/scripts/bcl_recreate_cmake_lists_arch.csh -d ${CHECKOUT_DIR} -t ${TOOLCHAIN} -a ${SHORT_NAME} -b ${BUILD_TYPE_ALIAS}
    @ BUILD_TYPE_INDEX++
  end
  @ TOOLCHAIN_INDEX++
end 

echo "Completed building all cmake files at "`date`
