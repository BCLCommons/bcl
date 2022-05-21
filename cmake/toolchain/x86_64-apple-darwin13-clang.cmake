INCLUDE(CMakeForceCompiler)
# this one is important
SET( CMAKE_SYSTEM_NAME Darwin)
SET( CMAKE_SYSTEM_VERSION 10.9)
SET( CMAKE_SYSTEM_PROCESSOR x86_64)
# this one not so much
#SET( DARWIN_MAJOR_VERSION 10)
#SET( DARWIN_MINOR_VERSION 9)
SET( TARGET_PLATFORM ${CMAKE_SYSTEM_PROCESSOR}-apple-darwin13)

# specify the cross compiler
INCLUDE( MacroWrapCompiler)

  MACRO_WRAP_COMPILER( CMAKE_C_COMPILER ${TARGET_PLATFORM}-clang)
  MACRO_WRAP_COMPILER( CMAKE_CXX_COMPILER ${TARGET_PLATFORM}-clang++)
  SET( COMPILER_BASENAME ${TARGET_PLATFORM}-clang++)

GET_FILENAME_COMPONENT( COMPILER_ABS_PATH ${${COMPILER_BASENAME}_PATH} REALPATH)
GET_FILENAME_COMPONENT( COMPILER_DIR_ABS_PATH ${COMPILER_ABS_PATH} PATH)

SET(
  APPLE_SDK_PATH
  ${COMPILER_DIR_ABS_PATH}/../SDK/MacOSX${CMAKE_SYSTEM_VERSION}.sdk
  CACHE
  INTERNAL
  "Location of apple SDKs"
)
SET( 
  APPLE_FRAMEWORKS_PATH 
  ${APPLE_SDK_PATH}/System/Library/Frameworks
  ${APPLE_SDK_PATH}/usr/include
  CACHE INTERNAL "Location of apple frameworks and related header files"
)

#SET( CMAKE_OSX_SYSROOT "/ebio/abt1/nwoetzel/src/osxcross/target/")

# where is the target environment 
#SET( OSX_DEVELOPER_ROOT "/ebio/abt1/nwoetzel/src/osxcross/target/")
#SET( CMAKE_OSX_SYSROOT ${OSX_DEVELOPER_ROOT})
#SET( CMAKE_FIND_ROOT_PATH ${OSX_DEVELOPER_ROOT} ${OSX_DEVELOPER_ROOT}/bin ${OSX_DEVELOPER_ROOT}/usr/bin)
#INCLUDE_DIRECTORIES( ${OSX_DEVELOPER_ROOT}/usr/lib/gcc/x86_64-apple-darwin11/4.2.1/include/)
#INCLUDE_DIRECTORIES( ${OSX_DEVELOPER_ROOT}/usr/include/c++/4.2.1)
#LINK_DIRECTORIES( ${OSX_DEVELOPER_ROOT}/usr/lib/gcc/x86_64-apple-darwin11/4.2.1/)
#ADD_DEFINITIONS( "-mmacosx-version-min=10.6")
#SET( CMAKE_EXE_LINKER_FLAGS "-mmacosx-version-min=10.6")
#SET( CMAKE_SHARED_LINKER_FLAGS "-mmacosx-version-min=10.6")
#SET( CMAKE_MODULE_LINKER_FLAGS "-mmacosx-version-min=10.6")

# search for programs in the build host directories
#SET( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER) 

# for libraries and headers in the target directories
#SET( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
#SET( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
