# 32-bit macs (intel core solo or duo)
# this one is important
SET( CMAKE_SYSTEM_NAME Darwin)
SET( CMAKE_SYSTEM_VERSION 10.4)
SET( CMAKE_SYSTEM_PROCESSOR i686)
# this one not so much
SET( TARGET_PLATFORM ${CMAKE_SYSTEM_PROCESSOR}-apple-darwin8)

# determine cross-compiler location and write the compiler wrapper
INCLUDE( MacroWrapCompiler)
MACRO_WRAP_COMPILER( CMAKE_C_COMPILER ${TARGET_PLATFORM}-gcc)
MACRO_WRAP_COMPILER( CMAKE_CXX_COMPILER ${TARGET_PLATFORM}-g++)

# where is the target environment 
#SET( CMAKE_FIND_ROOT_PATH /blue/meilerlab/apps/Linux2/x86_64/apple-darwin/2011.03.09/)
SET(
  APPLE_SDK_PATH
  /blue/meilerlab/apps/Linux2/x86_64/apple-darwin/2011.03.09/SDKs/MacOSX10.4u.sdk
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
