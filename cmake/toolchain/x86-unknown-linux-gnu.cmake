# this one is important
SET( CMAKE_SYSTEM_NAME Linux)
SET( CMAKE_SYSTEM_PROCESSOR x86)
# this one not so much
SET( CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler
SET( CMAKE_C_COMPILER   "gcc")
SET( CMAKE_CXX_COMPILER "g++")

# where is the target environment 
#SET( CMAKE_FIND_ROOT_PATH /blue/meilerlab/apps/Linux2/x86_64/mingw-w64/2010.09.07/x86_64-w64-mingw32/)

# search for programs in the build host directories
#SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)

# for libraries and headers in the target directories
#SET( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
#SET( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
