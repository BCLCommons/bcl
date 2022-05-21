# Try to find ati stream sdk
# Once done this will define
#  opencl_FOUND        - system has ati stream sdk
#  opencl_INCLUDE_DIRS - the ati include directory
#  opencl_LIBRARIES    - Link these to use ati

INCLUDE( LibFindMacros)

# this package
SET( opencl_PACKAGE_NAME opencl)
SET( opencl_PACKAGE_VERSION ${opencl_FIND_VERSION})

# include directory
FIND_PATH( ${opencl_PACKAGE_NAME}_INCLUDE_DIR NAMES CL/cl.hpp )
FIND_LIBRARY(${opencl_PACKAGE_NAME}_LIBRARIES OpenCL DOC "OpenCL lib for OSX")
FIND_PATH(${opencl_PACKAGE_NAME}_INCLUDE_DIRS OpenCL/cl.h DOC "Include for OpenCL on OSX")
    
# library
FIND_LIBRARY( ${opencl_PACKAGE_NAME}_LIBRARY NAMES OpenCL) 

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${opencl_PACKAGE_NAME}_PROCESS_INCLUDES ${opencl_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${opencl_PACKAGE_NAME}_PROCESS_LIBS ${opencl_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${opencl_PACKAGE_NAME})
