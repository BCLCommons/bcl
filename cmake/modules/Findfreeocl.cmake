# Try to find freeocl stream sdk
# Once done this will define
#  freeocl_FOUND        - system has freeocl stream sdk
#  freeocl_INCLUDE_DIRS - the freeocl include directory
#  freeocl_LIBRARIES    - Link these to use freeocl

INCLUDE( LibFindMacros)

# this package
SET( freeocl_PACKAGE_NAME freeocl)
SET( freeocl_PACKAGE_VERSION ${freeocl_FIND_VERSION})

# include directory
FIND_PATH( ${freeocl_PACKAGE_NAME}_INCLUDE_DIR NAMES CL/cl.hpp)

# library
FIND_LIBRARY( ${freeocl_PACKAGE_NAME}_LIBRARY NAMES OpenCL) 

MARK_AS_ADVANCED( ${freeocl_PACKAGE_NAME}_INCLUDE_DIR ${freeocl_PACKAGE_NAME}_LIBRARY)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${freeocl_PACKAGE_NAME}_PROCESS_INCLUDES ${freeocl_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${freeocl_PACKAGE_NAME}_PROCESS_LIBS ${freeocl_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${freeocl_PACKAGE_NAME})

IF( MINGW)
  FIND_LIBRARY( ${freeocl_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES OpenCL)
ELSEIF( UNIX)
  INCLUDE( MacroGetSoName)
  SET( ${freeocl_PACKAGE_NAME}_RUNTIME_LIBRARY ${${freeocl_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${freeocl_PACKAGE_NAME}")
  MACRO_GET_SO_NAME( ${${freeocl_PACKAGE_NAME}_LIBRARY})
  SET( ${freeocl_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${freeocl_PACKAGE_NAME}")
  UNSET( SO_FILE_NAME)
ENDIF()

#IF( ${freeocl_PACKAGE_NAME}_FOUND AND ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
#  ADD_DEFINITIONS( -D_WIN32)
#ENDIF()