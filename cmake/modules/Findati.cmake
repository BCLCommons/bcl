# Try to find ati stream sdk
# Once done this will define
#  ati_FOUND        - system has ati stream sdk
#  ati_INCLUDE_DIRS - the ati include directory
#  ati_LIBRARIES    - Link these to use ati

INCLUDE( LibFindMacros)

# this package
SET( ati_PACKAGE_NAME ati)
SET( ati_PACKAGE_VERSION ${ati_FIND_VERSION})

# include directory
FIND_PATH( ${ati_PACKAGE_NAME}_INCLUDE_DIR NAMES CL/cl.hpp)

# library
FIND_LIBRARY( ${ati_PACKAGE_NAME}_LIBRARY NAMES OpenCL) 

MARK_AS_ADVANCED( ${ati_PACKAGE_NAME}_INCLUDE_DIR ${ati_PACKAGE_NAME}_LIBRARY)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${ati_PACKAGE_NAME}_PROCESS_INCLUDES ${ati_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${ati_PACKAGE_NAME}_PROCESS_LIBS ${ati_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${ati_PACKAGE_NAME})

# runtime library - disabled right now, since there is no license that would allow redistribution
#IF( MINGW)
#  FIND_LIBRARY( ${ati_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES OpenCL.dll)
#ELSEIF( UNIX)
#  INCLUDE( MacroGetSoName)
#  SET( ${ati_PACKAGE_NAME}_RUNTIME_LIBRARY ${${ati_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${ati_PACKAGE_NAME}")
#  MACRO_GET_SO_NAME( ${${ati_PACKAGE_NAME}_LIBRARY})
#  SET( ${ati_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${ati_PACKAGE_NAME}")
#  UNSET( SO_FILE_NAME)
#ENDIF()
#
#IF( ${ati_PACKAGE_NAME}_FOUND AND ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
#  ADD_DEFINITIONS( -D_WIN32)
#ENDIF()