# Try to find openblas stream sdk
# Once done this will define
#  openblas_FOUND        - system has openblas stream sdk
#  openblas_INCLUDE_DIRS - the openblas include directory
#  openblas_LIBRARIES    - Link these to use openblas

INCLUDE( LibFindMacros)

# this package
SET( openblas_PACKAGE_NAME openblas)
SET( openblas_PACKAGE_VERSION ${openblas_FIND_VERSION})

# include directory
FIND_PATH( ${openblas_PACKAGE_NAME}_INCLUDE_DIR NAMES OpenBlas/cblas.h)

# library
FIND_LIBRARY( ${openblas_PACKAGE_NAME}_LIBRARY NAMES openblas) 

MARK_AS_ADVANCED( ${openblas_PACKAGE_NAME}_INCLUDE_DIR ${openblas_PACKAGE_NAME}_LIBRARY)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${openblas_PACKAGE_NAME}_PROCESS_INCLUDES ${openblas_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${openblas_PACKAGE_NAME}_PROCESS_LIBS ${openblas_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${openblas_PACKAGE_NAME})

# IF( MINGW)
#  FIND_LIBRARY( ${openblas_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES OpenBlas)
#ELSEIF( UNIX)
IF( UNIX)
  INCLUDE( MacroGetSoName)
  SET( ${openblas_PACKAGE_NAME}_RUNTIME_LIBRARY ${${openblas_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${openblas_PACKAGE_NAME}")
  MACRO_GET_SO_NAME( ${${openblas_PACKAGE_NAME}_LIBRARY})
  SET( ${openblas_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${openblas_PACKAGE_NAME}")
  UNSET( SO_FILE_NAME)
ENDIF()

#IF( ${openblas_PACKAGE_NAME}_FOUND AND ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
#  ADD_DEFINITIONS( -D_WIN32)
#ENDIF()