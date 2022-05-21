# - Try to find bzip2
# Once done this will define
#
#  bzip2_FOUND        - system has bzip2
#  bzip2_INCLUDE_DIRS - the bzip2 include directory
#  bzip2_LIBRARIES    - Link these to use bzip2
#  bzip2_NEED_PREFIX  - this is set if the functions are prefixed with BZ2_

INCLUDE( LibFindMacros)

# this package
SET( bzip2_PACKAGE_NAME bzip2)
SET( bzip2_PACKAGE_VERSION ${bzip2_FIND_VERSION})

# include directory
FIND_PATH( ${bzip2_PACKAGE_NAME}_INCLUDE_DIR NAMES bzlib.h )

# library
FIND_LIBRARY( ${bzip2_PACKAGE_NAME}_LIBRARY NAMES bz2 bzip2 bz2dll)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${bzip2_PACKAGE_NAME}_PROCESS_INCLUDES ${bzip2_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${bzip2_PACKAGE_NAME}_PROCESS_LIBS ${bzip2_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${bzip2_PACKAGE_NAME})

#IF( ${bzip2_PACKAGE_NAME}_FOUND)
#   INCLUDE( CheckLibraryExists)
#   CHECK_LIBRARY_EXISTS( ${${bzip2_PACKAGE_NAME}_LIBRARIES} BZ2_bzCompressInit "" BZIP2_NEED_PREFIX)
#ENDIF( ${bzip2_PACKAGE_NAME}_FOUND)

# runtime library
IF( MINGW)
  FIND_LIBRARY( ${bzip2_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES libbz2.dll)
ELSEIF( APPLE)
  FIND_LIBRARY( ${bzip2_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES libbz2.dylib)
ELSEIF( UNIX)
  INCLUDE( MacroGetSoName)
  SET( ${bzip2_PACKAGE_NAME}_RUNTIME_LIBRARY ${${bzip2_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${bzip2_PACKAGE_NAME}")
  MACRO_GET_SO_NAME( ${${bzip2_PACKAGE_NAME}_LIBRARY})
  SET( ${bzip2_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${bzip2_PACKAGE_NAME}")
  UNSET( SO_FILE_NAME)
ENDIF()
