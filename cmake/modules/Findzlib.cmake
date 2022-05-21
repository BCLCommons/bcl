# Try to find zlib
#  zlib_FOUND          - True if zlib found.
#  zlib_INCLUDE_DIRS   - where to find zlib.h, etc.
#  zlib_LIBRARIES      - List of libraries when using zlib.
#  zlib_VERSION_STRING - The version of zlib found (x.y.z)

INCLUDE( LibFindMacros)

# this package
SET( zlib_PACKAGE_NAME zlib)
SET( zlib_PACKAGE_VERSION 1.2.5)

# include directory
FIND_PATH( ${zlib_PACKAGE_NAME}_INCLUDE_DIR NAMES zlib.h)

# library
FIND_LIBRARY( ${zlib_PACKAGE_NAME}_LIBRARY NAMES z zlib zlib1 zdll)

# version
IF( ${zlib_PACKAGE_NAME}_INCLUDE_DIR AND EXISTS "${${zlib_PACKAGE_NAME}_INCLUDE_DIR}/zlib.h")
  FILE( READ "${${zlib_PACKAGE_NAME}_INCLUDE_DIR}/zlib.h" ZLIB_H)
  STRING( REGEX REPLACE ".*#define ZLIB_VERSION \"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\1.\\2.\\3" ${zlib_PACKAGE_NAME}_VERSION_STRING "${ZLIB_H}")
ENDIF()

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${zlib_PACKAGE_NAME}_PROCESS_INCLUDES ${zlib_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${zlib_PACKAGE_NAME}_PROCESS_LIBS ${zlib_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${zlib_PACKAGE_NAME})

# runtime library
IF( MINGW)
  FIND_LIBRARY( ${zlib_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES zlib1.dll)
ELSEIF( APPLE)
  FIND_LIBRARY( ${zlib_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES libz.dylib)
ELSEIF( UNIX)
  INCLUDE( MacroGetSoName)
  SET( ${zlib_PACKAGE_NAME}_RUNTIME_LIBRARY ${${zlib_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${zlib_PACKAGE_NAME}")
  MACRO_GET_SO_NAME( ${${zlib_PACKAGE_NAME}_LIBRARY})
  SET( ${zlib_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${zlib_PACKAGE_NAME}")
  UNSET( SO_FILE_NAME)
ENDIF()
