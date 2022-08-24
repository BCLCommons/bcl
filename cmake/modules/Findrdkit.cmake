# Try to find rdkit
#  rdkit_FOUND          - True if rdkit found.
#  rdkit_INCLUDE_DIRS   - where to find rdkit.h, etc.
#  rdkit_LIBRARIES      - List of libraries when using rdkit.
#  rdkit_VERSION_STRING - The version of rdkit found (x.y.z)

INCLUDE( LibFindMacros)

# this package
SET( rdkit_PACKAGE_NAME rdkit)
SET( rdkit_PACKAGE_VERSION 1.2022.09)

# include directory
FIND_PATH( ${rdkit_PACKAGE_NAME}_INCLUDE_DIR NAMES rdkit.h)

# library
FIND_LIBRARY( ${rdkit_PACKAGE_NAME}_LIBRARY NAMES rdkit)

# version
IF( ${rdkit_PACKAGE_NAME}_INCLUDE_DIR AND EXISTS "${${rdkit_PACKAGE_NAME}_INCLUDE_DIR}/rdkit.h")
  FILE( READ "${${rdkit_PACKAGE_NAME}_INCLUDE_DIR}/rdkit.h" RDKIT_H)
  STRING( REGEX REPLACE ".*#define RDKIT_VERSION \"([0-9]+)\\.([0-9]+)\\.([0-9]+)\".*" "\\1.\\2.\\3" ${rdkit_PACKAGE_NAME}_VERSION_STRING "${RDKIT_H}")
ENDIF()

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${rdkit_PACKAGE_NAME}_PROCESS_INCLUDES ${rdkit_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${rdkit_PACKAGE_NAME}_PROCESS_LIBS ${rdkit_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${rdkit_PACKAGE_NAME})

# runtime library
IF( MINGW)
  FIND_LIBRARY( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES rdkit1.dll)
ELSEIF( APPLE)
  FIND_LIBRARY( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES rdkit.dylib)
ELSEIF( UNIX)
  INCLUDE( MacroGetSoName)
  SET( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY ${${rdkit_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${rdkit_PACKAGE_NAME}")
  MACRO_GET_SO_NAME( ${${rdkit_PACKAGE_NAME}_LIBRARY})
  SET( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${rdkit_PACKAGE_NAME}")
  UNSET( SO_FILE_NAME)
ENDIF()
