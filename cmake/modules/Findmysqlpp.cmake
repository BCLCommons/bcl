# Try to find mysqlpp
# Once done this will define
#  mysql_FOUND         - system has mysql
#  mysql_INCLUDE_DIR   - the mysql include directory
#  mysql_LIBRARIES     - the libraries to link for mysql
#  mysqlpp_FOUND       - system has mysqlpp
#  mysqlpp_INCLUDE_DIR - the mysqlpp include directory
#  mysqlpp_LIBRARIES   - Link these to use mysqlpp

INCLUDE( LibFindMacros)

# this package
SET( mysqlpp_PACKAGE_NAME mysqlpp)
SET( mysqlpp_PACKAGE_VERSION ${mysqlpp_FIND_VERSION})

# dependencies
SET( mysql_PACKAGE_NAME mysql)

# include directory
FIND_PATH( ${mysqlpp_PACKAGE_NAME}_INCLUDE_DIR NAMES mysql++.h)

# library
FIND_LIBRARY( ${mysqlpp_PACKAGE_NAME}_LIBRARY NAMES mysqlpp) 

MARK_AS_ADVANCED( ${mysqlpp_PACKAGE_NAME}_INCLUDE_DIR ${mysqlpp_PACKAGE_NAME}_LIBRARY)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${mysqlpp_PACKAGE_NAME}_PROCESS_INCLUDES ${mysqlpp_PACKAGE_NAME}_INCLUDE_DIR ${mysql_PACKAGE_NAME}_INCLUDE_DIRS)
SET( ${mysqlpp_PACKAGE_NAME}_PROCESS_LIBS ${mysqlpp_PACKAGE_NAME}_LIBRARY ${mysql_PACKAGE_NAME}_LIBRARIES)
LIBFIND_PROCESS( ${mysqlpp_PACKAGE_NAME})

# if the mysql library was found, then mysql needs to be removed from the process li

# runtime library
IF( MINGW)
  IF( ${mysqlpp_PACKAGE_NAME}_FOUND)
    SET( ${mysqlpp_PACKAGE_NAME}_RUNTIME_LIBRARY ${${mysqlpp_PACKAGE_NAME}_LIBRARY})
  ENDIF()
ELSEIF( APPLE)
  IF( ${mysqlpp_PACKAGE_NAME}_FOUND)
    SET( ${mysqlpp_PACKAGE_NAME}_RUNTIME_LIBRARY ${${mysqlpp_PACKAGE_NAME}_LIBRARY})
    # only the mysqlpp runtime libraries needs to be added to the externally linked libraries,
    # so set the mysqlpp_LIBRARIES to the empty string, otherwise the
    # mysql libraries are linked twice, causing ld to emit warnings, and potentially other problems
    SET( ${mysqlpp_PACKAGE_NAME}_LIBRARIES ${${mysqlpp_PACKAGE_NAME}_LIBRARY})
  ENDIF()
ELSEIF( UNIX)
  IF( ${mysqlpp_PACKAGE_NAME}_FOUND)
    INCLUDE( MacroGetSoName)
    SET( ${mysqlpp_PACKAGE_NAME}_RUNTIME_LIBRARY ${${mysqlpp_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${mysqlpp_PACKAGE_NAME}")
    MACRO_GET_SO_NAME( ${${mysqlpp_PACKAGE_NAME}_LIBRARY})
    SET( ${mysqlpp_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${mysqlpp_PACKAGE_NAME}")
    UNSET( SO_FILE_NAME)
    # only the mysqlpp runtime libraries needs to be added to the externally linked libraries,
    # so set the mysqlpp_LIBRARIES to the empty string, otherwise the
    # mysql libraries are linked twice, causing ld to emit warnings, and potentially other problems
    SET( ${mysqlpp_PACKAGE_NAME}_LIBRARIES ${${mysqlpp_PACKAGE_NAME}_LIBRARY})
  ENDIF()
ENDIF()
