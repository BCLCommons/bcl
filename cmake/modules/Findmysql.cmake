# Try to find mysql
# Once done this will define
#  mysql_FOUND - system has MYSQLPP
#  mysql_INCLUDE_DIR - the MYSQLPP include directory
#  mysql_LIBRARIES - Link these to use MYSQLPP

INCLUDE( LibFindMacros)

# this package
SET( mysql_PACKAGE_NAME mysql)
SET( mysql_PACKAGE_VERSION ${mysql_FIND_VERSION})

# include directory
FIND_PATH( ${mysql_PACKAGE_NAME}_INCLUDE_DIR NAMES mysql.h)

# library; note, capitalization of mySQL is important since it must mesh with the actual windows dll in extern/
FIND_LIBRARY( ${mysql_PACKAGE_NAME}_LIBRARY NAMES mysqlclient_r mySQL mysql) 

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${mysql_PACKAGE_NAME}_PROCESS_INCLUDES ${mysql_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${mysql_PACKAGE_NAME}_PROCESS_LIBS ${mysql_PACKAGE_NAME}_LIBRARY)
LIBFIND_PROCESS( ${mysql_PACKAGE_NAME})

# runtime library
IF( MINGW)
  IF( ${mysql_PACKAGE_NAME}_FOUND)
    SET( ${mysql_PACKAGE_NAME}_RUNTIME_LIBRARY ${${mysql_PACKAGE_NAME}_LIBRARY})
    MESSAGE( STATUS 
      "Enabling --enable-stdcall-fixup so that linker doesn't cause make to return non-zero status code"
      "due to Warning: resolving _mysql_free_result@4 by linking to _mysql_free_result"
    )
    # define enable-stdcall-fixup so that warning:
    # Warning: resolving _mysql_free_result@4 by linking to _mysql_free_result
    # Use --enable-stdcall-fixup to disable these warnings
    # Use --disable-stdcall-fixup to disable these fixups 
    # does not cause cmake to return a non-zero error-code
    SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--enable-stdcall-fixup" CACHE INTERNAL "fixup necessary for linking mysql")
    SET( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--enable-stdcall-fixup" CACHE INTERNAL "fixup necessary for linking mysql")
    SET( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--enable-stdcall-fixup" CACHE INTERNAL "fixup necessary for linking mysql")
  ENDIF()
ELSEIF( APPLE)
  FIND_LIBRARY( ${mysql_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES libmysqlclient_r.dylib)
  SET( ${${mysql_PACKAGE_NAME}_LIBRARY} ${mysql_PACKAGE_NAME}_RUNTIME_LIBRARY CACHE INTERNAL "runtime library for ${mysql_PACKAGE_NAME}")
ELSEIF( UNIX)
  INCLUDE( MacroGetSoName)
  SET( ${mysql_PACKAGE_NAME}_RUNTIME_LIBRARY ${${mysql_PACKAGE_NAME}_LIBRARY} CACHE INTERNAL "runtime library for ${mysql_PACKAGE_NAME}")
  MACRO_GET_SO_NAME( ${${mysql_PACKAGE_NAME}_LIBRARY})
  SET( ${mysql_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${mysql_PACKAGE_NAME}")
  UNSET( SO_FILE_NAME)
ENDIF()
