# This file is only intended to be used to do the auto-configuration of various files for the Rosetta build
# It is *not* intended to be invoked for regular BCL compilation


# Minimal config to get things started.

INCLUDE( "../../bcl_version.cmake")
INCLUDE( "../../cmake/modules/BCLVersion.cmake" )
OPTION( BCL_LICENSE "compile executables with license information, e.g. an expiration time added to bclcommons published binaries" OFF)
OPTION( BCL_RELEASE "compile executables for release to " OFF)
SET( BCL_INSTALL_PREFIX "./")

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/../../source/bcl_version.cpp.in ${CMAKE_CURRENT_SOURCE_DIR}/../../source/bcl_version.cpp)
