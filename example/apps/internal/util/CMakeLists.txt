CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET(
  LIBRARY_SOURCES
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

SET( SUBDIRECTORY_LIBRARY_SOURCES ${LIBRARY_SOURCES} PARENT_SCOPE)
