CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplerandom)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_random_histogram_1d_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_random_histogram_2d_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_random_uniform_distribution.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
