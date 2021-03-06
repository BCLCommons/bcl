CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplesched)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sched.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sched_function_job_with_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sched_jobs_with_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sched_mutex.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sched_serial_scheduler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sched_sum_function.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
