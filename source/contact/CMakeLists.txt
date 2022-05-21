CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclcontact)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_aa_correlation_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_ann.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_calculate_correlations_mi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_calculate_correlations_sm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_correlation_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_correlation_storage_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_map.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_order.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_prediction_map.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_recovery.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_sse_prediction_map.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_statistics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_contact_types.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
