CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclscorestat)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_aa_count.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_aa_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_aa_distance_angle_contacts.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_aa_distance_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_contact_order.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_fold_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_loop_angle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_loop_closure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_loop_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_neighbor_count.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_neighbor_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_ols.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_phipsi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_protein_model_packing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_protein_model_sse_triplet_chirality.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_radius_of_gyration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_sheet_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_side_chain_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_sse_count.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_sse_membrane_alignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_sse_packing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_sspred_agreement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_scorestat_strand_alignment.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
