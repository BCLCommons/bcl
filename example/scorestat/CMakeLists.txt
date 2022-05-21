CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplescorestat)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_loop_closure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_loop_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_neighbor_count.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_neighbor_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_ols.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_phipsi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_protein_model_packing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_protein_model_sse_triplet_chirality.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_scorestat_radius_of_gyration.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
