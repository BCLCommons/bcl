CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplecluster)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_dendrogram.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_distances_euclidean.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_distances_stored.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_input_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_input_pairwise_list.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_input_table.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_linkage_average.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_linkage_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_linkage_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_linkage_single.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_linkage_total.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_node_colorer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_node_description_average.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_node_description_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_output_centers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_output_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_output_pymol.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_output_pymol_label_small_molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_output_pymol_label_string.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_cluster_output_rows.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
