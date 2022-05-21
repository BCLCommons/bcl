CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplegraph)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_common_subgraph_isomorphism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_common_subgraph_isomorphism_base.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_const_graph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_csi_substructure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_edge_cover_ring_perception.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_edge_with_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_exhaustive_ring_perception.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_graph_with_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_path.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_ring.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_subgraph.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_subgraph_isomorphism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_tree_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_undirected_edge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_graph_vertex_with_data.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
