CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplecoord)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_axes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_cyclic_coordinate_descent.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_cylinder_coordinates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_geometric_hash_storage_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_geometric_hash_storage_hash_map.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_geometric_hashing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_line_segment_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_line_segment_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_moment_of_inertia.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_movable_eccentric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_combine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_rotate_defined.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_rotate_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_rotate_random_external_reference.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_transform_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_translate_defined.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_translate_external_axis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_move_translate_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_point_cloud.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_point_to_key_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_point_to_key_spherical_radius.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_polygon.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_coord_sphere.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
