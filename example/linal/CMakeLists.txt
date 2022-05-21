CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplelinal)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_distance_geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix3x3.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_const_reference.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_inversion_cholesky.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_inversion_gauss_jordan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_inversion_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_inversion_moore_penrose.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_operations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_matrix_reference.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_operations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_principal_component_analysis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_symmetric_eigensolver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_2d_operations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_3d_operations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_const_reference.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_nd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_operations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_linal_vector_reference.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
