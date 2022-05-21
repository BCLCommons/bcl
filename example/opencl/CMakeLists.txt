CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

IF( NOT DEFINED ${OPENCL_VAR}_FOUND OR NOT ${${OPENCL_VAR}_FOUND})
  RETURN()
ENDIF()

SET( LIBRARY_NAME bclexampleopencl)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_approximator_resilient_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_approximator_sequential_minimial_optimization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_approximator_simple_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_arg_max.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_arg_min.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_dataset_min_max.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_device.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_euclidean_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_extension_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_feature_similarity_measures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_insertion_sort.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_kappa_nearest_neighbor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_matrix_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_matrix_multiply.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_matrix_transpose.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_neural_network.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_platform.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_quality_gdt.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_quality_lcs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_quality_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_saxs_debye.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_singular_value_decomposition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_support_vector_machine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_tools.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_transfer_function_gaussian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_transfer_function_sigmoid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opencl_vector_matrix_add.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
