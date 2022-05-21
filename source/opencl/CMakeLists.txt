CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

IF( NOT DEFINED ${OPENCL_VAR}_FOUND OR NOT ${OPENCL_VAR}_FOUND)
  RETURN()
ENDIF()

SET( LIBRARY_NAME bclopencl)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})
  
INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")
  
SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_approximator_resilient_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_approximator_sequential_minimial_optimization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_approximator_simple_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_buffer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_command_queue.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_context.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_density_fit_protein_minimizer_powell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_density_simulate_gaussian_sphere.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_device.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_extension_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_extensions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_feature_similarity_measures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_insertion_sort.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_kappa_nearest_neighbor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_kernel_source_alternative.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_kernel_source_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_kernel_source_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_kernel_source_string.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_kernel_sources.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_matrix3x3.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_neural_network.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_operations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_platform.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_protein_agreement_ccc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_quality_gdt.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_quality_lcs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_quality_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_saxs_debye.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_support_vector_machine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_opencl_tools.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
