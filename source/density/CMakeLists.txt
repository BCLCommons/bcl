CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bcldensity)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_connectivity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_fit_protein_minimizer_mc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_fit_protein_minimizer_powell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_fit_protein_minimizers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_map.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_map_cylindrical.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_mask_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_protein_agreement_ccc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_protein_agreement_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_protein_agreements.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_simulate_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_simulate_gaussian_sphere.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_simulate_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_density_simulators.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
