CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclappsinternalchemistry)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_add_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_alchemical_transformation_mapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_align_binding_poses.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_analyze_fld_results.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_build_conformer_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_build_fragment_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_build_rotamer_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_build_scaffold_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_conformer_generator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_evogen_analysis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_extract_fld_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_focused_library_design_old.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_focused_library_design_recombination_old.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_atom_environment_hashmap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_atom_hybridization_descriptors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_hierarchical_tree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_link_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_make_chimeric_molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_mutate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_react.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_quench_reactive_groups.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_react_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_reaction_combichem.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS $<TARGET_OBJECTS:${LIBRARY_OBJECTS}> PARENT_SCOPE)
