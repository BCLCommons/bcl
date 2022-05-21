CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclappsmolecule)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_align_to_scaffold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_evolutionary_generator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_extract_mcs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_focused_library_design.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_rosetta_ncaa_instructions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_make_fld_scaffold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_make_grow_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_map_params.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_match_reactions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_atom_pairwise_potential.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_coordinates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_extract_protein_pocket.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_features.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_filter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_fit.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_multialign.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_reorder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_molecule_unique.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_pharm_map.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_set_sample_by_parts_atoms.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS $<TARGET_OBJECTS:${LIBRARY_OBJECTS}> PARENT_SCOPE)
