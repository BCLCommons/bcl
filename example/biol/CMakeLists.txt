CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplebiol)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_back_bone.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_back_bone_completer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_base.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_ca_cb.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_sequence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_sequence_factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_sequence_flexibility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_sequence_phi_psi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_side_chain_factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_aa_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_align_by_aa_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_align_by_pdb_id.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_atom.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_atom_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_atom_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_blast_profile.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_blast_profile_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_chi_angle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_dssp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_environment_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_environment_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_exposure_prediction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_membrane.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_protein_charge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_protein_mutation_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_protein_params.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_ramachandran.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_rotamer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_rotamer_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_ss_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_biol_ss_types.cpp
)
ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
