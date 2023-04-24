CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplechemistry)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_aa_fragment_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_clash_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_configurational_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_conformational_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_constitutional_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_environment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_environment_bender.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_one_four_interaction_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_vdw_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_atom_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_bond_configurational.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_bond_conformational.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_bond_constitutional.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_bond_isometry_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_bond_lengths.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_collector_valence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_configuration_graph_converter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_configuration_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_configuration_set_same_constitution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_configurational_bond_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_configurational_bond_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_dihedral_bins.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_dihedrals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_property.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_real_space_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_substructure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_by_symmetry_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_property_field_correlation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_property_rmsd_x.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_psi_field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_comparison_psi_flex_field.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_graph_converter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_set_same_configuration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_conformation_set_same_constitution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_constitution_graph_converter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_constitution_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_constitutional_bond_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_constitutional_bond_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_coulombic_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_cs_code.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_descriptor_to_score_adaptor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_electron_configuration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_element_structure_factor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_element_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_element_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_configuration_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_conformation_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_connector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_constitution_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_cyclize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_feed_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_feed_from_stream.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_graph_marker.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_grow.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_molecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_probability_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_remove_bond.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_ring_swap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_ecfp_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_isolate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_largest_component.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_linear_fragments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_rigid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_rings.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_rings_with_unsaturated_substituents.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_scaffolds.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_split_unbridged_rings.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_fragment_swap_amide.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_has_properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_merge_fragment_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_mol_align_by_parts.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecular_configuration_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecular_conformation_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecular_constitution_shared.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecule_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecule_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecule_environment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecule_fragment_recombination.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_molecule_storage_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_mutate_chirality.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_mutate_fragment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_pick_atom_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_possible_atom_types_for_atom.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_priority_dihedral_angles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_reaction_complete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_reaction_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_reaction_search.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_reaction_structure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_reaction_worker.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_rotamer_dihedral_bond_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_rotamer_library_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_sample_conformations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_search_fragment_library_from_tree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_fragment_isomorphism.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_fragment_mapping.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_misc_properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_qsar_storage_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_string_properties_cached.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_string_properties_mapped.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_string_properties_name.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_string_properties_numeric.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_small_molecule_string_properties_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_solvents.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_stereocenters_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_sub_fragment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_chemistry_substituent_conformational.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
