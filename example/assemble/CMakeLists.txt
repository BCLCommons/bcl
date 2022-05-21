CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexampleassemble)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_aa_neighbor_list.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_aa_neighbor_list_container.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_aa_neighbor_list_container_generator_protein_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_analyze_chi_angle_pair_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_analyze_chi_angle_recovery.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_analyze_protein_ensemble_aa_neighborhood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_biomolecule.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_chain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_chain_multiplier.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_aa_specified.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_aa_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_all_possible_domains.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_common_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_protein_model_conformation_by_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_protein_model_conformations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_sheet.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_sse_paired.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_sse_size.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_sse_unpaired.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_topology_combined.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_collector_topology_sheet.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_fold_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_fold_template_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_atom.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_atom_coordinates_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_chain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_domain_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_domain_specified.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_domain_sse_pool_overlapping.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sse_from_protein_model_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sse_furthest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sse_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sse_terminus_residue.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sse_unpaired.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sses_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_locator_sub_domain_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_pick_protein_model_conformation_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_pick_sse_furthest_euclidean.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_pick_sse_furthest_euclidean_center.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_pick_sse_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_pick_sse_short_loops.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_pick_sses_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_printer_protein_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_printer_protein_model_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_printer_protein_model_movie.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_printer_protein_model_multimer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_printer_tracker_history.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model_inverter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model_moment_of_inertia.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model_multiplier.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model_with_cache.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_model_with_mutations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_storage_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_with_cache_dataset_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_with_cache_storage_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_protein_with_mutations_dataset_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_quality.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sheet_template_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_compare_extent.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_compare_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_factories.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_factory_conformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_factory_mc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_packer_all_fragment_pairs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_packer_best_fragment_pair.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_packer_best_fragment_pairs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_packing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_packing_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_packing_criteria.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_geometry_phi_psi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pair_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool_agreement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool_insert_coil_into_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool_join_sses.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool_move_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool_mutate_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_pool_split_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_sse_transformer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_topology.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_topology_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_assemble_voxel_grid_aa.cpp
)
ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
	ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
