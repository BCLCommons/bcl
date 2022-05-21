CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclfold)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_add_parabolic_loops.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_collector_loop_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_collector_loop_domain_all_non_rigid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_collector_loop_domain_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_collector_unconnected_sses.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_default_flags.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_default_mutates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_default_scores.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_handler_locator_loop_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_handler_locator_loop_segment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_locator_loop.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_locator_loop_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_locator_loop_segment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_locator_missing_coordinates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_locator_unconnected_segments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_loop_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_loop_domain_c_to_n.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_loop_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_loop_parameters.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_loop_segment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_loop_segment_sequence_order.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_aa_phi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_aa_psi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_aa_rotate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_aa_sequence_grow.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_aa_set_phi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_aa_set_psi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_domain_flip.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_domain_merge_consecutive_ss_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_domain_shuffle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_domain_sse_pair_trim.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_domain_sse_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_domain_transformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_loop_domain_dihedral.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_membrane.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_membrane_chain_move.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_multimer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_ensemble_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_ensemble_remove.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_add_sheet_from_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_chain_move.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_compress.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_domain_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_filter_conformations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_fix_loop_closure_wrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_grow_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_loop_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_loop_domain_ccd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_loop_domain_grow.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_loop_resize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_move_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_multiple_geometries.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_pair_strands.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_replicate_conformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_add_multiple.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_move.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_pair.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_pair_align_and_pull.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_pair_clash.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_pair_fix_loop_closure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_pair_hinge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_remove.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_resize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_seed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_swap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_swap_body.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_swap_multimer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_swap_with_pool.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_sse_swap_with_pool_overlap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_strand_switch_sheet.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_switch_conformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_protein_model_thread_sequence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_cycle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_divide.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_fit_to_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_order.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_register_fix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_register_shift.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_sort.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sheet_twist.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sse_bend_ramachandran.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sse_bend_random.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sse_bend_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_sse_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutate_tree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_mutation_residue.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_phi_psi_generator_ccd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_phi_psi_generator_ramachandran.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_domain_using_fold_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_sse_distance_restraint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_sse_into_body.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_sse_next_to_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_sse_short_loop.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_placement_strand_next_to_sheet.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protein_geometry.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_assembly.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_create.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_dock.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_em.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_ensemble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_ensemble_filter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_ensemble_replicate_conformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_ensemble_switch_conformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_loop_close.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_loop_coordinate_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_membrane.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_multimer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_refinement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_restraint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocol_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_protocols.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_score_weight_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_scores.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_setup.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_fold_stage_factory.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
