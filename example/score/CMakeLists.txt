CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplescore)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_blast_profile.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_blosum.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_identity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_pam.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_phat.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_property.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_assignment_ss_prediction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_neighborhood_distances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_neighborhood_exposure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_neighborhood_exposure_prediction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_atom_clash.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_clash.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_contact.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_contact_energy.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_distance_fitted_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_distance_smooth.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_hi_res_clash.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_pair_sidechain_interaction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_sequence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_aa_sequence_pair.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_accessibility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_accessibility_hydrophobic_moment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_accessibility_hydrophobic_moment_magnitude.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_alignment_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_alignment_gap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_alignment_quality.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_assignment_gap_simple.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_body_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_body_connectivity_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_body_extent_agreement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_body_extent_position_agreement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_consensus_enrichment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_contact_order.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_bipolar.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_coordinate_exclusion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_coordinate_triangulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_data_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_distance_change_magnitude.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_euclidian_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_residue_type_exclusion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_sequence_separation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_size.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_sse_center.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_sse_connection.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_sse_size.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_sse_term.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_data_set_pairwise_structural_exposure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_density_profile_sse_agreement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_energy_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_environment_predictions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_epr_accessibility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_fuzzy_logic_filter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_log_normal_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_loop.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_loop_closure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_phi_psi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_phi_psi_with_sspred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_pofr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_atom_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_aa_neighborhood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_completeness.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_defined_loops.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_gap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_inverted.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_loop_domain_closure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_membrane_topology.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_score_sum.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_sse_chirality.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_sse_completeness.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_sse_linear_loop_proximity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_sse_packing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_protein_model_sse_pairs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_radius_of_gyration.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_read_histograms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_residual_dipolar_coupling_histogram.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_residual_dipolar_coupling_q_value.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_atom_attraction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_atom_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_body_protein_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_distance_epr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_distance_spin_label.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_energy_well.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_noe_attraction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_noe_knowledge_based.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_residual_dipolar_coupling.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_saxs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_restraint_xlink.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_saxs_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_scores.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_slicelist_benchmark.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_membrane_alignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_pair_angle_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_pair_connectivity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_pair_contact.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_pair_packing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_pool_sses.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_sse_predictions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_score_strand_pairing.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
