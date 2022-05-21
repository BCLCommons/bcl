################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/score/bcl_score.cpp \
../source/score/bcl_score_aa_assignment_blast_profile.cpp \
../source/score/bcl_score_aa_assignment_blosum.cpp \
../source/score/bcl_score_aa_assignment_identity.cpp \
../source/score/bcl_score_aa_assignment_mean_similarity_matrix.cpp \
../source/score/bcl_score_aa_assignment_pam.cpp \
../source/score/bcl_score_aa_assignment_phat.cpp \
../source/score/bcl_score_aa_assignment_property.cpp \
../source/score/bcl_score_aa_assignment_ss_prediction.cpp \
../source/score/bcl_score_aa_assignments.cpp \
../source/score/bcl_score_aa_neighborhood_distances.cpp \
../source/score/bcl_score_aa_neighborhood_exposure.cpp \
../source/score/bcl_score_aa_neighborhood_exposure_prediction.cpp \
../source/score/bcl_score_aa_pair_atom_clash.cpp \
../source/score/bcl_score_aa_pair_clash.cpp \
../source/score/bcl_score_aa_pair_contact.cpp \
../source/score/bcl_score_aa_pair_contact_energy.cpp \
../source/score/bcl_score_aa_pair_distance.cpp \
../source/score/bcl_score_aa_pair_distance_fitted_function.cpp \
../source/score/bcl_score_aa_pair_distance_smooth.cpp \
../source/score/bcl_score_aa_pair_hi_res_clash.cpp \
../source/score/bcl_score_aa_pair_sidechain_interaction.cpp \
../source/score/bcl_score_aa_sequence.cpp \
../source/score/bcl_score_aa_sequence_pair.cpp \
../source/score/bcl_score_accessibility.cpp \
../source/score/bcl_score_accessibility_hydrophobic_moment.cpp \
../source/score/bcl_score_accessibility_hydrophobic_moment_magnitude.cpp \
../source/score/bcl_score_alignment_quality.cpp \
../source/score/bcl_score_assignment_gap_simple.cpp \
../source/score/bcl_score_body_assignment.cpp \
../source/score/bcl_score_body_connectivity_density.cpp \
../source/score/bcl_score_body_extent_agreement.cpp \
../source/score/bcl_score_body_extent_position_agreement.cpp \
../source/score/bcl_score_consensus_enrichment.cpp \
../source/score/bcl_score_contact_order.cpp \
../source/score/bcl_score_data_set_pairwise_bipolar.cpp \
../source/score/bcl_score_data_set_pairwise_coordinate_exclusion.cpp \
../source/score/bcl_score_data_set_pairwise_coordinate_triangulation.cpp \
../source/score/bcl_score_data_set_pairwise_data_density.cpp \
../source/score/bcl_score_data_set_pairwise_distance_change_magnitude.cpp \
../source/score/bcl_score_data_set_pairwise_euclidian_distance.cpp \
../source/score/bcl_score_data_set_pairwise_residue_type_exclusion.cpp \
../source/score/bcl_score_data_set_pairwise_sequence_separation.cpp \
../source/score/bcl_score_data_set_pairwise_size.cpp \
../source/score/bcl_score_data_set_pairwise_sse_center.cpp \
../source/score/bcl_score_data_set_pairwise_sse_connection.cpp \
../source/score/bcl_score_data_set_pairwise_sse_size.cpp \
../source/score/bcl_score_data_set_pairwise_sse_term.cpp \
../source/score/bcl_score_data_set_pairwise_structural_exposure.cpp \
../source/score/bcl_score_density_profile_sse_agreement.cpp \
../source/score/bcl_score_energy_distribution.cpp \
../source/score/bcl_score_environment_predictions.cpp \
../source/score/bcl_score_epr_accessibility.cpp \
../source/score/bcl_score_fuzzy_logic_filter.cpp \
../source/score/bcl_score_log_normal_distribution.cpp \
../source/score/bcl_score_loop.cpp \
../source/score/bcl_score_loop_angle.cpp \
../source/score/bcl_score_loop_closure.cpp \
../source/score/bcl_score_phi_psi.cpp \
../source/score/bcl_score_phi_psi_with_sspred.cpp \
../source/score/bcl_score_porf.cpp \
../source/score/bcl_score_protein_atom_density.cpp \
../source/score/bcl_score_protein_model.cpp \
../source/score/bcl_score_protein_model_aa_neighborhood.cpp \
../source/score/bcl_score_protein_model_aa_neighborhood_docking.cpp \
../source/score/bcl_score_protein_model_completeness.cpp \
../source/score/bcl_score_protein_model_defined_loops.cpp \
../source/score/bcl_score_protein_model_fragment_topology.cpp \
../source/score/bcl_score_protein_model_gap.cpp \
../source/score/bcl_score_protein_model_inverted.cpp \
../source/score/bcl_score_protein_model_loop_domain_closure.cpp \
../source/score/bcl_score_protein_model_membrane_topology.cpp \
../source/score/bcl_score_protein_model_score_sum.cpp \
../source/score/bcl_score_protein_model_sse.cpp \
../source/score/bcl_score_protein_model_sse_chirality.cpp \
../source/score/bcl_score_protein_model_sse_completeness.cpp \
../source/score/bcl_score_protein_model_sse_linear_loop_proximity.cpp \
../source/score/bcl_score_protein_model_sse_neighbors.cpp \
../source/score/bcl_score_protein_model_sse_packing.cpp \
../source/score/bcl_score_protein_model_sse_pairs.cpp \
../source/score/bcl_score_protein_model_topology.cpp \
../source/score/bcl_score_protein_model_wrapper.cpp \
../source/score/bcl_score_radius_of_gyration.cpp \
../source/score/bcl_score_read_histograms.cpp \
../source/score/bcl_score_residual_dipolar_coupling_histogram.cpp \
../source/score/bcl_score_residual_dipolar_coupling_q_value.cpp \
../source/score/bcl_score_restraint_atom_attraction.cpp \
../source/score/bcl_score_restraint_atom_distance.cpp \
../source/score/bcl_score_restraint_body_protein_model.cpp \
../source/score/bcl_score_restraint_distance_epr.cpp \
../source/score/bcl_score_restraint_distance_spin_label.cpp \
../source/score/bcl_score_restraint_energy_well.cpp \
../source/score/bcl_score_restraint_nmr_distance_interface.cpp \
../source/score/bcl_score_restraint_noe_attraction.cpp \
../source/score/bcl_score_restraint_noe_knowledge_based.cpp \
../source/score/bcl_score_restraint_pofr.cpp \
../source/score/bcl_score_restraint_residual_dipolar_coupling.cpp \
../source/score/bcl_score_restraint_saxs.cpp \
../source/score/bcl_score_restraint_xlink.cpp \
../source/score/bcl_score_sas_type.cpp \
../source/score/bcl_score_scores.cpp \
../source/score/bcl_score_sse_membrane_alignment.cpp \
../source/score/bcl_score_sse_pair_angle_distance.cpp \
../source/score/bcl_score_sse_pair_clash.cpp \
../source/score/bcl_score_sse_pair_connectivity.cpp \
../source/score/bcl_score_sse_pair_contact.cpp \
../source/score/bcl_score_sse_pair_gap.cpp \
../source/score/bcl_score_sse_pair_packing.cpp \
../source/score/bcl_score_sse_pairs_fragments.cpp \
../source/score/bcl_score_sse_pool_sses.cpp \
../source/score/bcl_score_sse_predictions.cpp \
../source/score/bcl_score_strand_pairing.cpp \
../source/score/bcl_score_symmetry.cpp 

OBJS += \
./source/score/bcl_score.o \
./source/score/bcl_score_aa_assignment_blast_profile.o \
./source/score/bcl_score_aa_assignment_blosum.o \
./source/score/bcl_score_aa_assignment_identity.o \
./source/score/bcl_score_aa_assignment_mean_similarity_matrix.o \
./source/score/bcl_score_aa_assignment_pam.o \
./source/score/bcl_score_aa_assignment_phat.o \
./source/score/bcl_score_aa_assignment_property.o \
./source/score/bcl_score_aa_assignment_ss_prediction.o \
./source/score/bcl_score_aa_assignments.o \
./source/score/bcl_score_aa_neighborhood_distances.o \
./source/score/bcl_score_aa_neighborhood_exposure.o \
./source/score/bcl_score_aa_neighborhood_exposure_prediction.o \
./source/score/bcl_score_aa_pair_atom_clash.o \
./source/score/bcl_score_aa_pair_clash.o \
./source/score/bcl_score_aa_pair_contact.o \
./source/score/bcl_score_aa_pair_contact_energy.o \
./source/score/bcl_score_aa_pair_distance.o \
./source/score/bcl_score_aa_pair_distance_fitted_function.o \
./source/score/bcl_score_aa_pair_distance_smooth.o \
./source/score/bcl_score_aa_pair_hi_res_clash.o \
./source/score/bcl_score_aa_pair_sidechain_interaction.o \
./source/score/bcl_score_aa_sequence.o \
./source/score/bcl_score_aa_sequence_pair.o \
./source/score/bcl_score_accessibility.o \
./source/score/bcl_score_accessibility_hydrophobic_moment.o \
./source/score/bcl_score_accessibility_hydrophobic_moment_magnitude.o \
./source/score/bcl_score_alignment_quality.o \
./source/score/bcl_score_assignment_gap_simple.o \
./source/score/bcl_score_body_assignment.o \
./source/score/bcl_score_body_connectivity_density.o \
./source/score/bcl_score_body_extent_agreement.o \
./source/score/bcl_score_body_extent_position_agreement.o \
./source/score/bcl_score_consensus_enrichment.o \
./source/score/bcl_score_contact_order.o \
./source/score/bcl_score_data_set_pairwise_bipolar.o \
./source/score/bcl_score_data_set_pairwise_coordinate_exclusion.o \
./source/score/bcl_score_data_set_pairwise_coordinate_triangulation.o \
./source/score/bcl_score_data_set_pairwise_data_density.o \
./source/score/bcl_score_data_set_pairwise_distance_change_magnitude.o \
./source/score/bcl_score_data_set_pairwise_euclidian_distance.o \
./source/score/bcl_score_data_set_pairwise_residue_type_exclusion.o \
./source/score/bcl_score_data_set_pairwise_sequence_separation.o \
./source/score/bcl_score_data_set_pairwise_size.o \
./source/score/bcl_score_data_set_pairwise_sse_center.o \
./source/score/bcl_score_data_set_pairwise_sse_connection.o \
./source/score/bcl_score_data_set_pairwise_sse_size.o \
./source/score/bcl_score_data_set_pairwise_sse_term.o \
./source/score/bcl_score_data_set_pairwise_structural_exposure.o \
./source/score/bcl_score_density_profile_sse_agreement.o \
./source/score/bcl_score_energy_distribution.o \
./source/score/bcl_score_environment_predictions.o \
./source/score/bcl_score_epr_accessibility.o \
./source/score/bcl_score_fuzzy_logic_filter.o \
./source/score/bcl_score_log_normal_distribution.o \
./source/score/bcl_score_loop.o \
./source/score/bcl_score_loop_angle.o \
./source/score/bcl_score_loop_closure.o \
./source/score/bcl_score_phi_psi.o \
./source/score/bcl_score_phi_psi_with_sspred.o \
./source/score/bcl_score_porf.o \
./source/score/bcl_score_protein_atom_density.o \
./source/score/bcl_score_protein_model.o \
./source/score/bcl_score_protein_model_aa_neighborhood.o \
./source/score/bcl_score_protein_model_aa_neighborhood_docking.o \
./source/score/bcl_score_protein_model_completeness.o \
./source/score/bcl_score_protein_model_defined_loops.o \
./source/score/bcl_score_protein_model_fragment_topology.o \
./source/score/bcl_score_protein_model_gap.o \
./source/score/bcl_score_protein_model_inverted.o \
./source/score/bcl_score_protein_model_loop_domain_closure.o \
./source/score/bcl_score_protein_model_membrane_topology.o \
./source/score/bcl_score_protein_model_score_sum.o \
./source/score/bcl_score_protein_model_sse.o \
./source/score/bcl_score_protein_model_sse_chirality.o \
./source/score/bcl_score_protein_model_sse_completeness.o \
./source/score/bcl_score_protein_model_sse_linear_loop_proximity.o \
./source/score/bcl_score_protein_model_sse_neighbors.o \
./source/score/bcl_score_protein_model_sse_packing.o \
./source/score/bcl_score_protein_model_sse_pairs.o \
./source/score/bcl_score_protein_model_topology.o \
./source/score/bcl_score_protein_model_wrapper.o \
./source/score/bcl_score_radius_of_gyration.o \
./source/score/bcl_score_read_histograms.o \
./source/score/bcl_score_residual_dipolar_coupling_histogram.o \
./source/score/bcl_score_residual_dipolar_coupling_q_value.o \
./source/score/bcl_score_restraint_atom_attraction.o \
./source/score/bcl_score_restraint_atom_distance.o \
./source/score/bcl_score_restraint_body_protein_model.o \
./source/score/bcl_score_restraint_distance_epr.o \
./source/score/bcl_score_restraint_distance_spin_label.o \
./source/score/bcl_score_restraint_energy_well.o \
./source/score/bcl_score_restraint_nmr_distance_interface.o \
./source/score/bcl_score_restraint_noe_attraction.o \
./source/score/bcl_score_restraint_noe_knowledge_based.o \
./source/score/bcl_score_restraint_pofr.o \
./source/score/bcl_score_restraint_residual_dipolar_coupling.o \
./source/score/bcl_score_restraint_saxs.o \
./source/score/bcl_score_restraint_xlink.o \
./source/score/bcl_score_sas_type.o \
./source/score/bcl_score_scores.o \
./source/score/bcl_score_sse_membrane_alignment.o \
./source/score/bcl_score_sse_pair_angle_distance.o \
./source/score/bcl_score_sse_pair_clash.o \
./source/score/bcl_score_sse_pair_connectivity.o \
./source/score/bcl_score_sse_pair_contact.o \
./source/score/bcl_score_sse_pair_gap.o \
./source/score/bcl_score_sse_pair_packing.o \
./source/score/bcl_score_sse_pairs_fragments.o \
./source/score/bcl_score_sse_pool_sses.o \
./source/score/bcl_score_sse_predictions.o \
./source/score/bcl_score_strand_pairing.o \
./source/score/bcl_score_symmetry.o 

CPP_DEPS += \
./source/score/bcl_score.d \
./source/score/bcl_score_aa_assignment_blast_profile.d \
./source/score/bcl_score_aa_assignment_blosum.d \
./source/score/bcl_score_aa_assignment_identity.d \
./source/score/bcl_score_aa_assignment_mean_similarity_matrix.d \
./source/score/bcl_score_aa_assignment_pam.d \
./source/score/bcl_score_aa_assignment_phat.d \
./source/score/bcl_score_aa_assignment_property.d \
./source/score/bcl_score_aa_assignment_ss_prediction.d \
./source/score/bcl_score_aa_assignments.d \
./source/score/bcl_score_aa_neighborhood_distances.d \
./source/score/bcl_score_aa_neighborhood_exposure.d \
./source/score/bcl_score_aa_neighborhood_exposure_prediction.d \
./source/score/bcl_score_aa_pair_atom_clash.d \
./source/score/bcl_score_aa_pair_clash.d \
./source/score/bcl_score_aa_pair_contact.d \
./source/score/bcl_score_aa_pair_contact_energy.d \
./source/score/bcl_score_aa_pair_distance.d \
./source/score/bcl_score_aa_pair_distance_fitted_function.d \
./source/score/bcl_score_aa_pair_distance_smooth.d \
./source/score/bcl_score_aa_pair_hi_res_clash.d \
./source/score/bcl_score_aa_pair_sidechain_interaction.d \
./source/score/bcl_score_aa_sequence.d \
./source/score/bcl_score_aa_sequence_pair.d \
./source/score/bcl_score_accessibility.d \
./source/score/bcl_score_accessibility_hydrophobic_moment.d \
./source/score/bcl_score_accessibility_hydrophobic_moment_magnitude.d \
./source/score/bcl_score_alignment_quality.d \
./source/score/bcl_score_assignment_gap_simple.d \
./source/score/bcl_score_body_assignment.d \
./source/score/bcl_score_body_connectivity_density.d \
./source/score/bcl_score_body_extent_agreement.d \
./source/score/bcl_score_body_extent_position_agreement.d \
./source/score/bcl_score_consensus_enrichment.d \
./source/score/bcl_score_contact_order.d \
./source/score/bcl_score_data_set_pairwise_bipolar.d \
./source/score/bcl_score_data_set_pairwise_coordinate_exclusion.d \
./source/score/bcl_score_data_set_pairwise_coordinate_triangulation.d \
./source/score/bcl_score_data_set_pairwise_data_density.d \
./source/score/bcl_score_data_set_pairwise_distance_change_magnitude.d \
./source/score/bcl_score_data_set_pairwise_euclidian_distance.d \
./source/score/bcl_score_data_set_pairwise_residue_type_exclusion.d \
./source/score/bcl_score_data_set_pairwise_sequence_separation.d \
./source/score/bcl_score_data_set_pairwise_size.d \
./source/score/bcl_score_data_set_pairwise_sse_center.d \
./source/score/bcl_score_data_set_pairwise_sse_connection.d \
./source/score/bcl_score_data_set_pairwise_sse_size.d \
./source/score/bcl_score_data_set_pairwise_sse_term.d \
./source/score/bcl_score_data_set_pairwise_structural_exposure.d \
./source/score/bcl_score_density_profile_sse_agreement.d \
./source/score/bcl_score_energy_distribution.d \
./source/score/bcl_score_environment_predictions.d \
./source/score/bcl_score_epr_accessibility.d \
./source/score/bcl_score_fuzzy_logic_filter.d \
./source/score/bcl_score_log_normal_distribution.d \
./source/score/bcl_score_loop.d \
./source/score/bcl_score_loop_angle.d \
./source/score/bcl_score_loop_closure.d \
./source/score/bcl_score_phi_psi.d \
./source/score/bcl_score_phi_psi_with_sspred.d \
./source/score/bcl_score_porf.d \
./source/score/bcl_score_protein_atom_density.d \
./source/score/bcl_score_protein_model.d \
./source/score/bcl_score_protein_model_aa_neighborhood.d \
./source/score/bcl_score_protein_model_aa_neighborhood_docking.d \
./source/score/bcl_score_protein_model_completeness.d \
./source/score/bcl_score_protein_model_defined_loops.d \
./source/score/bcl_score_protein_model_fragment_topology.d \
./source/score/bcl_score_protein_model_gap.d \
./source/score/bcl_score_protein_model_inverted.d \
./source/score/bcl_score_protein_model_loop_domain_closure.d \
./source/score/bcl_score_protein_model_membrane_topology.d \
./source/score/bcl_score_protein_model_score_sum.d \
./source/score/bcl_score_protein_model_sse.d \
./source/score/bcl_score_protein_model_sse_chirality.d \
./source/score/bcl_score_protein_model_sse_completeness.d \
./source/score/bcl_score_protein_model_sse_linear_loop_proximity.d \
./source/score/bcl_score_protein_model_sse_neighbors.d \
./source/score/bcl_score_protein_model_sse_packing.d \
./source/score/bcl_score_protein_model_sse_pairs.d \
./source/score/bcl_score_protein_model_topology.d \
./source/score/bcl_score_protein_model_wrapper.d \
./source/score/bcl_score_radius_of_gyration.d \
./source/score/bcl_score_read_histograms.d \
./source/score/bcl_score_residual_dipolar_coupling_histogram.d \
./source/score/bcl_score_residual_dipolar_coupling_q_value.d \
./source/score/bcl_score_restraint_atom_attraction.d \
./source/score/bcl_score_restraint_atom_distance.d \
./source/score/bcl_score_restraint_body_protein_model.d \
./source/score/bcl_score_restraint_distance_epr.d \
./source/score/bcl_score_restraint_distance_spin_label.d \
./source/score/bcl_score_restraint_energy_well.d \
./source/score/bcl_score_restraint_nmr_distance_interface.d \
./source/score/bcl_score_restraint_noe_attraction.d \
./source/score/bcl_score_restraint_noe_knowledge_based.d \
./source/score/bcl_score_restraint_pofr.d \
./source/score/bcl_score_restraint_residual_dipolar_coupling.d \
./source/score/bcl_score_restraint_saxs.d \
./source/score/bcl_score_restraint_xlink.d \
./source/score/bcl_score_sas_type.d \
./source/score/bcl_score_scores.d \
./source/score/bcl_score_sse_membrane_alignment.d \
./source/score/bcl_score_sse_pair_angle_distance.d \
./source/score/bcl_score_sse_pair_clash.d \
./source/score/bcl_score_sse_pair_connectivity.d \
./source/score/bcl_score_sse_pair_contact.d \
./source/score/bcl_score_sse_pair_gap.d \
./source/score/bcl_score_sse_pair_packing.d \
./source/score/bcl_score_sse_pairs_fragments.d \
./source/score/bcl_score_sse_pool_sses.d \
./source/score/bcl_score_sse_predictions.d \
./source/score/bcl_score_strand_pairing.d \
./source/score/bcl_score_symmetry.d 


# Each subdirectory must supply rules for building sources it contributes
source/score/%.o: ../source/score/%.cpp source/score/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


