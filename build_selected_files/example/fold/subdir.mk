################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../example/fold/example_fold.cpp \
../example/fold/example_fold_add_parabolic_loops.cpp \
../example/fold/example_fold_collector_loop_domain_all_non_rigid.cpp \
../example/fold/example_fold_collector_loop_domain_random.cpp \
../example/fold/example_fold_collector_unconnected_sses.cpp \
../example/fold/example_fold_default_mutates.cpp \
../example/fold/example_fold_default_scores.cpp \
../example/fold/example_fold_handler_locator_loop_domain.cpp \
../example/fold/example_fold_handler_locator_loop_segment.cpp \
../example/fold/example_fold_locator_loop.cpp \
../example/fold/example_fold_locator_loop_domain.cpp \
../example/fold/example_fold_locator_loop_segment.cpp \
../example/fold/example_fold_locator_missing_coordinates.cpp \
../example/fold/example_fold_locator_unconnected_segments.cpp \
../example/fold/example_fold_loop_domain.cpp \
../example/fold/example_fold_loop_domain_c_to_n.cpp \
../example/fold/example_fold_loop_library.cpp \
../example/fold/example_fold_loop_parameters.cpp \
../example/fold/example_fold_loop_segment.cpp \
../example/fold/example_fold_loop_segment_sequence_order.cpp \
../example/fold/example_fold_mutate_aa_phi.cpp \
../example/fold/example_fold_mutate_aa_psi.cpp \
../example/fold/example_fold_mutate_aa_rotate.cpp \
../example/fold/example_fold_mutate_aa_sequence_grow.cpp \
../example/fold/example_fold_mutate_aa_set_phi.cpp \
../example/fold/example_fold_mutate_aa_set_psi.cpp \
../example/fold/example_fold_mutate_domain_flip.cpp \
../example/fold/example_fold_mutate_domain_merge_consecutive_ss_types.cpp \
../example/fold/example_fold_mutate_domain_shuffle.cpp \
../example/fold/example_fold_mutate_domain_sse_pair_trim.cpp \
../example/fold/example_fold_mutate_domain_sse_split.cpp \
../example/fold/example_fold_mutate_domain_transformation.cpp \
../example/fold/example_fold_mutate_loop_domain_dihedral.cpp \
../example/fold/example_fold_mutate_membrane.cpp \
../example/fold/example_fold_mutate_membrane_chain_move.cpp \
../example/fold/example_fold_mutate_multimer.cpp \
../example/fold/example_fold_mutate_protein_ensemble_add.cpp \
../example/fold/example_fold_mutate_protein_ensemble_remove.cpp \
../example/fold/example_fold_mutate_protein_model.cpp \
../example/fold/example_fold_mutate_protein_model_compress.cpp \
../example/fold/example_fold_mutate_protein_model_domain.cpp \
../example/fold/example_fold_mutate_protein_model_domain_add.cpp \
../example/fold/example_fold_mutate_protein_model_filter_conformations.cpp \
../example/fold/example_fold_mutate_protein_model_grow_sse.cpp \
../example/fold/example_fold_mutate_protein_model_loop_domain.cpp \
../example/fold/example_fold_mutate_protein_model_loop_domain_ccd.cpp \
../example/fold/example_fold_mutate_protein_model_loop_domain_grow.cpp \
../example/fold/example_fold_mutate_protein_model_loop_resize.cpp \
../example/fold/example_fold_mutate_protein_model_multiple_geometries.cpp \
../example/fold/example_fold_mutate_protein_model_pair_strands.cpp \
../example/fold/example_fold_mutate_protein_model_replicate_conformation.cpp \
../example/fold/example_fold_mutate_protein_model_sse.cpp \
../example/fold/example_fold_mutate_protein_model_sse_add.cpp \
../example/fold/example_fold_mutate_protein_model_sse_add_multiple.cpp \
../example/fold/example_fold_mutate_protein_model_sse_move.cpp \
../example/fold/example_fold_mutate_protein_model_sse_pair.cpp \
../example/fold/example_fold_mutate_protein_model_sse_pair_align_and_pull.cpp \
../example/fold/example_fold_mutate_protein_model_sse_pair_clash.cpp \
../example/fold/example_fold_mutate_protein_model_sse_pair_hinge.cpp \
../example/fold/example_fold_mutate_protein_model_sse_remove.cpp \
../example/fold/example_fold_mutate_protein_model_sse_resize.cpp \
../example/fold/example_fold_mutate_protein_model_sse_split.cpp \
../example/fold/example_fold_mutate_protein_model_sse_swap.cpp \
../example/fold/example_fold_mutate_protein_model_sse_swap_body.cpp \
../example/fold/example_fold_mutate_protein_model_sse_swap_multimer.cpp \
../example/fold/example_fold_mutate_protein_model_sse_swap_with_pool.cpp \
../example/fold/example_fold_mutate_protein_model_sse_swap_with_pool_overlap.cpp \
../example/fold/example_fold_mutate_protein_model_strand_switch_sheet.cpp \
../example/fold/example_fold_mutate_protein_model_switch_conformation.cpp \
../example/fold/example_fold_mutate_protein_model_thread_sequence.cpp \
../example/fold/example_fold_mutate_sheet_cycle.cpp \
../example/fold/example_fold_mutate_sheet_divide.cpp \
../example/fold/example_fold_mutate_sheet_fit_to_template.cpp \
../example/fold/example_fold_mutate_sheet_order.cpp \
../example/fold/example_fold_mutate_sheet_register_fix.cpp \
../example/fold/example_fold_mutate_sheet_register_shift.cpp \
../example/fold/example_fold_mutate_sheet_sort.cpp \
../example/fold/example_fold_mutate_sheet_twist.cpp \
../example/fold/example_fold_mutate_sse_bend_ramachandran.cpp \
../example/fold/example_fold_mutate_sse_bend_random.cpp \
../example/fold/example_fold_mutate_sse_bend_template.cpp \
../example/fold/example_fold_mutate_sse_type.cpp \
../example/fold/example_fold_mutate_tree.cpp \
../example/fold/example_fold_mutation_residue.cpp \
../example/fold/example_fold_phi_psi_generator_ccd.cpp \
../example/fold/example_fold_phi_psi_generator_ramachandran.cpp \
../example/fold/example_fold_placement_domain_using_fold_template.cpp \
../example/fold/example_fold_placement_sse_distance_restraint.cpp \
../example/fold/example_fold_placement_sse_into_body.cpp \
../example/fold/example_fold_placement_sse_next_to_sse.cpp \
../example/fold/example_fold_placement_sse_short_loop.cpp \
../example/fold/example_fold_placement_strand_next_to_sheet.cpp \
../example/fold/example_fold_protein_geometry.cpp 

OBJS += \
./example/fold/example_fold.o \
./example/fold/example_fold_add_parabolic_loops.o \
./example/fold/example_fold_collector_loop_domain_all_non_rigid.o \
./example/fold/example_fold_collector_loop_domain_random.o \
./example/fold/example_fold_collector_unconnected_sses.o \
./example/fold/example_fold_default_mutates.o \
./example/fold/example_fold_default_scores.o \
./example/fold/example_fold_handler_locator_loop_domain.o \
./example/fold/example_fold_handler_locator_loop_segment.o \
./example/fold/example_fold_locator_loop.o \
./example/fold/example_fold_locator_loop_domain.o \
./example/fold/example_fold_locator_loop_segment.o \
./example/fold/example_fold_locator_missing_coordinates.o \
./example/fold/example_fold_locator_unconnected_segments.o \
./example/fold/example_fold_loop_domain.o \
./example/fold/example_fold_loop_domain_c_to_n.o \
./example/fold/example_fold_loop_library.o \
./example/fold/example_fold_loop_parameters.o \
./example/fold/example_fold_loop_segment.o \
./example/fold/example_fold_loop_segment_sequence_order.o \
./example/fold/example_fold_mutate_aa_phi.o \
./example/fold/example_fold_mutate_aa_psi.o \
./example/fold/example_fold_mutate_aa_rotate.o \
./example/fold/example_fold_mutate_aa_sequence_grow.o \
./example/fold/example_fold_mutate_aa_set_phi.o \
./example/fold/example_fold_mutate_aa_set_psi.o \
./example/fold/example_fold_mutate_domain_flip.o \
./example/fold/example_fold_mutate_domain_merge_consecutive_ss_types.o \
./example/fold/example_fold_mutate_domain_shuffle.o \
./example/fold/example_fold_mutate_domain_sse_pair_trim.o \
./example/fold/example_fold_mutate_domain_sse_split.o \
./example/fold/example_fold_mutate_domain_transformation.o \
./example/fold/example_fold_mutate_loop_domain_dihedral.o \
./example/fold/example_fold_mutate_membrane.o \
./example/fold/example_fold_mutate_membrane_chain_move.o \
./example/fold/example_fold_mutate_multimer.o \
./example/fold/example_fold_mutate_protein_ensemble_add.o \
./example/fold/example_fold_mutate_protein_ensemble_remove.o \
./example/fold/example_fold_mutate_protein_model.o \
./example/fold/example_fold_mutate_protein_model_compress.o \
./example/fold/example_fold_mutate_protein_model_domain.o \
./example/fold/example_fold_mutate_protein_model_domain_add.o \
./example/fold/example_fold_mutate_protein_model_filter_conformations.o \
./example/fold/example_fold_mutate_protein_model_grow_sse.o \
./example/fold/example_fold_mutate_protein_model_loop_domain.o \
./example/fold/example_fold_mutate_protein_model_loop_domain_ccd.o \
./example/fold/example_fold_mutate_protein_model_loop_domain_grow.o \
./example/fold/example_fold_mutate_protein_model_loop_resize.o \
./example/fold/example_fold_mutate_protein_model_multiple_geometries.o \
./example/fold/example_fold_mutate_protein_model_pair_strands.o \
./example/fold/example_fold_mutate_protein_model_replicate_conformation.o \
./example/fold/example_fold_mutate_protein_model_sse.o \
./example/fold/example_fold_mutate_protein_model_sse_add.o \
./example/fold/example_fold_mutate_protein_model_sse_add_multiple.o \
./example/fold/example_fold_mutate_protein_model_sse_move.o \
./example/fold/example_fold_mutate_protein_model_sse_pair.o \
./example/fold/example_fold_mutate_protein_model_sse_pair_align_and_pull.o \
./example/fold/example_fold_mutate_protein_model_sse_pair_clash.o \
./example/fold/example_fold_mutate_protein_model_sse_pair_hinge.o \
./example/fold/example_fold_mutate_protein_model_sse_remove.o \
./example/fold/example_fold_mutate_protein_model_sse_resize.o \
./example/fold/example_fold_mutate_protein_model_sse_split.o \
./example/fold/example_fold_mutate_protein_model_sse_swap.o \
./example/fold/example_fold_mutate_protein_model_sse_swap_body.o \
./example/fold/example_fold_mutate_protein_model_sse_swap_multimer.o \
./example/fold/example_fold_mutate_protein_model_sse_swap_with_pool.o \
./example/fold/example_fold_mutate_protein_model_sse_swap_with_pool_overlap.o \
./example/fold/example_fold_mutate_protein_model_strand_switch_sheet.o \
./example/fold/example_fold_mutate_protein_model_switch_conformation.o \
./example/fold/example_fold_mutate_protein_model_thread_sequence.o \
./example/fold/example_fold_mutate_sheet_cycle.o \
./example/fold/example_fold_mutate_sheet_divide.o \
./example/fold/example_fold_mutate_sheet_fit_to_template.o \
./example/fold/example_fold_mutate_sheet_order.o \
./example/fold/example_fold_mutate_sheet_register_fix.o \
./example/fold/example_fold_mutate_sheet_register_shift.o \
./example/fold/example_fold_mutate_sheet_sort.o \
./example/fold/example_fold_mutate_sheet_twist.o \
./example/fold/example_fold_mutate_sse_bend_ramachandran.o \
./example/fold/example_fold_mutate_sse_bend_random.o \
./example/fold/example_fold_mutate_sse_bend_template.o \
./example/fold/example_fold_mutate_sse_type.o \
./example/fold/example_fold_mutate_tree.o \
./example/fold/example_fold_mutation_residue.o \
./example/fold/example_fold_phi_psi_generator_ccd.o \
./example/fold/example_fold_phi_psi_generator_ramachandran.o \
./example/fold/example_fold_placement_domain_using_fold_template.o \
./example/fold/example_fold_placement_sse_distance_restraint.o \
./example/fold/example_fold_placement_sse_into_body.o \
./example/fold/example_fold_placement_sse_next_to_sse.o \
./example/fold/example_fold_placement_sse_short_loop.o \
./example/fold/example_fold_placement_strand_next_to_sheet.o \
./example/fold/example_fold_protein_geometry.o 

CPP_DEPS += \
./example/fold/example_fold.d \
./example/fold/example_fold_add_parabolic_loops.d \
./example/fold/example_fold_collector_loop_domain_all_non_rigid.d \
./example/fold/example_fold_collector_loop_domain_random.d \
./example/fold/example_fold_collector_unconnected_sses.d \
./example/fold/example_fold_default_mutates.d \
./example/fold/example_fold_default_scores.d \
./example/fold/example_fold_handler_locator_loop_domain.d \
./example/fold/example_fold_handler_locator_loop_segment.d \
./example/fold/example_fold_locator_loop.d \
./example/fold/example_fold_locator_loop_domain.d \
./example/fold/example_fold_locator_loop_segment.d \
./example/fold/example_fold_locator_missing_coordinates.d \
./example/fold/example_fold_locator_unconnected_segments.d \
./example/fold/example_fold_loop_domain.d \
./example/fold/example_fold_loop_domain_c_to_n.d \
./example/fold/example_fold_loop_library.d \
./example/fold/example_fold_loop_parameters.d \
./example/fold/example_fold_loop_segment.d \
./example/fold/example_fold_loop_segment_sequence_order.d \
./example/fold/example_fold_mutate_aa_phi.d \
./example/fold/example_fold_mutate_aa_psi.d \
./example/fold/example_fold_mutate_aa_rotate.d \
./example/fold/example_fold_mutate_aa_sequence_grow.d \
./example/fold/example_fold_mutate_aa_set_phi.d \
./example/fold/example_fold_mutate_aa_set_psi.d \
./example/fold/example_fold_mutate_domain_flip.d \
./example/fold/example_fold_mutate_domain_merge_consecutive_ss_types.d \
./example/fold/example_fold_mutate_domain_shuffle.d \
./example/fold/example_fold_mutate_domain_sse_pair_trim.d \
./example/fold/example_fold_mutate_domain_sse_split.d \
./example/fold/example_fold_mutate_domain_transformation.d \
./example/fold/example_fold_mutate_loop_domain_dihedral.d \
./example/fold/example_fold_mutate_membrane.d \
./example/fold/example_fold_mutate_membrane_chain_move.d \
./example/fold/example_fold_mutate_multimer.d \
./example/fold/example_fold_mutate_protein_ensemble_add.d \
./example/fold/example_fold_mutate_protein_ensemble_remove.d \
./example/fold/example_fold_mutate_protein_model.d \
./example/fold/example_fold_mutate_protein_model_compress.d \
./example/fold/example_fold_mutate_protein_model_domain.d \
./example/fold/example_fold_mutate_protein_model_domain_add.d \
./example/fold/example_fold_mutate_protein_model_filter_conformations.d \
./example/fold/example_fold_mutate_protein_model_grow_sse.d \
./example/fold/example_fold_mutate_protein_model_loop_domain.d \
./example/fold/example_fold_mutate_protein_model_loop_domain_ccd.d \
./example/fold/example_fold_mutate_protein_model_loop_domain_grow.d \
./example/fold/example_fold_mutate_protein_model_loop_resize.d \
./example/fold/example_fold_mutate_protein_model_multiple_geometries.d \
./example/fold/example_fold_mutate_protein_model_pair_strands.d \
./example/fold/example_fold_mutate_protein_model_replicate_conformation.d \
./example/fold/example_fold_mutate_protein_model_sse.d \
./example/fold/example_fold_mutate_protein_model_sse_add.d \
./example/fold/example_fold_mutate_protein_model_sse_add_multiple.d \
./example/fold/example_fold_mutate_protein_model_sse_move.d \
./example/fold/example_fold_mutate_protein_model_sse_pair.d \
./example/fold/example_fold_mutate_protein_model_sse_pair_align_and_pull.d \
./example/fold/example_fold_mutate_protein_model_sse_pair_clash.d \
./example/fold/example_fold_mutate_protein_model_sse_pair_hinge.d \
./example/fold/example_fold_mutate_protein_model_sse_remove.d \
./example/fold/example_fold_mutate_protein_model_sse_resize.d \
./example/fold/example_fold_mutate_protein_model_sse_split.d \
./example/fold/example_fold_mutate_protein_model_sse_swap.d \
./example/fold/example_fold_mutate_protein_model_sse_swap_body.d \
./example/fold/example_fold_mutate_protein_model_sse_swap_multimer.d \
./example/fold/example_fold_mutate_protein_model_sse_swap_with_pool.d \
./example/fold/example_fold_mutate_protein_model_sse_swap_with_pool_overlap.d \
./example/fold/example_fold_mutate_protein_model_strand_switch_sheet.d \
./example/fold/example_fold_mutate_protein_model_switch_conformation.d \
./example/fold/example_fold_mutate_protein_model_thread_sequence.d \
./example/fold/example_fold_mutate_sheet_cycle.d \
./example/fold/example_fold_mutate_sheet_divide.d \
./example/fold/example_fold_mutate_sheet_fit_to_template.d \
./example/fold/example_fold_mutate_sheet_order.d \
./example/fold/example_fold_mutate_sheet_register_fix.d \
./example/fold/example_fold_mutate_sheet_register_shift.d \
./example/fold/example_fold_mutate_sheet_sort.d \
./example/fold/example_fold_mutate_sheet_twist.d \
./example/fold/example_fold_mutate_sse_bend_ramachandran.d \
./example/fold/example_fold_mutate_sse_bend_random.d \
./example/fold/example_fold_mutate_sse_bend_template.d \
./example/fold/example_fold_mutate_sse_type.d \
./example/fold/example_fold_mutate_tree.d \
./example/fold/example_fold_mutation_residue.d \
./example/fold/example_fold_phi_psi_generator_ccd.d \
./example/fold/example_fold_phi_psi_generator_ramachandran.d \
./example/fold/example_fold_placement_domain_using_fold_template.d \
./example/fold/example_fold_placement_sse_distance_restraint.d \
./example/fold/example_fold_placement_sse_into_body.d \
./example/fold/example_fold_placement_sse_next_to_sse.d \
./example/fold/example_fold_placement_sse_short_loop.d \
./example/fold/example_fold_placement_strand_next_to_sheet.d \
./example/fold/example_fold_protein_geometry.d 


# Each subdirectory must supply rules for building sources it contributes
example/fold/%.o: ../example/fold/%.cpp example/fold/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	distcc g++ -I../include -I../example -I../apps -I../extern/noarch/mysql/5.1.48/include -I../extern/noarch/mysqlpp/3.1.0/include -I../extern/noarch/bzip2/1.0.5/include -I../extern/noarch/zlib/1.2.5/include -I../extern/noarch/ati/2.5/include -O2 -Wall -c -fmessage-length=0 -Wno-deprecated -fno-pretty-templates -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


