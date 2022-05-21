CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplemodel)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_applicability_domain_kohonen.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_decision_tree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_kappa_nearest_neighbor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_kohonen_network.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_leverage_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_linear_regression.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_neural_network.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_neural_network_selective.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_restricted_boltzmann_machine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_support_vector_machine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_approximator_support_vector_machine_multi_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_collect_features_above.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_collect_features_top.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_cross_validation_info.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_multiplied.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_reduced_to_cluster_centers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_reduced_to_k_means.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_reduced_to_principal_components.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_select_columns.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_data_set_statistics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_decision_tree.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_descriptor_selection_backward_elimination.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_descriptor_selection_exhaustive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_descriptor_selection_feature_forward.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_dtree_binary_partition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_dtree_gini_index_data_partition_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_dtree_information_gain_data_partition_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_dtree_roc_data_partition_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_dtree_sequence_data_partition_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_feature_data_reference.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_feature_data_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_feature_reference.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_feature_result_and_state.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_interface_retrieve_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_interface_store_in_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_kappa_nearest_neighbor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_kohonen_network_applicability_domain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_kohonen_network_average.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_kohonen_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_leverage_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_meta_data_storage_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_multiple_linear_regression.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_perturb_attenuate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_perturb_max_norm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_selective_backpropagation_accuracy.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_selective_backpropagation_adaptive_tolerance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_selective_backpropagation_balanced.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_selective_backpropagation_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_selective_backpropagation_hybrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_selective_backpropagation_tolerance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_update_weights_bounded_simple_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_update_weights_resilient_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_neural_network_update_weights_simple_propagation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_accuracy.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_accuracy_with_excluded_range.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_auc_roc_curve.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_binary_operation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_categorical_max.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_constant.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_contingency_matrix_measure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_enrichment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_enrichment_average.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_information_gain_ratio.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_integral_precision_fraction_predicted.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_integral_tnr_tpr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_mae.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_partial.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_segment_overlap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_objective_function_wrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_pretrain_neural_network_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_pretrain_stacked_auto_encoder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_rescale_feature_data_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_restricted_boltzmann_machine_layer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_balanced.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_bootstrap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_by_feature.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_by_id.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_by_result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_chunk.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_combined.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_encoded_by_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_from_delimited_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_from_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_randomized.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_rescaled.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_rows.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_data_set_yscramble.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_retrieve_dataset_subset.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_f_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_input_sensitivity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_input_sensitivity_discrete.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_neural_network_input_sensitivity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_non_redundant.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_partition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_score_dataset_pearson_correlation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_support_vector_kernel_polynomial.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_support_vector_kernel_rbf.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_support_vector_machine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_support_vector_machine_multi_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_train_restricted_boltzmann_machine_layer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_transfer_gaussian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_transfer_linear.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_transfer_rectifier.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_model_transfer_sigmoid.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
