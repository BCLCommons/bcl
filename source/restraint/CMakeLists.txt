CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclrestraint)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_accessibility_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_accessibility_aa_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_accessibility_profile.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_accessibility_profile_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_accessibility_change.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_atom_distance_heatmap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_atom_distance_mean_sd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_atom_distance_pymol.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_atom_distance_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_atom_distance_score_heatmap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_coordinate_distance_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_per_residue_rmsd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_per_residue_rmsd_between_ensembles.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_analyze_sas.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_atom_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_atom_distance_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_body.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_cone_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_contact_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_contains_body_origin.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_data_pairwise.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_data_set_pairwise.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_epr_accessibility_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_epr_decay.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_epr_decay_simulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_epr_distance_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_accessibility_aa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_atom_distance_assigned.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_body.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_data_set_pairwise_identifiers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_data_set_pairwise_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_epr_decay.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_handler_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_locator_coordinates_first_side_chain_atom.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_locator_coordinates_hydrogen.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_filter_aa_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_filter_coordinate_exclusion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_filter_euclidian_distance.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_filter_exposure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_filter_sse_size.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_filter_triangulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_data_set_pairwise_remove.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_transformation_matrix_3d_null.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_mutate_transformation_matrix_3d_rotate.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_noe_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_piesa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_pofr_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_pre_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_rdc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_rdc_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_rdc_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_analysis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_data_parameters.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_debye.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_density_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_distance_density_point.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_experimental_and_calculated_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_experimental_and_calculated_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_optimization.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_pofr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_scattering_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_scattering_point.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_sas_transformation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_saxs_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_saxs_data_reduction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_saxs_opti_result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_restraint_xlink_data.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
	ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
