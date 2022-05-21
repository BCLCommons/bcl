CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclappsinternalbiology)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_align_to_msa.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_analyze_loops.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_analyze_restraint_agreement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_benchmark_sse_pool.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_calculate_cone_model_parameters.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_distribution_plot_from_tables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_exposure_prediction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_extract_sites.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_fold_analysis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_fold_recognition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_aa_pair_statistics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_alignment_description.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_exposure_description.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_jufo_description.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_random_orientation_in_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_sse_pair_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_generate_wordbased_aligner_statistics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_loop_template.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_msa2pssm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_nmr_file_convert.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_optimize_data_set_pairwise.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_pdb_to_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_perturb_protein.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_protein_comparison_gallery.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_saxs_statistics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_scatterplot_from_tables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_simulate_accessibility_restraints.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_simulate_distance_restraints.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_simulate_rdc_restraints.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_sspred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_statistic_noe.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_statistic_spin_label.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_template_model.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
	ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS $<TARGET_OBJECTS:${LIBRARY_OBJECTS}> PARENT_SCOPE)
