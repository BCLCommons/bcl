CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexampleopti)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_evolution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_golden_section.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_nelder_mead.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_powell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_root_bisect.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_root_newton.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_root_regula_falsi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_approximator_root_secant.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_all.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_combine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_convergence_argument.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_convergence_result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_divergence_argument.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_elapsed_time.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_improvement_ratio.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_n_step.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_number_iterations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_phase.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_rejected.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_result_changed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_result_threshold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_skipped_steps.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_criterion_unimproved.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_ensemble_filter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_ensemble_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_evolution_operation_select.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_evolution_population.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_evolution_population_member.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_improvement_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_optimization_identity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_phase.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_pipeline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_printer_argument_to_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_printer_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_printer_with_criterion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_step_status.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_tracker.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_opti_tracker_with_history.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
