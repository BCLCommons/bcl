CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplemath)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_angle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_assignment_unary_standard.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_assignments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_bicubic_spline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_binary_function_bind_first.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_binary_function_bind_second.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_binary_function_cached.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_combination.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_comparisons.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_const_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_contingency_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_contingency_matrix_measures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_cubic_spline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_cubic_spline_damped.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_cubic_spline_variable_delta.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_discrete_set_selector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_function_cached.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_gaussian_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_gnuplot_heatmap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_gnuplot_multiplot.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_histogram.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_histogram_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_histogram_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_identity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_kernel_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_limits.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_linear_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_linear_least_squares.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_log_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_combine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_decision_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_move_wrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_perturbation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_result.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_terminate_dependent.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_transformation_matrix_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_mutate_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_object_probability_distribution.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_object_stochastic_selector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_piecewise_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_polynomial.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_quadratic_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_quaternion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_range.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_range_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_roc_curve.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_rotation_matrix_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_rotation_matrix_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_running_average.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_running_average_sd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_running_min_max.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_running_sum.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_smooth_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_statistics.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_sum_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_sum_function_mixin.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_tensor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_transformation_matrix_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_tricubic_spline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_trigonometric_transition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_math_z_score.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
