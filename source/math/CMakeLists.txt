CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclmath)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_angle.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_assign.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_assignment_by_comparison.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_assignment_unary_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_assignment_unary_standard.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_assignments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_bicubic_spline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_comparisons.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_contingency_matrix_measures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_cubic_spline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_cubic_spline_damped.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_cubic_spline_variable_delta.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_discrete_set_selector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_divide_equals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_gaussian_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_gnuplot.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_gnuplot_heatmap.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_gnuplot_multiplot.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_histogram.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_histogram_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_histogram_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_kernel_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_linear_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_linear_least_squares.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_log_likelihood.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_minus_equals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_mod_equals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_mutate_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_piecewise_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_plus_equals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_polynomial.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_power_equals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_quadratic_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_quaternion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_range.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_range_set.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_roc_curve.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_rotation_matrix_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_rotation_matrix_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_template_instantiations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_tensor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_times_equals.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_transformation_matrix_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_tricubic_spline.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_trigonometric_transition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_math_z_score.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
