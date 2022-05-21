CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexampleutil)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_assert.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_call_stack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_class_descriptor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_color_gradient.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_colors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_cpp_data_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_cpu_benchmark_whetstone.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_enum.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_format.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_implementation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_logger_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_loggers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_memory_usage.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_message.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_object_data_label.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_object_data_label_tokenizer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_object_instances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_object_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_runtime_environment_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_runtime_environment_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_runtime_environments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_sh_ptr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_sh_ptr_list.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_sh_ptr_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_si_ptr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_si_ptr_list.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_si_ptr_vector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_stopwatch.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_string_functions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_string_replacement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_thunk_wrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_time.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_undefined.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_wrapper.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_util_wrapper_base.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
	ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
