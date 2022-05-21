CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclutil)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_assert.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_call_stack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_class_descriptor.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_class_name_standardizer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_cleanable_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_color_gradient.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_colors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_cpp_data_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_data_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_enumerated.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_enums_instances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_format.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_functional_type.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_implementation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_implementation_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_logger_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_logger_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_loggers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_memory_usage.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_message.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_object_data_label.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_object_data_label_tokenizer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_object_instances.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_object_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_ptr_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_runtime_environment_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_runtime_environment_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_runtime_environments.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_serializable_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_sh_ptr.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_static_initialization_fiasco_finder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_stopwatch.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_string_functions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_string_numeric_conversion.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_string_replacement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_template_instantiations.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_time.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_util_undefined.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
	ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
