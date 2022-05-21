CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclcommand)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_app_default_flags.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_command.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_command_line_writer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_command_state.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_default_flag_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_flag_dynamic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_flag_separator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_flag_static.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_flag_static_and_dynamic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_guesser.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_allowed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_allowed_non_const.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_extension.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_extensions_file_existence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_file_existence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_file_in_search_path.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_or.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_ranged.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_command_parameter_check_serializable.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
