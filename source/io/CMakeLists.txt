CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclio)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
	LIBRARY_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_binary_serialize.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_directory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_directory_entry.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_file.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_file_stream_buffer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_file_stream_buffer_encrypted.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_fixed_line_width_writer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_ifstream.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_ofstream.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serialization.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serialization_builtin.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serialization_interface.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serialization_with_check.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serialization_with_min_max.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serialize.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_serializer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_stream_buffer_classes.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_stream_buffer_interface.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_stream_interface.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/bcl_io_validation_result.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
	ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
