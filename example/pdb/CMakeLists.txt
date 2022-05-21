CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplepdb)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_entry_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_entry_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_head.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_ligand.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_line.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_line_criterium.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_line_type_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_line_types.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_biomatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_body_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_loop_closure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_membrane.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_quality_docking.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_quality_membrane.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_quality_multimer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_printer_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_residue.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_site.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_pdb_tail.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
