CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclmc)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/../../include/")

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_movie_printer_chimera.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_movie_printer_interface.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_movie_printer_pymol.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_movie_printers.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutate_loop_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutate_loop_add_resize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutate_loop_fragment_add.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutate_loop_fragment_replace.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutate_loop_remove.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutate_loop_replace.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_mutates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_optimization_ccd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_optimization_docking.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_optimization_mcm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_stage.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_temperature_accepted.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_temperature_default.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_temperature_exponential.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_temperature_linear.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_mc_template_instantiations.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
