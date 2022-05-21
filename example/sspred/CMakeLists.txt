CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplesspred)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_b2tmpred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_boctopus.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_conpred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_dssp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_dssp_stride.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_jufo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_jufo9d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_kaksi.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_masp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_method_handler.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_methods.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_octopus.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_palsse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_partifold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_pdb.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_profphd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_proftmb.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_psipred.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_sam.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_sse_factory_highest.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_sse_factory_threshold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_stride.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_talos.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_tmbetanet.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_tmhmm.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_sspred_tmmod.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
