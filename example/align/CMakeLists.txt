CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclexamplealign)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_aligner_dp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_aligner_dynamic_programming.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_aligner_merge.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_aligner_progressive.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_aligner_shift.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_aligner_wordbased.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_alignment_hit.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_alignment_leaf.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_alignment_node.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_alignment_word.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_assignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_handler_blc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_handler_blocked.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_handler_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_handler_fasta.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_handler_pir.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_handler_standard.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_multiple_aligner_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_pairwise_aligner_classes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_sequence.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_word_generator_high_scoring.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_align_word_generator_subsequences.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS ${LIBRARY_OBJECTS} PARENT_SCOPE)
