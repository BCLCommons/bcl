CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_descriptor_analyze.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_descriptor_generate_dataset.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_descriptor_refine_by_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_descriptor_score_dataset.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_descriptor_sequential_feature_selection.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

SET( SUBDIRECTORY_LIBRARY_SOURCES ${LIBRARY_SOURCES} PARENT_SCOPE)
