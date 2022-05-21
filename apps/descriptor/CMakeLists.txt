CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET( LIBRARY_NAME bclappsdescriptor)
SET( LIBRARY_OBJECTS ${LIBRARY_NAME}-objects)
PROJECT( lib${LIBRARY_NAME})

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_analyze.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_convert_code_object_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_dataset_similarity_measures.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_generate_dataset.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_generate_pca_eigenvectors.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_refine_by_score.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_score_dataset.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/bcl_app_descriptor_sequential_feature_selection.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

ADD_LIBRARY( ${LIBRARY_OBJECTS} OBJECT ${LIBRARY_SOURCES})

IF( LINK_NAMESPACE_ARCHIVES)
  ADD_LIBRARY( ${LIBRARY_NAME} STATIC $<TARGET_OBJECTS:${LIBRARY_OBJECTS}>)
ENDIF()

SET( SUBDIRECTORY_LIBRARY_NAME ${LIBRARY_NAME} PARENT_SCOPE)
SET( SUBDIRECTORY_LIBRARY_OBJECTS $<TARGET_OBJECTS:${LIBRARY_OBJECTS}> PARENT_SCOPE)
