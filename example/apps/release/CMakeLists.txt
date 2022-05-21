CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_alignment.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_cluster.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_contact_prediction.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_create_sse_pool.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_fit_in_density.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_fit_in_density_minimize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_fold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_fusion_protein.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_jufo.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_optimize.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_pdb_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_pdb_convert.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_score_protein.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

SET( SUBDIRECTORY_LIBRARY_SOURCES ${LIBRARY_SOURCES} PARENT_SCOPE)
