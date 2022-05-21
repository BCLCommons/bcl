CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

SET(
  LIBRARY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_align_to_scaffold.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_build_scaffold_library.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_compare.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_coordinates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_filter.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_properties.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_reorder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_split.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/example_app_molecule_unique.cpp
)
MACRO_CHECK_CMAKE_DIRECTORY_LIST( ${CMAKE_CURRENT_SOURCE_DIR} "${LIBRARY_SOURCES}")

SET( SUBDIRECTORY_LIBRARY_SOURCES ${LIBRARY_SOURCES} PARENT_SCOPE)
