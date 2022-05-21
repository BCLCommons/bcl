MACRO( MACRO_CHECK_CMAKE_DIRECTORY_LIST DIRECTORY SOURCES_INCLUDED)
  # this macro is obsolete when cmake is used from within eclipse
  # because the builder automatically updates the cmake lists first before building them
  FOREACH( FILE_NAME ${SOURCES_INCLUDED})
    IF( NOT(EXISTS ${FILE_NAME}))
      GET_FILENAME_COMPONENT( FILE_NAME_WITHOUT_PATH ${FILE_NAME} NAME)
      MESSAGE( "Non-existing ${FILE_NAME_WITHOUT_PATH} should be removed from ${DIRECTORY}/CMakeLists.txt") 
    ENDIF()
  ENDFOREACH()
  
  FILE(
    GLOB
    EXPECTED_LIBRARY_SOURCES
    ${DIRECTORY}/*.cpp
  )

  FOREACH( FILE_NAME ${EXPECTED_LIBRARY_SOURCES})
    SET( FOUND 0)
    FOREACH( INCLUDED_FILE_NAME ${SOURCES_INCLUDED})
      IF( ${INCLUDED_FILE_NAME} STREQUAL ${FILE_NAME})
        SET( FOUND 1)
      ENDIF()
    ENDFOREACH()
    IF( FOUND EQUAL 0)
      GET_FILENAME_COMPONENT( FILE_NAME_WITHOUT_PATH ${FILE_NAME} NAME)
      MESSAGE( "Missing ${FILE_NAME_WITHOUT_PATH} should be included in ${DIRECTORY}/CMakeLists.txt") 
    ENDIF()
  ENDFOREACH()
ENDMACRO()
