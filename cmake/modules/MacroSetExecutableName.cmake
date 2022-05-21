MACRO( MACRO_SET_EXECUTABLE_NAME)
  SET( EXECUTABLE_BASENAME "bcl")
  # gnucxx does not add .exe automatically
  IF( UNIX)
    SET( CMAKE_EXECUTABLE_SUFFIX_CXX ".exe")
  ELSE( UNIX)
  ENDIF( )
ENDMACRO()
