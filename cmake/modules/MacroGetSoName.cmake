# get the so name of a given shared library with extension so
MACRO( MACRO_GET_SO_NAME LIB_FILE_NAME)
  IF( UNIX)
    FIND_PROGRAM( OBJDUMP_CMD objdump)
    IF( OBJDUMP_CMD)
      EXECUTE_PROCESS( COMMAND ${OBJDUMP_CMD} ${LIB_FILE_NAME} -p OUTPUT_VARIABLE OBJDUMP_OUT RESULT_VARIABLE SO_NAME_RESULT)
      IF( NOT SO_NAME_RESULT)
        STRING( REGEX MATCH "SONAME[^\n]+" SO_FILE_NAME ${OBJDUMP_OUT})
        STRING( REGEX REPLACE "SONAME[ ]+" "" SO_FILE_NAME ${SO_FILE_NAME})
      ENDIF()
      UNSET( OBJDUMP_OUT)
      UNSET( SO_NAME_RESULT)
    ENDIF()
  ENDIF()
ENDMACRO()
