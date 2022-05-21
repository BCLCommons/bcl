# Write a compiler wrapper using distcc, if it is available and set the variable referenced by
# WRAPPER_CMAKE_VARIABLE_NAME to COMPILER_NAME
MACRO( MACRO_WRAP_COMPILER WRAPPER_CMAKE_VARIABLE_NAME COMPILER_NAME)
  # if distcc is available in $PATH, automatically use it
  # Note that if you have distcc installed but don't want to use it, it is sufficient to comment out the next line
  # however, its probably a better idea to move distcc out of the path, since the check_pump_make script will try to
  # run with all available distcc cores anyway
  FIND_PROGRAM( DISTCC_CMD distcc PATHS ENV PATH NO_DEFAULT_PATH)
  SET( DISTCC_CMD ${DISTCC_CMD} CACHE FILEPATH "distcc command for distributed compiling")
  FIND_PROGRAM( ${COMPILER_NAME}_PATH ${COMPILER_NAME} PATHS ENV PATH NO_DEFAULT_PATH)
  # dereference any symlinks in the path to either file
  GET_FILENAME_COMPONENT( ${COMPILER_NAME}_PATH ${${COMPILER_NAME}_PATH} REALPATH)
  IF( NOT DISTCC_CMD)
    # no distcc
    # Use the compiler in $PATH directly, which is marginally faster/better than calling a shell script all the time
    SET( ${WRAPPER_CMAKE_VARIABLE_NAME} ${${COMPILER_NAME}_PATH})
  ELSE()
    GET_FILENAME_COMPONENT( DISTCC_CMD ${DISTCC_CMD} REALPATH)
    # configure compiler wrapper
    # use distcc; write a compiler wrapper for it
    SET( DISTCC_COMPILER_WRAPPER_DIR ${CMAKE_BINARY_DIR}/distcc_compiler_wrapper)
    # all compiler wrappers go into a folder within the build tree
    IF( NOT EXISTS ${DISTCC_COMPILER_WRAPPER_DIR})
      MESSAGE( STATUS "creating directory for distcc compiler wrappers: " ${DISTCC_COMPILER_WRAPPER_DIR})
      FILE( MAKE_DIRECTORY ${DISTCC_COMPILER_WRAPPER_DIR})
    ENDIF()
    FILE( WRITE ${DISTCC_COMPILER_WRAPPER_DIR}/${COMPILER_NAME} "#!/bin/sh\nexport DISTCC_DIR=/tmp/$USER\n[ -d $DISTCC_DIR ] || mkdir -p $DISTCC_DIR\n" ${DISTCC_CMD} " " ${${COMPILER_NAME}_PATH} " $@")
    FIND_PROGRAM( CHMOD chmod)
    EXECUTE_PROCESS( COMMAND ${CHMOD} +x ${DISTCC_COMPILER_WRAPPER_DIR}/${COMPILER_NAME})
    # set the desired cmake variable
    SET( ${WRAPPER_CMAKE_VARIABLE_NAME} ${DISTCC_COMPILER_WRAPPER_DIR}/${COMPILER_NAME})
    
    # Note: as of 2.8.11.2, cmake has a bug that prevents it from understanding any compiler
    # with a space in the name, e.g. SET( CMAKE_CXX_COMPILER "distcc g++") will fail, claiming that the compiler was not
    # found.  Once this is fixed, it would speed up builds a little, and help ensure that the system doesn't run out of
    # process ids (leading to system failed vfork errors) if this configuration file is removed and we just do
    # SET( ${WRAPPER_CMAKE_VARIABLE_NAME} "${DISTCC_CMD} ${${COMPILER_NAME}_PATH}")
  ENDIF()
ENDMACRO()
