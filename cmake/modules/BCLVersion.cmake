# Update the version
# The BCL version number.
SET( BCL_SVN_REVISION 0)

MACRO( GET_RELEASE_REVISION_FROM_SVN REV_VARIABLE)
  SET( ${REV_VARIABLE} 0)
  FIND_PROGRAM( SVNVERSION_CMD svnversion HINTS ENV PATH)
  IF( SVNVERSION_CMD)
    MESSAGE( STATUS "Calling svnversion on source root for BCL_RELEASE. This will take a little while.")
    EXECUTE_PROCESS( COMMAND ${SVNVERSION_CMD} --no-newline WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE ${REV_VARIABLE})
    STRING( COMPARE EQUAL "${${REV_VARIABLE}}" "exported" BCL_SVN_REVISION_NOT_SET) # if not a working copy, svnversion returns 'exported'
    IF( BCL_SVN_REVISION_NOT_SET)
      SET( BCL_SVN_REVISION 0)
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO( GET_RELEASE_REVISION_FROM_GIT REV_VARIABLE)
  SET( ${REV_VARIABLE} 0)
  FIND_PROGRAM( GIT_CMD git) # check for git
  IF( GIT_CMD)
    # run git and check return code
    EXECUTE_PROCESS(
      COMMAND ${GIT_CMD} svn info
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE INFO_RESULT
      OUTPUT_VARIABLE INFO_STDOUT
      ERROR_VARIABLE INFO_STDERR # use error_variable to avoid error messages not being captured
    )
    IF( NOT INFO_RESULT)
      STRING( REGEX MATCH "Revision: [0-9]+" ${REV_VARIABLE} ${INFO_STDOUT})
      STRING( REGEX REPLACE "Revision: " "" ${REV_VARIABLE} ${${REV_VARIABLE}})

      # get last svn commit hash
      EXECUTE_PROCESS( 
        COMMAND ${GIT_CMD} log --grep=^git-svn-id: --first-parent -1 --oneline 
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        OUTPUT_VARIABLE LAST_SVN_COMMIT
      )
      STRING( SUBSTRING "${LAST_SVN_COMMIT}" 0 7 LAST_SVN_COMMIT)
      # extract commits since svn commit hash
      EXECUTE_PROCESS( COMMAND ${GIT_CMD} log ${LAST_SVN_COMMIT}..HEAD --oneline WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE LOCAL_COMMITS)
      STRING( LENGTH "${LOCAL_COMMITS}" HAS_LOCAL_COMMITS)
      IF( HAS_LOCAL_COMMITS)
        SET( ${REV_VARIABLE} ${${REV_VARIABLE}}C) # add "C" if local commits exist
      ENDIF()

      # extract local changes since last commit
      # it might be necessary to add "build/" to .git/info/exclude
      EXECUTE_PROCESS( COMMAND ${GIT_CMD} status --porcelain WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE LOCAL_CHANGES)
      STRING( LENGTH ${LOCAL_CHANGES} HAS_LOCAL_CHANGES)
      IF( HAS_LOCAL_CHANGES)
        SET( ${REV_VARIABLE} ${${REV_VARIABLE}}M) # add "M" if uncommited local modifications exist
      ENDIF()
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO( GET_REVISION_FROM_SVN REV_VARIABLE)
  SET( ${REV_VARIABLE} 0)
  FIND_PROGRAM( SVN_CMD svn)
  IF( SVN_CMD) # if svn executable exists
    # run svn command to determine revision
    EXECUTE_PROCESS( 
      COMMAND ${SVN_CMD} info
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE INFO_RESULT
      OUTPUT_VARIABLE INFO_STDOUT
      ERROR_VARIABLE INFO_STDERR # use error_variable to avoid error messages not being captured
    )

    # if svn command successful
    IF( NOT INFO_RESULT)
      STRING( REGEX MATCH "Revision: [0-9]+" ${REV_VARIABLE} ${INFO_STDOUT})
      STRING( REGEX REPLACE "Revision: " "" ${REV_VARIABLE} ${${REV_VARIABLE}})
    ENDIF()
  ENDIF()
ENDMACRO()

MACRO( GET_REVISION_FROM_GIT REV_VARIABLE)
  SET( ${REV_VARIABLE} 0)
  FIND_PROGRAM( GIT_CMD git)
  IF( GIT_CMD) # if git executable exists
    # run git command to determine revision
    EXECUTE_PROCESS(
      COMMAND ${GIT_CMD} svn info
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE INFO_RESULT
      OUTPUT_VARIABLE INFO_STDOUT
      ERROR_VARIABLE INFO_STDERR # use error_variable to avoid error messages not being captured
    )

    # if git command successful
    IF( NOT INFO_RESULT)
      STRING( REGEX MATCH "Revision: [0-9]+" ${REV_VARIABLE} ${INFO_STDOUT})
      STRING( REGEX REPLACE "Revision: " "" ${REV_VARIABLE} ${${REV_VARIABLE}})

      # append the abbreviated commit hash, so that the commit can found quickly
      EXECUTE_PROCESS(
        COMMAND ${GIT_CMD} log -1 --format=\"%h\"
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        RESULT_VARIABLE INFO_RESULT
        OUTPUT_VARIABLE INFO_STDOUT
        ERROR_VARIABLE INFO_STDERR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      STRING( REGEX REPLACE "\"" "" INFO_STDOUT ${INFO_STDOUT}) # remove all quotes
      SET( ${REV_VARIABLE} "${${REV_VARIABLE}}-git${INFO_STDOUT}")
    ENDIF()
  ENDIF()
ENDMACRO()

# Releases/licenses
IF( NOT BCL_LICENSE)
  SET( BCL_LICENSE false)
ELSE()
  SET( BCL_LICENSE true)
ENDIF()

IF( NOT BCL_RELEASE)
  SET( BCL_RELEASE false)
ELSE()
  SET( BCL_RELEASE true)
ENDIF()

# find the svn revision
# TODO: Write a script that actually parses svn status and determines whether local modifications are relevant to the 
# release.  As it is, necessary changes to line endings for windows, for example, would result in an incorrect release
# number when compiling that release.  For now, just trust the person that runs the installation script to not make
# important changes to the source code
#IF( BCL_RELEASE)
  # slow but adds character if there was a modification
#  GET_RELEASE_REVISION_FROM_SVN( BCL_SVN_REVISION)
#  IF( NOT BCL_SVN_REVISION)
#    GET_RELEASE_REVISION_FROM_GIT( BCL_SVN_REVISION)
#  ENDIF()
#ELSE()
  # fast but does not reflect local changes 
  GET_REVISION_FROM_SVN( BCL_SVN_REVISION)
  IF( NOT BCL_SVN_REVISION)
    GET_REVISION_FROM_GIT( BCL_SVN_REVISION)
  ENDIF()
#ENDIF()
MESSAGE( STATUS "current revision: ${BCL_SVN_REVISION}")

# Compute the full version string.
SET( BCL_VERSION ${BCL_VERSION_MAJOR}.${BCL_VERSION_MINOR}.${BCL_VERSION_PATCH})
MESSAGE( STATUS "BCL_VERSION ${BCL_VERSION}")
