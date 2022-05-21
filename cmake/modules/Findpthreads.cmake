# Try to find pthreads
#  pthreads_FOUND        - system has pthreads
#  pthreads_INCLUDE_DIRS - the pthreads include directory
#  pthreads_LIBRARIES    - Link these to use pthreads

INCLUDE( LibFindMacros)

# this package
SET( pthreads_PACKAGE_NAME pthreads)
SET( pthreads_PACKAGE_VERSION ${pthreads_FIND_VERSION})

# always prefer pthreads if available
SET( CMAKE_THREAD_PREFER_PTHREAD TRUE)

# dependencies
SET( Threads_PACKAGE_NAME Threads)
LIBFIND_PACKAGE( ${pthreads_PACKAGE_NAME} ${Threads_PACKAGE_NAME})

SET( ${Threads_PACKAGE_NAME}_FOUND ${THREADS_FOUND})
IF( ${Threads_PACKAGE_NAME}_FOUND AND NOT CMAKE_USE_WIN32_THREADS_INIT)
  SET( ${Threads_PACKAGE_NAME}_INCLUDE_DIRS "./")
  SET( ${Threads_PACKAGE_NAME}_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
  SET( ${pthreads_PACKAGE_NAME}_INCLUDE_DIR "./")
  SET( ${pthreads_PACKAGE_NAME}_LIBRARY ";")
ELSE()
  SET( ${Threads_PACKAGE_NAME}_LIBRARIES ";")
  SET( ${Threads_PACKAGE_NAME}_INCLUDE_DIRS "./")
  # include directory
  FIND_PATH( ${pthreads_PACKAGE_NAME}_INCLUDE_DIR NAMES pthread.h)
  
  # library
  FIND_LIBRARY( ${pthreads_PACKAGE_NAME}_LIBRARY NAMES pthreadGC2-w32 pthreadGC2-w64) 
ENDIF()

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET( ${pthreads_PACKAGE_NAME}_PROCESS_INCLUDES ${pthreads_PACKAGE_NAME}_INCLUDE_DIR ${Threads_PACKAGE_NAME}_INCLUDE_DIRS)
SET( ${pthreads_PACKAGE_NAME}_PROCESS_LIBS ${pthreads_PACKAGE_NAME}_LIBRARY ${Threads_PACKAGE_NAME}_LIBRARIES)
LIBFIND_PROCESS( ${pthreads_PACKAGE_NAME})

# runtime library
IF( MINGW)
  FIND_LIBRARY( ${pthreads_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES pthreadGC2-w32.dll pthreadGC2-w64.dll)
ELSEIF( APPLE)
  FIND_LIBRARY( ${pthreads_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES libpthread.dylib)
ENDIF()