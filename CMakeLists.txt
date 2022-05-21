# minimum cmake versions
CMAKE_MINIMUM_REQUIRED( VERSION 2.8.8)

# option for license and release
OPTION( BCL_LICENSE "compile executables with license information, e.g. an expiration time added to bclcommons published binaries" OFF)
OPTION( BCL_RELEASE "compile executables for release to " OFF)
OPTION( BCL_EXTERNALS "explicitly set a list of externals to be used instead of the default list ( e.g. 'freeglut;ati')" OFF)
OPTION( USE_GOLD_LINKER "set this option to use the gold linker; gold sometimes produces invalid executables (particularly on debug builds), so use this option only for development" OFF)
OPTION( USE_ATI_OCL "Use ATI's libOpenCL, if available" OFF)
# another option to ignoring signal handling entirely for Java JNI-interfacing code is to preload jsig, but this introduces yet another shared library into the mix, at no real benefit
OPTION( NO_OS_SIGNAL_HANDLING "Defer signal handling; necessary when compiling shared or static libs that will interface with java via JNI, which crashes if signals are already being handled by the BCL" OFF)
SET( ADDITIONAL_COMPILER_OPTIONS CACHE STRING "Additional compiler options, example usage: cmake -DADDITIONAL_COMPILER_OPTIONS=\"-D__USE_XOPEN2K8\"")

# By this option, all objects from each namespace is linked into an archive ({namespace}.a) which is just the collection
# of all .o files. These archives can then be linked instead of all object files. This does work with gnu toolchains but
# not for visual studio since these archives (.lib) are libs and hide all previously public symbols.
OPTION( LINK_NAMESPACE_ARCHIVES "set this option to link namespace objects into archives" OFF)

#cache the cmake build type, if it was given over the command line
# this must be placed before the PROJECT() command (see http://permalink.gmane.org/gmane.comp.programming.tools.cmake.user/15952)
IF( DEFINED CMAKE_BUILD_TYPE)
	 SET( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
	 SET( CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

# module dir; must be specified before PROJECT( ...), since PROJECT( calls the toolchain files, and the module path
# may be specified in them
LIST( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/")

# project name
PROJECT( bcl-project C CXX)

# add any additional compiler options as definitions; adds nothing if an empty string or no option was set
ADD_DEFINITIONS( ${ADDITIONAL_COMPILER_OPTIONS})

IF( NO_OS_SIGNAL_HANDLING)
	SET( SIGNAL_HANDLING_SWITCH "-DBCL_NO_OS_SIGNAL_HANDLING" CACHE STRING "preprocessor flag that commands to bcl to be built without os-signal handling")
ELSE()
	SET( SIGNAL_HANDLING_SWITCH "" CACHE STRING "preprocessor flag that commands to bcl to be built without os-signal handling")
ENDIF()
ADD_DEFINITIONS( ${SIGNAL_HANDLING_SWITCH})

# update the forward headers
FIND_PROGRAM( PYTHON_CMD python2)
IF( PYTHON_CMD)
	EXECUTE_PROCESS(
		COMMAND ${PYTHON_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/code/CreateNamespaceForwardHeaders.py
		${CMAKE_CURRENT_SOURCE_DIR} -o -f
		OUTPUT_VARIABLE FWD_HDR_OUT
		RESULT_VARIABLE FWD_HDR_RES
	)
	IF( FWD_HDR_RES)
		MESSAGE( STATUS "Error building forward headers: " ${FWD_HDR_OUT})
	ELSEIF( ${FWD_HDR_OUT})
		MESSAGE( STATUS "Rebuilt forward headers: " ${FWD_HDR_OUT})
	ENDIF()
	UNSET( FWD_HDR_OUT)
	UNSET( FWD_HDR_RES)
ENDIF()

# cache variables defined by all toolchain files; this way if cmake reruns itself (commonly because a cmakelist was
# updated), the last value used in the cache will be used
SET( CMAKE_SYSTEM_NAME ${CMAKE_SYSTEM_NAME} CACHE INTERNAL "name of the system to build for")
SET( CMAKE_SYSTEM_PROCESSOR ${CMAKE_SYSTEM_PROCESSOR} CACHE INTERNAL "architecture that will be compiled for")
SET( CMAKE_SYSTEM_VERSION ${CMAKE_SYSTEM_VERSION} CACHE INTERNAL "version of the system that is being compiled for")

# options
# disable external libraries
# cmake build type

###########################
# compiler and linker flags
###########################

# flag for 32bit compile
IF( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "powerpc" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i686")
	ADD_DEFINITIONS( "-m32")
	SET( CMAKE_EXE_LINKER_FLAGS "-m32 ${CMAKE_EXE_LINKER_FLAGS}")
	SET( CMAKE_MODULE_LINKER_FLAGS "-m32 ${CMAKE_MODULE_LINKER_FLAGS}")
	SET( CMAKE_SHARED_LINKER_FLAGS "-m32 ${CMAKE_SHARED_LINKER_FLAGS}")
ENDIF()

# if bcl license option was given
IF( BCL_LICENSE)
	SET( BCL_RELEASE ON)
	ADD_DEFINITIONS( "-DNOCALLSTACK")
ENDIF()

# install prefix
SET( BCL_INSTALL_PREFIX "./")
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

# version
INCLUDE( "bcl_version.cmake")
INCLUDE( BCLVersion)

# gnu cxx flags
IF( CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# -O3 is not recommended, since the optimizations are supposedly too aggressive
	IF( CMAKE_CXX_FLAGS_RELEASE)
		STRING( REGEX REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
	ENDIF()
	IF( CMAKE_C_FLAGS_RELEASE)
		STRING( REGEX REPLACE "-O3" "-O2" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
	ENDIF()
	ADD_DEFINITIONS( "-Wno-deprecated")
	ADD_DEFINITIONS( "-Wall")

	# enable C++11 standard if the compiler supports it
	INCLUDE( CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG( "-std=c++11" COMPILER_SUPPORTS_CXX11)
	IF( COMPILER_SUPPORTS_CXX11)
		SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11")
		IF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND APPLE)
		  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
          SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -stdlib=libc++")
		ENDIF()
	ELSE()
		MESSAGE( FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support.")
	ENDIF()

	IF( CMAKE_CXX_COMPILER_ID MATCHES "Clang")
		# these warnings are ones that clang is either incorrect about or which are absolutely irrelevant and should be
		# ignored by all but the most pedantic of programmers
		ADD_DEFINITIONS( "-Wno-unused-private-field")
		ADD_DEFINITIONS( "-Wno-overloaded-virtual")
		ADD_DEFINITIONS( "-Wno-unneeded-internal-declaration")
		ADD_DEFINITIONS( "-Wno-unnamed-type-template-args")
		ADD_DEFINITIONS( "-Qunused-arguments")
	ENDIF()

	INCLUDE( CheckCXXCompilerFlag)
	IF( NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
		ADD_DEFINITIONS( "-fPIC")

		# this can be a flag, so that runtime libraries do not to be installed - but they should not be packaged into the installer then
		# and no other external needs to be linked, since they might still require our most recent libstdc++.so!!
		# ADD_DEFINITIONS( "-static-libgcc")
		# ADD_DEFINITIONS( "-static-libstdc++")
#    IF( BCL_LICENSE)
#      SET( CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
#      SET( CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
#    ENDIF()
	ELSE()
		Check_CXX_Compiler_Flag( "-Wno-virtual-move-assign" CXX_CAN_TURN_OFF_VIRTUAL_MOVE_ASSIGN)
		IF( CXX_CAN_TURN_OFF_VIRTUAL_MOVE_ASSIGN)
			ADD_DEFINITIONS( "-Wno-virtual-move-assign")
		ENDIF()
		Check_CXX_Compiler_Flag( "-fvisibility-inlines-hidden" CXX_HAS_VISIBILITY_INLINES_HIDDEN_FLAG)
        IF( CXX_HAS_VISIBILITY_INLINES_HIDDEN_FLAG)
            ADD_DEFINITIONS( "-fvisibility=hidden")
            ADD_DEFINITIONS( "-fvisibility-inlines-hidden")
        ENDIF()
	ADD_DEFINITIONS("-fno-inline-functions-called-once")
        #IF( CMAKE_CXX_FLAGS_RELEASE)
        #        STRING( REGEX REPLACE "-O2" "-O1" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
        #ENDIF()
        #IF( CMAKE_C_FLAGS_RELEASE)
        #        STRING( REGEX REPLACE "-O2" "-O1" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
        #ENDIF()
	ENDIF()

	# test if compiler would produce pretty template names - this needs to be disabled for now, so that the Read and Write work properly

	Check_CXX_Compiler_Flag( "-fno-pretty-templates" CXX_HAS_PRETTY_TEMPLATES_FLAG)
	IF( NOT ( CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
		IF( CXX_HAS_PRETTY_TEMPLATES_FLAG)
		 ADD_DEFINITIONS( "-fno-pretty-templates")
		ENDIF()
	ENDIF()
	IF( APPLE)
		SET( BCL_LINK_WHOLE_ARCHIVE "-Wl,-force_load")

		# -g will typically default to using STABS debug output format.
		# There is a bug (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=28837) in the current cross-compiler that causes
		#   messages of this form when used in debug mode on apple
		# ld warning: can't find atom for N_GSYM stabs (some static variable):G(0,some number) in some_source_file.so
		# Outputting dwarf symbols (-gdwarf-2) will fix the N_GSYM warnings (http://gcc.gnu.org/bugzilla/show_bug.cgi?id=34719)
		# Using dwarf symbols while linking to shared libraries that do not have dwarf symbols raises warnings while
		# debugging
		# Debug level 1 is necessary for building the entire bcl (the examples work fine with -g) because ld will
		# consume more than 4GB, and the current cross-compiler is only 32 bit.
		STRING( REGEX REPLACE "-g" "-gdwarf-2 -g1" CMAKE_CXX_FLAGS_DEBUG ${CMAKE_CXX_FLAGS_DEBUG})
		STRING( REGEX REPLACE "-g" "-gdwarf-2 -g1" CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
	ELSE()
		SET( BCL_LINK_WHOLE_ARCHIVE "-Wl,-whole-archive")
		SET( BCL_LINK_NO_WHOLE_ARCHIVE "-Wl,-no-whole-archive")
	ENDIF()

	IF( USE_GOLD_LINKER)
		# use gold linker on linux; gold provides 3-5x faster link times but only works for linux x86/x86_64
		IF( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86" OR ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
			IF( ${CMAKE_SYSTEM_NAME} STREQUAL "Linux")

				# add -B to force g++ to search the linker directory first when looking for ld/collect2, which are symlinks to gold
				SET( LINKER_FLAGS "-B/blue/meilerlab/apps/Linux2/x86_64/gold_linker/")

				IF( NOT ${CMAKE_BUILD_TYPE} STREQUAL "Release")
					# other flags necessary to keep gold from producing a binary that segfaults
					# icf,none prevents omission of various constructors that gold doesn't think will be used, but may be needed for
					#          the debugger
					# -E export all dynamic symbols
					SET( LINKER_FLAGS "${LINKER_FLAGS} -Wl,--icf,none -Wl,-E")
				ELSE()
					# do not add a callstack to release builds
					ADD_DEFINITIONS( "-DNOCALLSTACK")
				ENDIF()
				SET( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LINKER_FLAGS}" CACHE INTERNAL "linker flags to use gold")
				SET( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${LINKER_FLAGS}" CACHE INTERNAL "linker flags to use gold")
				SET( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${LINKER_FLAGS} " CACHE INTERNAL "linker flags to use gold")
			ENDIF()
		ENDIF()
	ENDIF()
ENDIF()

# use, i.e. don't skip the full RPATH for the build tree
#SET( CMAKE_SKIP_BUILD_RPATH TRUE)

# when building, don't use the install RPATH already
# (but later on when installing)
#SET( CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# the RPATH to be used when installing
SET( CMAKE_INSTALL_RPATH ".")

IF(APPLE)
	SET(CMAKE_INSTALL_NAME_DIR "@executable_path")
	SET( CMAKE_INSTALL_RPATH ".")
	SET( CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
ENDIF()

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# includes
INCLUDE_DIRECTORIES(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/apps
)

# library of useful macros
INCLUDE( MacroLibrary)

#######
# cpack
#######

INCLUDE( "bcl_package.cmake")
# also includes CPack which provides the cpack_add_component macro, that is used to define components and groups
# it also defines some install groups, that can than be used for easier packing

###############################
# external dependency detection
###############################

INCLUDE( "bcl_externals.cmake")

#############################################
# runtime libraries that need to be installed
#############################################

INCLUDE( "bcl_runtime_libraries.cmake")

#################################
# special PTHREAD case
# hopefully will be obsolete soon
#################################

IF( NOT pthreads_FOUND AND NOT CMAKE_HAVE_PTHREAD_H)
	ADD_DEFINITIONS( "-DBCL_WITHOUT_PTHREAD") # do not compile pthread dependent code in source/sched/bcl_sched_mutex.cpp
ENDIF()

#######################
## settings for eclipse
#######################

# it is necessary to set the message length to 0 so that errors are all put on one line
IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

################
# global targets
################

# statically linked
SET( BCL_TARGET_STATIC static)
MACRO_ADD_CUSTOM_TARGET(
	${BCL_TARGET_STATIC}
	COMMENT "applications statically linked to bcl"
)

# dynamically linked
SET( BCL_TARGET_SHARED shared)
MACRO_ADD_CUSTOM_TARGET(
	${BCL_TARGET_SHARED}
	COMMENT "applications dynamically linked to the bcl library"
)
SET_TARGET_PROPERTIES( ${BCL_TARGET_SHARED} PROPERTIES BUILD_WITH_INSTALL_RPATH true)

# apps
SET( BCL_TARGET_APPS apps)
MACRO_ADD_CUSTOM_TARGET(
	${BCL_TARGET_APPS}
	COMMENT "all applications"
)

################
# subdirectories
################

MACRO_ADD_SUBDIRECTORY( source)
MACRO_ADD_SUBDIRECTORY( apps)
MACRO_ADD_SUBDIRECTORY( example)

#################
# other resources
#################

# histograms
INSTALL(
	DIRECTORY histogram
	DESTINATION .
	COMPONENT BclReleaseAll
)

# (jufo) models
INSTALL(
	DIRECTORY model/jufo
	DESTINATION model
	COMPONENT BclReleaseAll
)

# opencl kernels
INSTALL(
	DIRECTORY opencl_kernels
	DESTINATION .
	COMPONENT BclReleaseAll
)

# rotamer library
INSTALL(
	DIRECTORY rotamer_library
	DESTINATION .
	COMPONENT BclReleaseAll
)

# license
INSTALL(
	FILES "installer/License.txt"
	DESTINATION .
	COMPONENT BclLicense
)

# ReadMe
INSTALL(
	FILES "installer/ReadMe_${CMAKE_SYSTEM_NAME}.txt"
	DESTINATION .
	COMPONENT BclReleaseAll
	RENAME "ReadMe.txt"
)

# including cpack itself finalizes the installation process for the release
INCLUDE( CPack)
