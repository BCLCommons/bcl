# 64-bit mac systems; intel core 2 duo or higher
# this one is important
SET( CMAKE_SYSTEM_NAME Darwin)
SET( CMAKE_SYSTEM_PROCESSOR x86_64)
# this one not so much
SET( CMAKE_SYSTEM_VERSION 1)

# specify the cross compiler (apple cross compiler handles both i686 and x86_64 architectures, so use that one again here)
INCLUDE( MacroWrapCompiler)
WRAP_COMPILER( CMAKE_C_COMPILER "i686-apple-darwin8-gcc")
WRAP_COMPILER( CMAKE_CXX_COMPILER "i686-apple-darwin8-g++")

# where is the target environment
#SET( CMAKE_FIND_ROOT_PATH /blue/meilerlab/apps/Linux2/x86_64/apple-darwin/2011.03.09/)
SET(
	APPLE_SDK_PATH
	/blue/meilerlab/apps/Linux2/x86_64/apple-darwin/2011.03.09/SDKs/MacOSX10.4u.sdk
	CACHE
	INTERNAL
	"Location of apple SDKs"
)
SET(
	APPLE_FRAMEWORKS_PATH
	${APPLE_SDK_PATH}/System/Library/Frameworks
	${APPLE_SDK_PATH}/usr/include
	CACHE INTERNAL "Location of apple frameworks and related header files"
)

# require at least version 10.4; this probably nets us around 98% of running macs; and certainly most of the 2% we miss
# have systems that wouldn't run the bcl for other reasons.  Apple stopped updating 10.3 in 2005.
# see, for example, this poll from 2009: http://hints.macworld.com/polls/index.php?pid=snowingyet
# also, http://update.omnigroup.com/ (click major version)
# Note: mysql may actually require mac os 10.5 or above
#SET_PROPERTY( GLOBAL PROPERTY MAC_OS_X_VERSION_MIN_REQUIRED 1040)

# search for programs in the build host directories
#SET( CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)

# for libraries and headers in the target directories
#SET( CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
#SET( CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
