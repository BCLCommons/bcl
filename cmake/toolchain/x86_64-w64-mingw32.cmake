# this one is important
SET( CMAKE_SYSTEM_NAME Windows)
SET( CMAKE_SYSTEM_PROCESSOR x86_64)
# this one not so much
SET( CMAKE_SYSTEM_VERSION 1)
SET( TARGET_PLATFORM ${CMAKE_SYSTEM_PROCESSOR}-w64-mingw32)

# cmake looks for a compiler called "rc" if compiling for windows but generating UNIX makefiles
# This command can be used to convert .rc (text resource) files into .res (binary) or .coff (object/exe)
# Though the bcl does not currently use this executable, omitting the RC compiler will cause cmake to fail to generate
# a unix makefile.  Ditto for archiver
FIND_PROGRAM( CMAKE_RC_COMPILER ${TARGET_PLATFORM}-windres PATHS ENV PATH NO_DEFAULT_PATH)
FIND_PROGRAM( CMAKE_AR ${TARGET_PLATFORM}-ar PATHS ENV PATH NO_DEFAULT_PATH)

# determine mxe compiler location and write the compiler wrapper
INCLUDE( MacroWrapCompiler)
MACRO_WRAP_COMPILER( CMAKE_C_COMPILER ${TARGET_PLATFORM}-gcc)
MACRO_WRAP_COMPILER( CMAKE_CXX_COMPILER ${TARGET_PLATFORM}-g++)
