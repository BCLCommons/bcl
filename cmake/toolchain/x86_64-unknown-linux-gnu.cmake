# this one is important
SET( CMAKE_SYSTEM_NAME Linux)
SET( CMAKE_SYSTEM_PROCESSOR x86_64)
# this one not so much
SET( CMAKE_SYSTEM_VERSION 1)

# specify the compiler
INCLUDE( MacroWrapCompiler)

MACRO_WRAP_COMPILER( CMAKE_C_COMPILER "gcc")
MACRO_WRAP_COMPILER( CMAKE_CXX_COMPILER "c++")
MACRO_WRAP_COMPILER( CMAKE_CXX_COMPILER "g++")

