#!/bin/sh

# bcl root is the starting path
BCL_ROOT=`pwd -P`

# get bcl version from bcl_version.cmake
BCL_VERSION=`grep MAJOR bcl_version.cmake | tr -d '[A-Z_() ]'`.`grep MINOR bcl_version.cmake | tr -d '[A-Z_() ]'`.`grep PATCH bcl_version.cmake | tr -d '[A-Z_() ]'`

# components that will be added to every installation 
DEFAULT_COMPONENTS="bcl_release_static;StdRuntimeLibraries;BclLicense;BclReadMe;BclMLModels;BclHistogram;bzip2;bzip2_LICENSE;zlib;zlib_LICENSE"
