// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_DEFINES_H_
#define BCL_DEFINES_H_

#if defined (__GNUC__)
#elif defined (_MSC_VER) //this block exist to avoid unnecessary warnings with Visual C++
  #pragma warning ( disable : 4146 ) // disable MDE warning: "unary minus operator applied to unsigned type, result still unsigned"
  #pragma warning ( disable : 4996 ) // disable MDR warning: "Function call with parameters that may be unsafe - this call relies on the caller to check that the passed values are correct"
  #pragma warning ( disable : 4251 ) // disable MDR warning: "needs to have dll-interface to be used by clients of class"
  #pragma warning ( disable : 4661 ) // disable MDR warning: "no suitable definition provided for explicit template instantiation request"
  // min and max macro defined for ms compiler - need to undefine to use the std::min and std::max
  #define NOMINMAX
  #define __PRETTY_FUNCTION__ __FUNCSIG__ //__PRETTY_FUNCTION__ is defined by gcc and is called __FUNCSIG__ in VS and has a slightly different look
#endif

// this define is needed to get the UINT32_C and UINT64_C macros, defined in stdint.h
// These macros add the correct suffix (U, UL, or ULL) to compile-time constants.
// These macros should always be used whenever defining a constant of a particular size
#define __STDC_CONSTANT_MACROS

#if !defined(__MINGW32__) && ( defined(__WIN32__) || defined(_WIN32))
  #define PATH_SEPARATOR '\\'
  #if defined(BCL_DLL_EXPORT)
    // all BCL_API tagged symbols are exported into the import library
    #define BCL_API
 //__declspec(dllexport)
    #define BCL_API_REAL __declspec(dllexport)
    // http://support.microsoft.com/kb/168958
    #define BCL_EXPIMP_TEMPLATE
  #elif defined(BCL_DLL_IMPORT)
    // all BCL_API tagged symbols need to be imported, when linking against the dll
    #define BCL_API __declspec(dllimport)
    #define BCL_API_REAL
    #define BCL_EXPIMP_TEMPLATE extern
    #pragma warning ( disable : 4231 ) // disable MDR warning: "nonstandard extension used : 'extern' before template explicit instantiation"
  #else
    // no dll nor linkage against the import library
    #define BCL_API
    #define BCL_API_REAL
    #define BCL_EXPIMP_TEMPLATE
  #endif
#else
  // not a windows system
  #define BCL_API
  #define BCL_API_REAL
  #define BCL_EXPIMP_TEMPLATE extern
  #define PATH_SEPARATOR '/'
#endif

// MSC prior to 2010 does not have stdint.h
// currently, we only use a typedef and a macro from this header, so put them in directly
// MSC_VER gives the compiler major and minor version concatenated, which is 1600 as of VS2010
#if defined( _MSC_VER) and _MSC_VER < 1600
  typedef unsigned __int64 uint64_t; // unsigned 64 bit integer type
  #define UINT64_C( x) x##ui64       // declares x to be an unsigned 64 bit constant
#else

  #ifndef __STDC_CONSTANT_MACROS
    // this define is needed to get the UINT64_C macro from stdint.h
    // This macros adds the correct suffix (U, UL, or ULL) to compile-time constants.
    // UINT64_C should be used whenever defining a constant that must be 64 bits
    #define __STDC_CONSTANT_MACROS
  #endif

  // use the stdint header directly
  #include <stdint.h>
#endif

// declare iostream before any BCL-includes to ensure that std::cout is available during static initialization
#include <iostream>

#endif // BCL_DEFINES_H_
