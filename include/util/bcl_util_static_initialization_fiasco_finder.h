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

#ifndef BCL_UTIL_STATIC_INITIALIZATION_FIASCO_FINDER_H_
#define BCL_UTIL_STATIC_INITIALIZATION_FIASCO_FINDER_H_

//#define ENABLE_FIASCO_FINDER
#ifdef ENABLE_FIASCO_FINDER

// include forward header of this class

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
// do not include bcl_util.h here!! otherwise this mechanism does not work!!

// external includes - sorted alphabetically
#include <iostream>
#include <string>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_util_static_initialization_fiasco_finder.h
  //! @brief functions for debugging static initialization
  //!
  //! @remarks example unnecessary
  //! @author woetzen
  //! @date Jan 26, 2011
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace util
  {
    //! @brief write the current file name and the number of the file that is initialized
    //! gives a hint about the linking order and the static initialization order, to detect where a segfault might arise
    //! from if it happens before the entry into main()
    //! @author woetzen
    //! @date Jan 25, 2011
    //! @remarks example unnecessary
    //! @see http://www.parashift.com/c++-faq-lite/ctors.html#faq-10.12 @endsee
    //! To enable this feature, define ENABLE_FIASCO_FINDER, rebuild, and run
    //! #define ENABLE_FIASCO_FINDER
    //! on the to of every cpp, this is to be included and BCL_StaticInitializationFiascoFinder has to be added:
    //! #include util/bcl_util_static_initialization_fiasco_finder.h
    //! BCL_StaticInitializationFiascoFinder
    //! @return bool, true which can be used to initialize the BCL_StaticInitializationFiascoFinder

    bool WriteFiasco( const std::string &FILE_NAME);

  } // namespace util
} // namespace bcl

// If you get a name collision on the following line, your usage is likely incorrect - only include this header in cpps!
#define BCL_StaticInitializationFiascoFinder static const bool g_PsuedoUniqueName = bcl::util::WriteFiasco(__FILE__);

#else // ENABLE_FIASCO_FINDER
// do nothing
#define BCL_StaticInitializationFiascoFinder
#endif // ENABLE_FIASCO_FINDER

#endif // BCL_UTIL_STATIC_INITIALIZATION_FIASCO_FINDER_H_
