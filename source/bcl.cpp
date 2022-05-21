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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "bcl.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl //! global namespace of the biochemistry library
{

  //! @brief identifier for the name space
  //! @return the name of the namespace
  const std::string &GetNamespaceIdentifier()
  {
    // please do leave that line, until the visual studio linker error is fixed. The ExtractNamespaceIdentifier function
    // is a BCL_API function and bcl.cpp is in the bcl.dll project and cannot import symbols, so you would get linker error
    static const std::string *s_namespace_name( new std::string( "bcl"));
//    static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
    return *s_namespace_name;
  }

  //! @brief get the copyright
  //! @return the copyright
  const std::string &GetCopyright()
  {
    static const std::string s_copyright( "BioChemistry Library (BCL), Copyright Meiler Lab 2006-2014, www.meilerlab.org");
    return s_copyright;
  }

} // namespace bcl
