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

#ifndef BCL_ALIGN_H_
#define BCL_ALIGN_H_

// include the namespace forward header
#include "bcl_align.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically
#include <string>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_align.h
  //! @brief This namespace collects all alignment related functionality.
  //!
  //! @see @link example_align.cpp @endlink
  //! @author heinzes1
  //! @date Oct 30 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace align
  {
    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    //! @brief template function to extract a character id from a t_Member object;
    //! has to be specialized for every type used with align::Handler*
    //! @param MEMBER the t_Member object to extract the character id from
    //! @return the character id
    template< typename t_Member>
    char GetCharId( const t_Member &MEMBER);

    //! @brief template function to extract complete identifier from a t_Member object;
    //! has to be specialized for every type used with align::Handler*
    //! @param MEMBER the t_Member object to extract the complete identifier from
    //! @return the complete identifier
    template< typename t_Member>
    std::string GetCompleteId( const t_Member &MEMBER);

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_H_
