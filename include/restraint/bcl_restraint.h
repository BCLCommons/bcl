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

#ifndef BCL_RESTRAINT_H_
#define BCL_RESTRAINT_H_

// include the namespace forward header
#include "bcl_restraint.fwd.hh"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_restraint.h
  //! @brief TODO: document
  //! @details TODO: document
  //!
  //! @see @link example_restraint.cpp @endlink
  //! @author alexanns
  //! @date Jul 22, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace restraint
  {

    //! @brief return command line flag for specifying desired restraint types to use
    //! @return command line flag for specifying desired restraint types to use
    BCL_API const util::ShPtr< command::FlagInterface> &GetFlagRestraintsTypes();

    //! @brief return command line flag for specifying the prefix prepended to the each restraint's postfix
    //! @return command line flag for specifying the prefix prepended to the each restraint's postfix
    BCL_API util::ShPtr< command::FlagInterface> &GetFlagRestraintsFilePrefix();

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_H_
