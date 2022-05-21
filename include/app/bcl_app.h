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

#ifndef BCL_APP_H_
#define BCL_APP_H_

// include the namespace forward header
#include "bcl_app.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_app.h
  //! @namespace namespace with the application interface, enum, and all application classes
  //!
  //! @see @link example_app.cpp @endlink
  //! @author mendenjl
  //! @date Dec 10, 2012
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace app
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API const std::string &GetNamespaceIdentifier();

  } // namespace app
} // namespace bcl

#endif // BCL_APP_H_
