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

#ifndef BCL_OPENBLAS_H_
#define BCL_OPENBLAS_H_

// include the namespace forward header
#include "bcl_openblas.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_openblas.h
  //! @brief This namespace is for classes that depend on openblas
  //! @details All classes that use openblas should be declared and used only in source/openblas/
  //! The classes can then be linked to an enum, e.g. io::StreamBufferClasses,
  //! that provides similar functionality in the event the user does not have this library
  //!
  //! @see @link example_openblas.cpp @endlink
  //! @author mendenjl
  //! @date Nov 10, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace openblas
  {
    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

  } // namespace openblas
} // namespace bcl

#endif // BCL_OPENBLAS_H_
