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

#ifndef BCL_SCHED_H_
#define BCL_SCHED_H_

// include the namespace forward header
#include "bcl_sched.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_sched.h
  //! @brief sched for scheduling processes on different cpus using different programming interfaces
  //! @details those programming models can include but are not limited to
  //! \li pthread
  //! \li mpi
  //! \li cpu (serial)
  //! \li cuda
  //! \li opencl
  //!
  //! @see @link example_sched.cpp @endlink
  //! @author woetzen
  //! @date Nov 10, 2009
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace sched
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

  } // namespace sched
} // namespace bcl

#endif // BCL_SCHED_H_
