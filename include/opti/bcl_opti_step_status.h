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

#ifndef BCL_OPTI_STEP_STATUS_H_
#define BCL_OPTI_STEP_STATUS_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_opti_step_status.h
  //! @brief the step status of an mc::Approximator
  //!
  //! @see @link example_opti_step_status.cpp @endlink
  //! @author mendenjl
  //! @date Sep 03, 2013
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace opti
  {

    //! enumerator for status of mc::Approximator steps, whether it is improved, accepted by metropolis or rejected
    enum StepStatus
    {
      e_Improved,
      e_Accepted,
      e_Rejected,
      e_Skipped,
      s_NumberStepStatus
    };

    //! @brief conversion to a string from a StepStatus
    //! @param STEP_STATUS the step status to get a string for
    //! @return a string representing that step status
    const std::string &GetStepStatusName( const StepStatus &STEP_STATUS);

    //! @brief enum class wrapper for Orientation
    typedef util::WrapperEnum< StepStatus, &GetStepStatusName, s_NumberStepStatus> StepStatusEnum;

  } // namespace opti
} // namespace bcl

#endif //BCL_OPTI_STEP_STATUS_H_
