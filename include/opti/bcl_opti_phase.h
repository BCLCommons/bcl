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

#ifndef BCL_OPTI_PHASE_H_
#define BCL_OPTI_PHASE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util_wrapper_enum.h"

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_opti_phase.h
    //! @brief enum wrapper indicating phase of optimization (at-start, end, or during optimization)
    //!
    //! @see @link example_opti_phase.cpp @endlink
    //! @author mendenjl
    //! @date Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //! enumerator for change types that indicate an improvement
    enum Phase
    {
      e_Start,           //!< Before any optimization has taken place
      e_Iteration,       //!< During iterations of the optimization
      e_End,             //!< At the conclusion of optimization
      e_Always,          //!< meta-phase meaning regardless of optimization phase
      s_NumberPhases
    };

    //! @brief returns Phase as string
    //! @param TYPE the improvement type
    //! @return the string for the improvement type
    const std::string &GetPhaseName( const Phase &TYPE);

    //! SSEInfoTypeEnum simplifies the usage of the SSEInfoType enum of this class
    typedef util::WrapperEnum< Phase, &GetPhaseName, s_NumberPhases> PhaseEnum;

    //! @brief test whether one phase logically equals another
    //! @param A, B the phases to compare
    //! @return true if A == B
    bool PhasesEqual( const Phase &A, const Phase &B);

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_PHASE_H_
