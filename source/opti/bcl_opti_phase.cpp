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
#include "opti/bcl_opti_phase.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

    //! @brief returns Phase as string
    //! @param TYPE the improvement type
    //! @return the string for the improvement type
    const std::string &GetPhaseName( const Phase &TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "Start",
        "Iteration",
        "End",
        "Always",
        GetStaticClassName< Phase>()
      };

      return s_descriptors[ TYPE];
    }

    //! @brief test whether one phase logically equals another
    //! @param A, B the phases to compare
    //! @return true if A == B
    bool PhasesEqual( const Phase &A, const Phase &B)
    {
      return A == B || A == e_Always || B == e_Always;
    }

  } // namespace opti
} // namespace bcl
