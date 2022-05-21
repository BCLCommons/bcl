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
#include "math/bcl_math_assignment_unary_interface.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_assignment_unary_standard.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief get all the different unary assignments
    //! @return vector of all different unary assignment types
    const storage::Vector< util::Implementation< AssignmentUnaryInterface> >
      &AssignmentUnaryInterface::GetUnaryAssignments()
    {
      static const storage::Vector< util::Implementation< AssignmentUnaryInterface> > s_impls
      (
        storage::Vector< util::Implementation< AssignmentUnaryInterface> >::Create
        (
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( cos)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( sin)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( log)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( log10)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( exp)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( std::abs)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( std::sqrt)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( AssignmentUnaryStandard::e_Square)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( AssignmentUnaryStandard::e_Not)),
          util::Implementation< AssignmentUnaryInterface>( new AssignmentUnaryStandard( AssignmentUnaryStandard::e_Negative))
        )
      );
      return s_impls;
    }

  } // namespace math
} // namespace bcl
