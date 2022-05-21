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
#include "math/bcl_math_assignment_by_comparison.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    // register the class with the enumerated class instance
    template< typename t_ArgumentType, template< typename> class t_ComparisonType>
    const util::SiPtr< const AssignmentOperationInterface< t_ArgumentType> >
    AssignmentByComparison< t_ArgumentType, t_ComparisonType>::s_Instance
    (
      util::Enumerated< AssignmentOperationInterface< t_ArgumentType> >::AddInstance
      (
        new AssignmentByComparison< t_ArgumentType, t_ComparisonType>()
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API AssignmentByComparison< float, std::less>;
    template class BCL_API AssignmentByComparison< float, std::less_equal>;
    template class BCL_API AssignmentByComparison< float, std::greater>;
    template class BCL_API AssignmentByComparison< float, std::greater_equal>;
    template class BCL_API AssignmentByComparison< float, std::equal_to>;
    template class BCL_API AssignmentByComparison< float, std::not_equal_to>;

  } // namespace math
} // namespace bcl
