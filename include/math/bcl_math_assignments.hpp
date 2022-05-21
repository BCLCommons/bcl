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

#ifndef BCL_MATH_ASSIGNMENTS_HPP_
#define BCL_MATH_ASSIGNMENTS_HPP_

// include header of this class
#include "bcl_math_assignments.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_assignment_by_comparison.h"
#include "bcl_math_divide_equals.h"
#include "bcl_math_minus_equals.h"
#include "bcl_math_mod_equals.h"
#include "bcl_math_plus_equals.h"
#include "bcl_math_power_equals.h"
#include "bcl_math_times_equals.h"

// external includes - sorted alphabetically
#include <string>

namespace bcl
{
  namespace math
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief get all the different assignments
    //! @return vector of all different assignment types
    template< typename t_DataType>
    const storage::Vector< util::Implementation< AssignmentOperationInterface< t_DataType> > > &
      Assignments< t_DataType>::GetAssignments()
    {
      // construct all different assignments
      static const storage::Vector< util::Implementation< AssignmentOperationInterface< t_DataType> > >
        s_implementations
        (
          storage::Vector< util::Implementation< AssignmentOperationInterface< t_DataType> > >::Create
          (
            util::Implementation< AssignmentOperationInterface< t_DataType> >( PlusEquals< t_DataType>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( MinusEquals< t_DataType>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( DivideEquals< t_DataType>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( TimesEquals< t_DataType>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( PowerEquals< t_DataType>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( ModEquals< t_DataType>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( AssignmentByComparison< t_DataType, std::less>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( AssignmentByComparison< t_DataType, std::less_equal>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( AssignmentByComparison< t_DataType, std::greater>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( AssignmentByComparison< t_DataType, std::greater_equal>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( AssignmentByComparison< t_DataType, std::equal_to>()),
            util::Implementation< AssignmentOperationInterface< t_DataType> >( AssignmentByComparison< t_DataType, std::not_equal_to>())
          )
        );
      return s_implementations;
    }

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGNMENTS_HPP_
