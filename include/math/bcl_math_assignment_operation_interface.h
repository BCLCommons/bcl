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

#ifndef BCL_MATH_ASSIGNMENT_OPERATION_INTERFACE_H_
#define BCL_MATH_ASSIGNMENT_OPERATION_INTERFACE_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AssignmentOperationInterface
    //! @details Interface for binary function objects that perform simple assignment operations ( e.g. =, +=, -=, *=)
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Apr 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType>
    class AssignmentOperationInterface :
      public util::SerializableInterface
    {

    public:

      typedef t_ArgumentType       first_argument_type;  //!< type of the first argument
      typedef const t_ArgumentType second_argument_type; //!< type of the first argument
      typedef void                 result_type;          //!< type of the return type

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      virtual AssignmentOperationInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief test whether assignment operation yield the same results independent of argument order
      //! @return by default false
      virtual bool IsSymmetric() const
      {
        return false;
      }

      //! @brief test whether assignment operation with zero LHS implies that the outcome will also be zero
      //!        in this case the right hand side need not be evaluated
      //! @return by default false
      virtual bool HasZeroProperty() const
      {
        return false;
      }

      //! @brief get a verb for this operation type
      //! @return a verb for this operation type
      virtual const std::string &GetVerb() const = 0;

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param LHS the argument on the left hand side of the equation (e.g. x = y, x is the LHS)
      //! @param RHS the argument on the left hand side of the equation (e.g. x = y, y is the RHS)
      virtual void operator()( t_ArgumentType &LHS, const t_ArgumentType &RHS) const = 0;

    }; // class AssignmentOperationInterface

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGNMENT_OPERATION_INTERFACE_H_
