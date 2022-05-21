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

#ifndef BCL_MATH_ASSIGNMENT_BY_COMPARISON_H_
#define BCL_MATH_ASSIGNMENT_BY_COMPARISON_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_assignment_operation_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AssignmentByComparison
    //! @brief this class adds instantiations for all comparison operations
    //!
    //! @see @link example_math_assignments.cpp @endlink
    //! @author butkiem1
    //! @date Feb 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, template< typename> class t_ComparisonType>
    class AssignmentByComparison :
      public AssignmentOperationInterface< t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! comparison object
      t_ComparisonType< t_ArgumentType> m_Comparison;

    public:

      //! single instance of that class
      static const util::SiPtr< const AssignmentOperationInterface< t_ArgumentType> > s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      AssignmentByComparison *Clone() const
      {
        return new AssignmentByComparison( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a verb for this operation type
      //! @return a verb for this operation type
      const std::string &GetVerb() const
      {
        static std::string s_name
        (
          m_Comparison( 0, 0)
          ? m_Comparison( 1, 0) // these have an == as a comparison component
            ? "GreaterEqual"
            : m_Comparison( 0, 1) ? "LessEqual" : "Equal"
          : // these do not have == as a component
            m_Comparison( 1, 0) // < || > || !=
            ? m_Comparison( 0, 1) ? "NotEqual" : "Greater"
            : "Less"
        );
        return s_name;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static std::string s_name
        (
          m_Comparison( 0, 0)
          ? m_Comparison( 1, 0) // these have an == as a comparison component
            ? ">="
            : m_Comparison( 0, 1) ? "<=" : "=="
          : // these do not have == as a component
            m_Comparison( 1, 0) // < || > || !=
            ? m_Comparison( 0, 1) ? "!=" : ">"
            : "<"
        );
        return s_name;
      }

      //! @brief test whether assignment operation yield the same results independent of argument order
      //! @return true if operation is symmetric - false otherwise
      bool IsSymmetric() const
      {
        return m_Comparison( 0, 1) == m_Comparison( 1, 0);
      }

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param LHS the argument on the left hand side of the equation (e.g. x = y, x is the LHS)
      //! @param RHS the argument on the left hand side of the equation (e.g. x = y, y is the RHS)
      virtual void operator()( t_ArgumentType &LHS, const t_ArgumentType &RHS) const
      {
        if( !util::IsNaN( LHS) && !util::IsNaN( RHS))
        {
          LHS = m_Comparison( LHS, RHS);
        }
        else
        {
          LHS = util::GetUndefined< t_ArgumentType>();
        }
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;
        parameters.SetClassDescription( "LHS " + GetAlias() + " RHS");
        return parameters;
      }

    }; // template class AssignmentByComparison

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API AssignmentByComparison< float, std::less>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AssignmentByComparison< float, std::less_equal>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AssignmentByComparison< float, std::greater>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AssignmentByComparison< float, std::greater_equal>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AssignmentByComparison< float, std::equal_to>;
    BCL_EXPIMP_TEMPLATE template class BCL_API AssignmentByComparison< float, std::not_equal_to>;

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGNMENT_BY_COMPARISON_H_
