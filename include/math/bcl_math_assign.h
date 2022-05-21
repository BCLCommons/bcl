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

#ifndef BCL_MATH_ASSIGN_H_
#define BCL_MATH_ASSIGN_H_

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
    //! @class Assign
    //! @brief binary function object that performs operator= on t_ArgumentTypes
    //!
    //! @see @link example_math_assignments.cpp @endlink
    //! @author mendenjl
    //! @date Apr 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_ArgumentType>
    class Assign :
      public AssignmentOperationInterface< t_ArgumentType>
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      Assign *Clone() const
      {
        return new Assign( *this);
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
        static const std::string s_name( "Assign");
        return s_name;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "=");
        return s_name;
      }

      //! @brief operator the implements the assignment operation on the two arguments returning a result
      //! @param LHS the argument on the left hand side of the equation (e.g. x = y, x is the LHS)
      //! @param RHS the argument on the left hand side of the equation (e.g. x = y, y is the RHS)
      virtual void operator()( t_ArgumentType &LHS, const t_ArgumentType &RHS) const
      {
        LHS = RHS;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer parameters;
        parameters.SetClassDescription( "LHS " + GetAlias() + " RHS");
        return parameters;
      }

    }; // template class Assign

    // register the class with the enumerated class instance
    template< typename t_ArgumentType>
    const util::SiPtr< const util::ObjectInterface> Assign< t_ArgumentType>::s_Instance
    (
      util::Enumerated< AssignmentOperationInterface< t_ArgumentType> >::AddInstance( new Assign< t_ArgumentType>())
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Assign< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Assign< float>;

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGN_H_
