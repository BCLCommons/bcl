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

#ifndef BCL_MATH_ASSIGNMENT_UNARY_STANDARD_H_
#define BCL_MATH_ASSIGNMENT_UNARY_STANDARD_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_assignment_unary_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AssignmentUnaryStandard
    //! @brief performs any of the standard template library functions
    //!
    //! @see @link example_math_assignment_unary_standard.cpp @endlink
    //! @author mendenjl
    //! @date Mar 11, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API AssignmentUnaryStandard :
      public AssignmentUnaryInterface
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! actual function to call
      double ( *m_Function)( double);

      //! Alias of the function
      std::string m_Alias;

      //! @brief helper function, get alias and description for a function
      //! @param FUNCTION the function of interest
      static std::pair< std::string, std::string> GetUnaryAssignmentAliasDescription( double ( *FUNCTION)( double));

      //! @brief function to add all the instances
      static util::SiPtr< const util::ObjectInterface> AddInstances();

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! Enum for additional functions that this class can use
      enum FunctionType
      {
        e_Not,
        e_Square,
        e_Negative
      };

      //! @brief constructor from function
      //! @param FUNCTION the function to use
      AssignmentUnaryStandard( double ( *FUNCTION)( double));

      //! @brief constructor from local enum
      AssignmentUnaryStandard( FunctionType FUNC);

      //! virtual copy constructor
      AssignmentUnaryStandard *Clone() const
      {
        return new AssignmentUnaryStandard( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        return m_Alias;
      }

      //! @brief operator the implements the operation on the argument
      //! @param VALUE the value to operate on
      void operator()( float &VALUE) const
      {
        VALUE = ( *m_Function)( VALUE);
      }

      //! @brief operator the implements the operation on the argument
      //! @param VALUE the value to operate on
      void operator()( double &VALUE) const
      {
        VALUE = ( *m_Function)( VALUE);
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AssignmentUnaryStandard

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_ASSIGNMENT_UNARY_STANDARD_H_
