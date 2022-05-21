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

#ifndef BCL_FUNCTION_MEMBER_H_
#define BCL_FUNCTION_MEMBER_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_function_unary_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace function
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Member
    //! @brief this class wraps member function into a function object, so that class members can be used in algorithms
    //!
    //! @tparam t_ClassType type of the class
    //! @tparam t_ReturnType return type of the member function
    //!
    //! @see @link example_function_member.cpp @endlink
    //! @author woetzen
    //! @date Aug 7, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ClassType, typename t_ResultType>
    class Member :
      public UnaryInterface< t_ClassType, t_ResultType>
    {

    public:

      //! @brief typedef for member function pointer
      typedef t_ResultType ( t_ClassType::*MemberFunctionPtr)();

    private:

    //////////
    // data //
    //////////

      MemberFunctionPtr m_FunctionPtr;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from member function pointer
      Member( MemberFunctionPtr MEMBER_FUNCTION_PTR) :
        m_FunctionPtr( MEMBER_FUNCTION_PTR)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Member< t_ClassType>
      Member< t_ClassType, t_ResultType> *Clone() const
      {
        return new Member< t_ClassType, t_ResultType>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief virtual operator taking a CLASS and returning a t_ResultType object
      //! @param CLASS class to be used to evaluate the function
      //! @return function value of the given class
      t_ResultType operator()( t_ClassType &CLASS) const
      {
        return ( CLASS.*m_FunctionPtr)();
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // template class Member

    //! @brief helper function to construct a Member class from a ptr to member function
    //! this can be done without the explicit use of template parameters
    //! @return Member object with the function pointer
    template< typename t_ClassType, typename t_ResultType>
    Member< t_ClassType, t_ResultType>
    MemberFunction( t_ResultType ( t_ClassType::*MEMBER_FUNCTION_PTR)())
    {
      return Member< t_ClassType, t_ResultType>( MEMBER_FUNCTION_PTR);
    }

  } // namespace function
} // namespace bcl

#endif // BCL_FUNCTION_MEMBER_H_ 
