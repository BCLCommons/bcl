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

#ifndef BCL_MATH_CONST_FUNCTION_H_
#define BCL_MATH_CONST_FUNCTION_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface_serializable.h"
#include "util/bcl_util_class_descriptor.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstFunction
    //! @brief This class is derived from function interface and is designed to give back
    //!        a constant value initialized by instantiation of function, also called 'constant function'.
    //!
    //! @tparam t_DataType template parameter for data type of passed argument
    //!
    //! @see @link example_math_const_function.cpp @endlink
    //! @author butkiem1
    //! @date Dec 2, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class ConstFunction :
      public FunctionInterfaceSerializable< t_ArgumentType, t_ResultType>
    {

    private:

    //////////
    // data //
    //////////

      //! @brief const result
      t_ResultType m_Result;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConstFunction() :
        m_Result()
      {
      }

      //! @brief default constructor
      //! @param RESULT the const result that will be returned by the operator
      ConstFunction( const t_ResultType &RESULT) :
        m_Result( RESULT)
      {
      }

      //! @brief Clone function
      //! @return pointer to new ConstFunction< t_ArgumentType>
      ConstFunction< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new ConstFunction< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const
      {
        static const std::string s_scheme( "ConstFunction");
        return s_scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator taking an ARGUMENT and returning a const result
      //! @param ARGUMENT Argument to be used to evaluate the function
      //! @return function value of the given argument
      t_ResultType operator()( const t_ArgumentType &ARGUMENT) const
      {
        // return the const result
        return m_Result;
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
        // read member
        io::Serialize::Read( m_Result, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Result, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class ConstFunction

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> ConstFunction< t_ArgumentType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new ConstFunction< t_ArgumentType, t_ResultType>())
    );

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_CONST_FUNCTION_H_
