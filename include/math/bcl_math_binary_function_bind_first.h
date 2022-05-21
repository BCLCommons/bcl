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

#ifndef BCL_MATH_BINARY_FUNCTION_BIND_FIRST_H_
#define BCL_MATH_BINARY_FUNCTION_BIND_FIRST_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryFunctionBindFirst
    //! @brief It binds a first argument to a binary function, so that it becomes just a function with one argument
    //!
    //! @see @link example_math_binary_function_bind_first.cpp @endlink
    //! @author woetzen
    //! @date 20.04.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    class BinaryFunctionBindFirst :
      public FunctionInterfaceSerializable< t_ArgumentType2, t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      //! the binary function the first argumetn is bound
      util::ShPtr< util::BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > m_BinaryFunction;

      //! the bound first argument to the binary function
      t_ArgumentType1 m_Argument1;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType>()
      {
      }

      //! construct from binary function and a bound 1st argument
      //! @param BINARY_FUNCTION binary function to which the first argument is bound to
      //! @param BOUND_ARGUMENT1 the bound argument
      BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType>
      (
        const util::BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> &BINARY_FUNCTION,
        const t_ArgumentType1 &BOUND_ARGUMENT1
      ) :
        m_BinaryFunction( BINARY_FUNCTION.Clone()),
        m_Argument1( BOUND_ARGUMENT1)
      {
      }

      //! construct from shptr to binary function and a bound 1st argument
      //! @param SP_BINARY_FUNCTION binary function to which the first argument is bound to
      //! @param BOUND_ARGUMENT1 the bound argument
      BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType>
      (
        const util::ShPtr< util::BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType> > &SP_BINARY_FUNCTION,
        const t_ArgumentType1 &BOUND_ARGUMENT1
      ) :
        m_BinaryFunction( SP_BINARY_FUNCTION),
        m_Argument1( BOUND_ARGUMENT1)
      {
      }

      //! copy constructor
      BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType> *Clone() const
      {
        return new BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType>( *this);
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

      //! @brief operator the implements the operation on one arguments returning a result
      //! @param ARGUMENT2 argument 2
      //! @return the Result of the operation
      t_ResultType operator()( const t_ArgumentType2 &ARGUMENT2) const
      {
        return m_BinaryFunction->operator()( m_Argument1, ARGUMENT2);
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
        // read members
        io::Serialize::Read( m_BinaryFunction, ISTREAM);
        io::Serialize::Read( m_Argument1, ISTREAM);

        // return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT indentation
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_BinaryFunction, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Argument1, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; //end template class BinaryFunctionBindFirst

    //! single instance of that class
    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new BinaryFunctionBindFirst< t_ArgumentType1, t_ArgumentType2, t_ResultType>())
    );

  } // namespace math
} // namespace bcl

#endif //BCL_MATH_BINARY_FUNCTION_BIND_FIRST_H_
