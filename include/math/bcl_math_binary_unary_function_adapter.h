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

#ifndef BCL_MATH_BINARY_UNARY_FUNCTION_ADAPTER_H_
#define BCL_MATH_BINARY_UNARY_FUNCTION_ADAPTER_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_binary_function_interface_serializable.h"
#include "bcl_math_function_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryUnaryFunctionAdapter
    //! @brief This class is a FunctionAdapter class
    //! @details It is supposed to be used as an adapter class for two functions a,b->c c->d to be combined into
    //!
    //! @tparam t_ArgumentType1 Type of the first argument to the function
    //! @tparam t_ArgumentType2 Type of the second argument to the function
    //! @tparam t_IntermediateType Type of the inermediat type
    //! @tparam t_ResulType     Type of the result of the function
    //!
    //! @see @link example_math_binary_unary_function_adapter.cpp @endlink
    //! @author woetzen
    //! @date Apr 21, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_IntermediateType, typename t_ResultType>
    class BinaryUnaryFunctionAdapter :
      public BinaryFunctionInterfaceSerializable< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    {

    private:

    //////////
    // data //
    //////////

      //! first function a,b->c
      util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_IntermediateType> > m_FunctionABC;

      //! second Function c->d
      util::ShPtr< FunctionInterfaceSerializable< t_IntermediateType, t_ResultType> > m_FunctionCD;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from binary and unary function
      //! @param SP_FUNCTION_ABC binary function a,b->c
      //! @param SP_FUNCTION_CD unary function c->d
      BinaryUnaryFunctionAdapter
      (
        const util::ShPtr< BinaryFunctionInterface< t_ArgumentType1, t_ArgumentType2, t_IntermediateType> > &SP_FUNCTION_ABC,
        const util::ShPtr< FunctionInterfaceSerializable< t_IntermediateType, t_ResultType> > &SP_FUNCTION_CD
      ) :
        m_FunctionABC( SP_FUNCTION_ABC),
        m_FunctionCD( SP_FUNCTION_CD)
      {
      }

      //! @brief Clone function
      //! @return pointer to new BinaryUnaryFunctionAdapter
      BinaryUnaryFunctionAdapter *Clone() const
      {
        return new BinaryUnaryFunctionAdapter< t_ArgumentType1, t_ArgumentType2, t_IntermediateType, t_ResultType>( *this);
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

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator a,b->d
      //! @param ARGUMENT1 first argument
      //! @param ARGUMENT2 second argument
      //! @return result of function evaluation
      t_ResultType operator()( const t_ArgumentType1 &ARGUMENT1, const t_ArgumentType2 &ARGUMENT2) const
      {
        return m_FunctionCD->operator()( m_FunctionABC->operator()( ARGUMENT1, ARGUMENT2));
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
        io::Serialize::Read( m_FunctionABC, ISTREAM);
        io::Serialize::Read( m_FunctionCD , ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_FunctionABC, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_FunctionCD , OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // class BinaryUnaryFunctionAdapter

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_BINARY_UNARY_FUNCTION_ADAPTER_H_
