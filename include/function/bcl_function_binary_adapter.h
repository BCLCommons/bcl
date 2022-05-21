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

#ifndef BCL_FUNCTION_BINARY_ADAPTER_H_
#define BCL_FUNCTION_BINARY_ADAPTER_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_function_binary_interface.h"
#include "bcl_function_unary_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace function
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BinaryAdapter
    //! @brief This class is a BinaryAdapter class
    //! @details It is supposed to be used as an adapter class for two functions a1,a2->b b->c to be combined into
    //! one a1,a2->c
    //!
    //! @see @link example_function_binary_adapter.cpp @endlink
    //! @author woetzen
    //! @date 08.06.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType1, typename t_ArgumentType2, typename t_IntermediateTypeLeft, typename t_IntermediateTypeRight, typename t_ResultType>
    class BinaryAdapter :
      public BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      //! first function a->b
      util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_IntermediateTypeLeft> > m_FunctionA1A2B;

      //! second Function b->c
      util::ShPtr< UnaryInterface< t_IntermediateTypeRight, t_ResultType> > m_FunctionBC;

      //! scheme to be used
      std::string m_Scheme;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from Function a1,a2->b and b->c
      //! @param SP_FUNCTION_A1A2B ShPtr to Function a1,a2->b
      //! @param SP_FUNCTION_BC ShPtr to Function b->c
      BinaryAdapter
      (
        const util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_IntermediateTypeLeft> > &SP_FUNCTION_A1A2B,
        const util::ShPtr< UnaryInterface< t_IntermediateTypeRight, t_ResultType> > &SP_FUNCTION_BC
      ) :
        m_FunctionA1A2B( SP_FUNCTION_A1A2B),
        m_FunctionBC( SP_FUNCTION_BC),
        m_Scheme( GetClassIdentifier())
      {
      }

      //! @brief construct from Function a1,a2->b and b->c
      //! @param SP_FUNCTION_A1A2B ShPtr to Function a1,a2->b
      //! @param SP_FUNCTION_BC ShPtr to Function b->c
      //! @param SCHEME Scheme to be used
      BinaryAdapter
      (
        const util::ShPtr< BinaryInterface< t_ArgumentType1, t_ArgumentType2, t_IntermediateTypeLeft> > &SP_FUNCTION_A1A2B,
        const util::ShPtr< UnaryInterface< t_IntermediateTypeRight, t_ResultType> > &SP_FUNCTION_BC,
        const std::string &SCHEME
      ) :
        m_FunctionA1A2B( SP_FUNCTION_A1A2B),
        m_FunctionBC( SP_FUNCTION_BC),
        m_Scheme( SCHEME)
      {
      }

      //! virtual copy constructor
      BinaryAdapter< t_ArgumentType1, t_ArgumentType2, t_IntermediateTypeLeft, t_IntermediateTypeRight, t_ResultType> *Clone() const
      {
        return new BinaryAdapter< t_ArgumentType1, t_ArgumentType2, t_IntermediateTypeLeft, t_IntermediateTypeRight, t_ResultType>( *this);
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

      //! @brief return scheme
      //! @return scheme
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

    ///////////////
    // operators //
    ///////////////

      //! operator taking an ARGUMENT1 and ARGUMENT2 and returning a t_ResultType object
      t_ResultType operator()( t_ArgumentType1 &ARGUMENT1, t_ArgumentType2 &ARGUMENT2) const
      {
        return m_FunctionBC->operator()( m_FunctionA1A2B->operator()( ARGUMENT1, ARGUMENT2));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_FunctionA1A2B, ISTREAM);
        io::Serialize::Read( m_FunctionBC   , ISTREAM);
        io::Serialize::Read( m_Scheme       , ISTREAM);

        // end
        return ISTREAM;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_FunctionA1A2B, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_FunctionBC   , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Scheme       , OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class BinaryAdapter

  } // namespace function
} // namespace bcl

#endif //BCL_FUNCTION_BINARY_ADAPTER_H_
