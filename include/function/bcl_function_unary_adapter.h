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

#ifndef BCL_FUNCTION_UNARY_ADAPTER_H_
#define BCL_FUNCTION_UNARY_ADAPTER_H_

// include the namespace header
#include "bcl_function.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_function_unary_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace function
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UnaryAdapter
    //! @brief This class is a UnaryAdapter class
    //! @details It is supposed to be used as an adapter class for two functions a->b b->c to be combined into
    //! one a->c
    //!
    //! @see @link example_function_unary_adapter.cpp @endlink
    //! @author woetzen
    //! @date 08.06.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_IntermediateType, typename t_ResultType>
    class UnaryAdapter :
      public UnaryInterface< t_ArgumentType, t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      //! first function a->b
      util::ShPtr< UnaryInterface< t_ArgumentType, t_IntermediateType> > m_FunctionAB;

      //! second Function b->c
      util::ShPtr< UnaryInterface< t_IntermediateType, t_ResultType> > m_FunctionBC;

      //! scheme to be used
      std::string m_Scheme;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from Function a->b and b->c
      //! @param SP_FUNCTION_AB ShPtr to Function a->b
      //! @param SP_FUNCTION_BC ShPtr to Function b->c
      UnaryAdapter
      (
        const util::ShPtr< UnaryInterface< t_ArgumentType    , t_IntermediateType> > &SP_FUNCTION_AB,
        const util::ShPtr< UnaryInterface< t_IntermediateType, t_ResultType> > &SP_FUNCTION_BC
      ) :
        m_FunctionAB( SP_FUNCTION_AB),
        m_FunctionBC( SP_FUNCTION_BC),
        m_Scheme( GetClassIdentifier())
      {
      }

      //! @brief construct from Function a->b and b->c
      //! @param SP_FUNCTION_AB ShPtr to Function a->b
      //! @param SP_FUNCTION_BC ShPtr to Function b->c
      //! @param SCHEME Scheme to be used
      UnaryAdapter
      (
        const util::ShPtr< UnaryInterface< t_ArgumentType    , t_IntermediateType> > &SP_FUNCTION_AB,
        const util::ShPtr< UnaryInterface< t_IntermediateType, t_ResultType> > &SP_FUNCTION_BC,
        const std::string &SCHEME
      ) :
        m_FunctionAB( SP_FUNCTION_AB),
        m_FunctionBC( SP_FUNCTION_BC),
        m_Scheme( SCHEME)
      {
      }

      //! virtual copy constructor
      UnaryAdapter< t_ArgumentType, t_IntermediateType, t_ResultType> *Clone() const
      {
        return new UnaryAdapter< t_ArgumentType, t_IntermediateType, t_ResultType>( *this);
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

      //! operator taking an ARGUMENT and returning a t_ResultType object
      t_ResultType operator()( t_ArgumentType &ARGUMENT) const
      {
        return m_FunctionBC->operator()( m_FunctionAB->operator()( ARGUMENT));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_FunctionAB, ISTREAM);
        io::Serialize::Read( m_FunctionBC, ISTREAM);
        io::Serialize::Read( m_Scheme    , ISTREAM);

        // end
        return ISTREAM;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_FunctionAB, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_FunctionBC, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Scheme    , OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class UnaryAdapter

  } // namespace function
} // namespace bcl

#endif //BCL_FUNCTION_UNARY_ADAPTER_H_
