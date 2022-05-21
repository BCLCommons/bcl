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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "score/bcl_score_protein_model_wrapper.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelWrapper::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelWrapper())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelWrapper::ProteinModelWrapper() :
      m_Function(),
      m_ScoreType( e_Undefined),
      m_ReadableScheme()
    {
    }

    //! @brief construct from member
    //! @param SP_FUNCTION function to wrap
    //! @brief construct from members
    //! @param SP_FUNCTION function to wrap
    //! @param TYPE score type
    //! @param READABLE_SCHEME readable scheme
    ProteinModelWrapper::ProteinModelWrapper
    (
      const util::ShPtr< math::FunctionInterfaceSerializable< assemble::ProteinModel, double> > &SP_FUNCTION,
      const Type &TYPE,
      const std::string &READABLE_SCHEME,
      const std::string &SCHEME
    ) :
      m_Function( SP_FUNCTION),
      m_ScoreType( TYPE),
      m_ReadableScheme( READABLE_SCHEME),
      m_Scheme( SCHEME)
    {
      if( m_ReadableScheme.empty())
      {
        m_ReadableScheme = m_Function->GetScheme();
      }
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelWrapper
    ProteinModelWrapper *ProteinModelWrapper::Clone() const
    {
      return new ProteinModelWrapper( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelWrapper::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelWrapper::GetScheme() const
    {
      return m_Scheme.empty() ? m_Function->GetScheme() : m_Scheme;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @return function value of the given argument
    double ProteinModelWrapper::operator()( const assemble::ProteinModel &ARGUMENT) const
    {
      return m_Function->operator ()( ARGUMENT);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelWrapper::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelWrapper::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

    //! @brief write detailed scheme and values to OSTREAM
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelWrapper::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &ARGUMENT,
      std::ostream &OSTREAM
    ) const
    {
      return m_Function->WriteDetailedSchemeAndValues( ARGUMENT, OSTREAM);
    }

  } // namespace score
} // namespace bcl
