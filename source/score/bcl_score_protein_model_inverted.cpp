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
#include "score/bcl_score_protein_model_inverted.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelInverted::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelInverted())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelInverted::ProteinModelInverted() :
      m_Score(),
      m_Inverter(),
      m_Scheme()
    {
    }

    //! @brief constructor from a score function, a inverter and a scheme
    //! @param SP_SCORE ShPtr to ProteinModel scoring function to be used
    //! @param SP_INVERTER ShPtr to ProteinModelInverter to be used
    //! @param SCHEME Scheme to be used
    //! @param SCORE_TYPE score type
    //! @param READABLE_SCHEME scheme that is more human readable
    ProteinModelInverted::ProteinModelInverted
    (
      const util::ShPtr< ProteinModel> &SP_SCORE,
      const util::ShPtr< assemble::ProteinModelInverter> &SP_INVERTER,
      const std::string &SCHEME,
      const ProteinModel::Type &SCORE_TYPE,
      const std::string &READABLE_SCHEME
    ) :
      m_Score( SP_SCORE),
      m_Inverter( SP_INVERTER),
      m_Scheme( SCHEME),
      m_ScoreType( SCORE_TYPE),
      m_ReadableScheme( READABLE_SCHEME)
    {
      if( m_ReadableScheme.empty())
      {
        m_ReadableScheme = GetScheme();
      }
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelInverted
    ProteinModelInverted *ProteinModelInverted::Clone() const
    {
      return new ProteinModelInverted( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelInverted::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief invert the given ProteinModel and score it
    //! @param PROTEIN_MODEL ProteinModel to be inverted and scored
    //! @return score of the inverted model
    double ProteinModelInverted::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // calculate the inverted model
      util::ShPtr< assemble::ProteinModel> sp_inverted_model( m_Inverter->GetInvertedModel( PROTEIN_MODEL));

      // make sure it is defined
      BCL_Assert( sp_inverted_model.IsDefined(), "ProteinModelInverter returned empty ShPtr");

      // calculate score and return it
      return m_Score->operator ()( *sp_inverted_model);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelInverted::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_Inverter, ISTREAM);
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelInverted::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Inverter, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
