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
#include "score/bcl_score_protein_model_sse_pairs.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSEPairs::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSEPairs())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    ProteinModelSSEPairs::ProteinModelSSEPairs() :
      m_SpScorePair(),
      m_Normalize( false)
    {
    }

    //! @brief construct from Pair function
    //! @param SP_PAIR_FUNCTION binary function to score a pair of sses
    //! @param NORMALIZE if true, final score will be normalized by the number of sses/sequences in the protein model
    //! @param SCORE_TYPE score type
    //! @param READABLE_SCHEME scheme that is more human readable
    ProteinModelSSEPairs::ProteinModelSSEPairs
    (
      const util::ShPtr< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> > &SP_PAIR_FUNCTION,
      const bool NORMALIZE,
      const ProteinModel::Type &SCORE_TYPE,
      const std::string &READABLE_SCHEME
    ) :
      m_SpScorePair( *SP_PAIR_FUNCTION),
      m_Normalize( NORMALIZE),
      m_ScoreType( SCORE_TYPE),
      m_ReadableScheme( READABLE_SCHEME)
    {
      if( m_ReadableScheme.empty())
      {
        m_ReadableScheme = GetScheme();
      }
    }

    //! virtual copy constructor
    ProteinModelSSEPairs *ProteinModelSSEPairs::Clone() const
    {
      return new ProteinModelSSEPairs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSEPairs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelSSEPairs::GetScheme() const
    {
      return m_SpScorePair->GetScheme();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores the Protein model
    //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
    //! @return score
    double ProteinModelSSEPairs::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      //instantiate the scoresum
      double score( 0.0);

      //instantiate a util::SiPtrVector to the elements in the ProteinModel
      const util::SiPtrVector< const assemble::SSE> all_elements( PROTEIN_MODEL.GetSSEs());

      //iterate over all possible pairs of elements
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator element_itr_a( all_elements.Begin()),
          element_itr_end( all_elements.End());
        element_itr_a != element_itr_end;
        ++element_itr_a
      )
      {
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator element_itr_b( element_itr_a + 1);
          element_itr_b != element_itr_end;
          ++element_itr_b
        )
        {
          //call the pairscore
          score += m_SpScorePair->operator()( **element_itr_a, **element_itr_b);
        }
      }

      // if protein model has sses and normalize flag is set
      if( m_Normalize && !all_elements.IsEmpty())
      {
        return score / all_elements.GetSize();
      }

      //end
      return score;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ProteinModelSSEPairs::GetAlias() const
    {
      static const std::string s_name( "ProteinModelSSEPairs");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSEPairs::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scores interactions between SSE pairs.");
      serializer.AddInitializer
      (
        "score function",
        "function to score residue-residue interactions between the two SSEs",
        io::Serialization::GetAgent( &m_SpScorePair)
      );
      serializer.AddInitializer
      (
        "normalize",
        "normalize score by number of residues",
        io::Serialization::GetAgent( &m_Normalize)
      );

      return serializer;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelSSEPairs::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      //instantiate a util::SiPtrVector to the elements in the ProteinModel
      const util::SiPtrVector< const assemble::SSE> all_elements( PROTEIN_MODEL.GetSSEs());

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator element_itr_a( all_elements.Begin()),
          element_itr_end( all_elements.End());
        element_itr_a != element_itr_end;
        ++element_itr_a
      )
      {
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator element_itr_b( element_itr_a + 1);
          element_itr_b != element_itr_end;
          ++element_itr_b
        )
        {
          //write the detailed scheme and value for the pair score
          m_SpScorePair->WriteDetailedSchemeAndValues( **element_itr_a, **element_itr_b, OSTREAM);
        }
      }

      //end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &ProteinModelSSEPairs::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SpScorePair, ISTREAM);
      io::Serialize::Read( m_Normalize, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &ProteinModelSSEPairs::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SpScorePair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Normalize, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
