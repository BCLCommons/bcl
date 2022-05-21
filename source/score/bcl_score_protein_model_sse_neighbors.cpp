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
#include "score/bcl_score_protein_model_sse_neighbors.h"

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
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSENeighbors::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinModelSSENeighbors())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelSSENeighbors::ProteinModelSSENeighbors()
    {
    }

    //! @brief construct from Pair function
    //! @param SP_PAIR_FUNCTION binary function to score a pair of sses
    //! @param NORMALIZE if true, final score will be normalized by the number of sse pairs in the protein model
    ProteinModelSSENeighbors::ProteinModelSSENeighbors
    (
      const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> > &SP_PAIR_FUNCTION,
      const bool NORMALIZE
    ) :
      m_ScorePair( SP_PAIR_FUNCTION),
      m_Normalize( NORMALIZE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelSSENeighbors
    ProteinModelSSENeighbors *ProteinModelSSENeighbors::Clone() const
    {
      return new ProteinModelSSENeighbors( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSENeighbors::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelSSENeighbors::GetScheme() const
    {
      return m_ScorePair->GetScheme();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores the chain
    //! @param CHAIN the chain for which all neighbor scores are calculated
    //! @return score
    double ProteinModelSSENeighbors::operator()( const assemble::Chain &CHAIN) const
    {
      // sum of all scores
      double score( 0.0);

      // need at least two sses in the chain
      if( CHAIN.GetNumberSSEs() < 2)
      {
        return score;
      }

      size_t count( 0);

      // iterate through the sses of "CHAIN" to calculate their score sum
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr_a( CHAIN.GetData().Begin()),
          sse_itr_b( ++CHAIN.GetData().Begin()),
          sse_itr_end( CHAIN.GetData().End());
        sse_itr_b != sse_itr_end;
        ++sse_itr_a, ++sse_itr_b, ++count
      )
      {
        score += m_ScorePair->operator ()( **sse_itr_a, **sse_itr_b);
      }

      // normalize
      if( m_Normalize)
      {
        score /= double( count);
      }

      // return the score sum
      return score;
    }

    //! @brief operator that scores the Protein model
    //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
    //! @return score
    double ProteinModelSSENeighbors::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // score sum
      double score( 0.0);

      // iterate through the chains of PROTEIN_MODEL
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        score += operator()( **chain_itr);
      }

      // return score
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelSSENeighbors::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ScorePair, ISTREAM);
      io::Serialize::Read( m_Normalize, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelSSENeighbors::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ScorePair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Normalize, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelSSENeighbors::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      // iterate through the chains of PROTEIN_MODEL
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        if( ( *chain_itr)->GetData().IsEmpty())
        {
          continue;
        }
        // iterate through the sses of "CHAIN" to calculate their score sum
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr_a( ( *chain_itr)->GetData().Begin()),
            sse_itr_b( ++( *chain_itr)->GetData().Begin()),
            sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr_b != sse_itr_end;
          ++sse_itr_a, ++sse_itr_b
        )
        {
          m_ScorePair->WriteDetailedSchemeAndValues( **sse_itr_a, **sse_itr_b, OSTREAM);
        }
      }

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
