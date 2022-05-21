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
#include "score/bcl_score_aa_sequence_pair.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_sse_pair.h"
#include "biol/bcl_biol_aa_base.h"
#include "fold/bcl_fold_default_flags.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASequencePair::s_Instance
    (
      util::Enumerated< math::BinaryFunctionInterfaceSerializable< assemble::SSE, assemble::SSE, double> >::AddInstance
      (
        new AASequencePair
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequencePair::AASequencePair() :
      m_ScoreAAPair(),
      m_Normalize(),
      m_AANeighborListContainerGenerator()
    {
    }

    //! @brief construct from ShPtr to amino acid pair potential
    //! @param AA_PAIR_POTENTIAL the amino acid potential to be used
    AASequencePair::AASequencePair
    (
      const AAPairDistanceInterface &AA_PAIR_POTENTIAL,
      const bool NORMALIZE
    ) :
      m_ScoreAAPair( AA_PAIR_POTENTIAL),
      m_Normalize( NORMALIZE),
      m_AANeighborListContainerGenerator
      (
        assemble::AANeighborListContainerGeneratorSSEPair::AANeighborListGenerator
        (
          AA_PAIR_POTENTIAL.GetDistanceCutoff(),
          AA_PAIR_POTENTIAL.GetMinimalSequenceSeparation(),
          AA_PAIR_POTENTIAL.GetConsiderDifferentChain(),
          false
        )
      )
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AASequencePair copied from this one
    AASequencePair *AASequencePair::Clone() const
    {
      return new AASequencePair( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequencePair::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AASequencePair::GetScheme() const
    {
      return m_ScoreAAPair->GetScheme();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &AASequencePair::GetAlias() const
    {
      static const std::string s_name( "AASequencePair");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AASequencePair::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Scores pairs of SSEs using pairwise residue interactions.");
      serializer.AddInitializer
      (
        "scoring function",
        "function to score the residue-residue interactions",
        io::Serialization::GetAgent( &m_ScoreAAPair)
      );
      serializer.AddInitializer
      (
        "normalize",
        "normalize by number of scored residue pairs",
        io::Serialization::GetAgent( &m_Normalize)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the members of this object from the given label
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool AASequencePair::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_AANeighborListContainerGenerator =
      (
        assemble::AANeighborListContainerGeneratorSSEPair::AANeighborListGenerator
        (
          m_ScoreAAPair->GetDistanceCutoff(),
          m_ScoreAAPair->GetMinimalSequenceSeparation(),
          m_ScoreAAPair->GetConsiderDifferentChain(),
          false
        )
       );

      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief score aa distances between all pairs of amino acids between the two SSEs (not within either one)
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @return overall interaction potential
    double AASequencePair::operator()
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B
    ) const
    {
      double score( 0);

      size_t number_pair_scores( 0);

      // create neighborlistcontainer between both sequences
      const assemble::AANeighborListContainer neighbor_list_container
      (
        m_AANeighborListContainerGenerator->operator ()( SSE_A, SSE_B)
      );

      // iterate over all amino acids
      for
      (
        assemble::AANeighborListContainer::const_iterator
          aa_itr( neighbor_list_container.Begin()), aa_itr_end( neighbor_list_container.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // current aa
        const biol::AABase &current_aa( *aa_itr->second.GetCenterAminoAcid());

        // iterate over all neighbors
        for
        (
          assemble::AANeighborList::const_iterator neigh_itr( aa_itr->second.Begin()), neigh_itr_end( aa_itr->second.End());
          neigh_itr != neigh_itr_end;
          ++neigh_itr
        )
        {
          const biol::AABase &neigh_aa( *neigh_itr->First());

          // score only one pair
          if( biol::AALessThanSeqID()( current_aa, neigh_aa))
          {
            continue;
          }

          // score distance between the two AAs
          const double current_score( m_ScoreAAPair->operator()( current_aa, neigh_aa, neigh_itr->Second()));

          if( util::IsDefined( current_score) && current_score != 0)
          {
            score += current_score;
            ++number_pair_scores;
          }
        }
      }

      // normalize
      if( m_Normalize && number_pair_scores > 0)
      {
        score /= double( number_pair_scores);
      }

      // end
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read from std::istream
    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AASequencePair::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ScoreAAPair, ISTREAM);
      io::Serialize::Read( m_Normalize, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @return returns the output stream
    std::ostream &AASequencePair::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ScoreAAPair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Normalize, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &AASequencePair::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {
      // create neighborlistcontainer between both sequences
      const assemble::AANeighborListContainer neighbor_list_container( m_AANeighborListContainerGenerator->operator ()( SSE_A, SSE_B));

      // iterate over all amino acids
      for
      (
        assemble::AANeighborListContainer::const_iterator
          aa_itr( neighbor_list_container.Begin()), aa_itr_end( neighbor_list_container.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // current aa
        const biol::AABase &current_aa( *aa_itr->second.GetCenterAminoAcid());

        // iterate over all neighbors
        for
        (
          assemble::AANeighborList::const_iterator neigh_itr( aa_itr->second.Begin()), neigh_itr_end( aa_itr->second.End());
          neigh_itr != neigh_itr_end;
          ++neigh_itr
        )
        {
          const biol::AABase &neigh_aa( *neigh_itr->First());

          // score only one pair
          if( biol::AALessThanSeqID()( current_aa, neigh_aa))
          {
            continue;
          }

          // score distance between the two AAs
          m_ScoreAAPair->WriteDetailedSchemeAndValues( current_aa, neigh_aa, OSTREAM);
        }
      }

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
