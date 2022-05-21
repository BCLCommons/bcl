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
#include "score/bcl_score_aa_sequence.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container_generator_sse.h"
#include "biol/bcl_biol_aa_base.h"
#include "fold/bcl_fold_default_flags.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASequence::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequence( util::ShPtr< AAPairDistanceInterface>(), true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequence::AASequence() :
      m_ScoreAAPair(),
      m_Normalize( false),
      m_AANeighborListContainerGenerator()
    {
    }

    //! @brief construct from ShPtr to amino acid pair potential
    //! @param SP_AA_PAIR_POTENTIAL ShPtr to the amino acid potential to be used
    AASequence::AASequence
    (
      const util::ShPtr< AAPairDistanceInterface> &SP_AA_PAIR_POTENTIAL,
      const bool NORMALIZE
    ) :
      m_ScoreAAPair( SP_AA_PAIR_POTENTIAL),
      m_Normalize( NORMALIZE),
      m_AANeighborListContainerGenerator
      (
        SP_AA_PAIR_POTENTIAL.IsDefined() ?
          assemble::AANeighborListContainerGeneratorSSE::AANeighborListGenerator
          (
            SP_AA_PAIR_POTENTIAL->GetDistanceCutoff(),
            SP_AA_PAIR_POTENTIAL->GetMinimalSequenceSeparation(),
            SP_AA_PAIR_POTENTIAL->GetConsiderDifferentChain(),
            true
          )
          :
          util::ShPtr< math::FunctionInterfaceSerializable< assemble::SSE, assemble::AANeighborListContainer> >()
      )
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AASequence copied from this one
    AASequence *AASequence::Clone() const
    {
      return new AASequence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &AASequence::GetScheme() const
    {
      return m_ScoreAAPair->GetScheme();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief score aa distances between all pairs of amino acids within a single SSE
    //! @param THIS_SSE SSE of interest
    //! @param MEMBRANE membrane object
    //! @return pair of overall interaction potential and number of scored entities
    storage::Pair< double, size_t> AASequence::operator()
    (
      const assemble::SSE &THIS_SSE,
      const biol::Membrane &MEMBRANE
    ) const
    {
      double score( 0);

      size_t number_pair_scores( 0);

      // create neighborlistcontainer between both sequences
      const assemble::AANeighborListContainer neighbor_list_container( m_AANeighborListContainerGenerator->operator ()( THIS_SSE));

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
      return storage::Pair< double, size_t>( score, number_pair_scores);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read from std::istream
    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &AASequence::Read( std::istream &ISTREAM)
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
    std::ostream &AASequence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ScoreAAPair, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Normalize, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param THIS_SSE SSE of interest
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &AASequence::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &THIS_SSE,
      std::ostream &OSTREAM
    ) const
    {
      // create neighborlistcontainer between both sequences
      const assemble::AANeighborListContainer neighbor_list_container( m_AANeighborListContainerGenerator->operator ()( THIS_SSE));

      biol::AALessThanSeqID aa_lt;
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
          if( aa_lt( current_aa, neigh_aa))
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
