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
#include "assemble/bcl_assemble_aa_neighbor_list.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AANeighborList::s_Instance
    (
      GetObjectInstances().AddInstance( new AANeighborList())
    );

    //! @brief access to default center aa for undefined list
    const biol::AABase &AANeighborList::GetDefaultCenterAA()
    {
      static const biol::AABase *s_default( new biol::AA());
      return *s_default;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief helper constructor to find a neighbor list in a AANeihborListCotainer from given amino acid
    //! @brief CENTER_AMINO_ACID amino acid for which the neighbor list is generated
    AANeighborList::AANeighborList( const biol::AABase &CENTER_AMINO_ACID) :
      m_CenterAminoAcid( &CENTER_AMINO_ACID),
      m_Neighbors(),
      m_DistanceCutoff( util::GetUndefined< double>()),
      m_MinimalSequenceSeparation( util::GetUndefined< size_t>()),
      m_ConsiderDifferentChain( true)
    {
    }

    //! @brief constructor from central AA, neighboring AAs and their distances
    //! @param CENTER_AMINO_ACID amino acid for which the neighbor list is generated
    //! @param AMINO_ACIDS_DISTANCES the list of amino acids in question along w/ their distances
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborList::AANeighborList
    (
      const biol::AABase &CENTER_AMINO_ACID,
      const NeighborContainerType &AMINO_ACIDS_DISTANCES,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_CenterAminoAcid( &CENTER_AMINO_ACID),
      m_Neighbors( AMINO_ACIDS_DISTANCES),
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
    }

    //! @brief construct from center amino acids, a list of amino acids, distance cutoff and minimal sequence separation
    //! @param CENTER_AMINO_ACID amino acid for which the neighbor list is generated
    //! @param AMINO_ACIDS the list of amino acids in question
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborList::AANeighborList
    (
      const biol::AABase &CENTER_AMINO_ACID,
      const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_CenterAminoAcid( &CENTER_AMINO_ACID),
      m_Neighbors(),
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
      // create reference to biol::Atom "atom_a" and initialize with the "ATOM_TYPE" atom of "AMINO_ACID"
      const biol::Atom &atom_a( CENTER_AMINO_ACID.GetFirstSidechainAtom());
      const char center_aa_chain_id( CENTER_AMINO_ACID.GetChainID());

      double square_distance_cutoff( math::Sqr( m_DistanceCutoff));
      // iterate through the residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr( AMINO_ACIDS.Begin()), aa_itr_end( AMINO_ACIDS.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        if( !m_ConsiderDifferentChain && ( center_aa_chain_id != ( *aa_itr)->GetChainID()))
        {
          continue;
        }

        // skip identical amino acids
        if( m_CenterAminoAcid == *aa_itr)
        {
          continue;
        }

        // check that amino acids are from the same chain and have a sufficient sequence distance
        const size_t seq_separation( biol::SequenceSeparation( CENTER_AMINO_ACID, **aa_itr));

        // seq separation can only be undefined if aa are from different chain or have the same seqid - which was checked
        if( util::IsDefined( seq_separation) && seq_separation < m_MinimalSequenceSeparation)
        {
          continue;
        }

        // create reference to biol::Atom "atom_b" and initialize with the "ATOM_TYPE" atom of "AA_ITR"
        const biol::Atom &atom_b( ( *aa_itr)->GetFirstSidechainAtom());

        // create double "distance" and initialize with the distance between "atom_a" and "atom_b"
        double sq_distance( biol::SquareDistance( atom_a, atom_b));

        // check that distance is defined and under cutoff
        if( sq_distance > square_distance_cutoff)
        {
          continue;
        }
        else if( sq_distance <= 0.0)
        {
          sq_distance = 0.0;
        }
        else if( !util::IsDefined( sq_distance))
        {
          continue;
        }

        // add the residue behind "AA_ITR" and "distance" to "current_neighbors"
        m_Neighbors.PushBack( storage::Pair< util::SiPtr< const biol::AABase>, double>( *aa_itr, math::Sqrt( sq_distance)));
      }
    }

    //! @brief construct from a neighbor list setup with parameters that'd make it a superset of the desired list
    //! @param PARENT the original neighborlist
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborList::AANeighborList
    (
      const AANeighborList &PARENT,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_CenterAminoAcid( PARENT.m_CenterAminoAcid),
      m_Neighbors(),
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
      if
      (
        m_DistanceCutoff == PARENT.m_DistanceCutoff
        && m_MinimalSequenceSeparation == PARENT.m_MinimalSequenceSeparation
        && ( m_ConsiderDifferentChain == PARENT.m_ConsiderDifferentChain)
      )
      {
        m_Neighbors = PARENT.m_Neighbors;
        return;
      }

      m_Neighbors.AllocateMemory( size_t( PARENT.m_Neighbors.GetSize() * m_DistanceCutoff / PARENT.m_DistanceCutoff));

      BCL_Assert
      (
        m_MinimalSequenceSeparation >= PARENT.m_MinimalSequenceSeparation
        && ( PARENT.m_ConsiderDifferentChain || !m_ConsiderDifferentChain)
        && m_DistanceCutoff <= PARENT.m_DistanceCutoff,
        "Cannot expand a neighborlist with more broad parameters: " +
        util::Format()( m_MinimalSequenceSeparation) + " min sep from "
        + util::Format()( PARENT.m_MinimalSequenceSeparation)
        + ( m_ConsiderDifferentChain && !PARENT.m_ConsiderDifferentChain
            ? ". Cannot get neighbors from other chain from neighborlist containing only same chain. " : " ")
        + ( m_DistanceCutoff >= PARENT.m_DistanceCutoff ? "Distance cutoff violation: from "
            + util::Format()( PARENT.m_DistanceCutoff) + " to " + util::Format()( m_DistanceCutoff) : ""
          )
      );
      const bool same_seq_sep( m_MinimalSequenceSeparation == PARENT.m_MinimalSequenceSeparation);
      const bool same_alt_chain( m_ConsiderDifferentChain == PARENT.m_ConsiderDifferentChain);
      const char chain_id( m_CenterAminoAcid->GetChainID());
      // iterate over all this neighbor list
      for( const_iterator itr( PARENT.m_Neighbors.Begin()), itr_end( PARENT.m_Neighbors.End()); itr != itr_end; ++itr)
      {
        if
        (
          itr->Second() <= m_DistanceCutoff
          && (
               same_seq_sep
               || biol::SequenceSeparation( *m_CenterAminoAcid, *itr->First()) >= m_MinimalSequenceSeparation
             )
          && ( same_alt_chain || itr->First()->GetChainID() == chain_id)
        )
        {
          m_Neighbors.PushBack( *itr);
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief push back another amino acid into the neighbor. Checks that the added AA/dist are appropriate for this list
    //!        should be done prior to this
    void AANeighborList::PushBack( const util::SiPtr< const biol::AABase> &AA, const double &DIST)
    {
      m_Neighbors.PushBack( storage::Pair< util::SiPtr< const biol::AABase>, double>( AA, DIST));
    }

    //! @brief return the intersect with the given neighbor list, the distances will indicate the differences from
    //! this list to the given list for the common residues
    //! @param AA_NEIGHBOR_LIST the neighbor list to calculate the intersection with
    //! @return Number of elements in the intersection of the lists
    size_t AANeighborList::IntersectionSize( const AANeighborList &AA_NEIGHBOR_LIST) const
    {
      size_t intersection_sz( 0);
      std::set< std::pair< char, int> > aa_pdb_id_set_a;
      // iterate over neighbors in this list
      for( const_iterator itr( m_Neighbors.Begin()), itr_end( m_Neighbors.End()); itr != itr_end; ++itr)
      {
        aa_pdb_id_set_a.insert( aa_pdb_id_set_a.end(), std::make_pair( itr->First()->GetChainID(), itr->First()->GetSeqID()));
      }
      for( const_iterator itr( AA_NEIGHBOR_LIST.m_Neighbors.Begin()), itr_end( AA_NEIGHBOR_LIST.m_Neighbors.End()); itr != itr_end; ++itr)
      {
        if
        (
          aa_pdb_id_set_a.find( std::make_pair( itr->First()->GetChainID(), itr->First()->GetSeqID()))
          != aa_pdb_id_set_a.end()
        )
        {
          ++intersection_sz;
        }
      }

      // end
      return intersection_sz;
    }

    //! @brief Remove neighbors with different chain id
    void AANeighborList::RemoveNeighborsWithDifferentChainID()
    {
      // nothing to do
      if( !m_ConsiderDifferentChain)
      {
        return;
      }

      const char chain_id( m_CenterAminoAcid->GetChainID());

      // iterate over all neighbors
      storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> > new_neigh;
      new_neigh.AllocateMemory( m_Neighbors.GetSize());

      // iterate over all neighbors
      for( auto itr( m_Neighbors.Begin()), itr_end( m_Neighbors.End()); itr != itr_end; ++itr)
      {
        if( itr->First()->GetChainID() == chain_id)
        {
          new_neigh.PushBack( *itr);
        }
      }
      m_Neighbors.InternalData().swap( new_neigh.InternalData());

      // update
      m_ConsiderDifferentChain = false;
    }

    //! @brief decrease the distance cutoff
    //! @param DISTANCE_CUTOFF
    //! @return true, on success, false if DISTANCE_CUTOFF is larger than current distance cutoff
    bool AANeighborList::DecreaseDistanceCutoff( const double DISTANCE_CUTOFF)
    {
      // distance cutoff can only be reduced, if the current one is larger
      if( m_DistanceCutoff < DISTANCE_CUTOFF)
      {
        return false;
      }

      // nothing to do
      if( m_DistanceCutoff == DISTANCE_CUTOFF)
      {
        return true;
      }

      // iterate over all neighbors
      storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> > new_neigh;
      new_neigh.AllocateMemory( m_Neighbors.GetSize());

      // iterate over all neighbors
      for( auto itr( m_Neighbors.Begin()), itr_end( m_Neighbors.End()); itr != itr_end; ++itr)
      {
        if( itr->Second() <= DISTANCE_CUTOFF)
        {
          new_neigh.PushBack( *itr);
        }
      }

      // update
      m_DistanceCutoff = DISTANCE_CUTOFF;
      new_neigh.InternalData().swap( m_Neighbors.InternalData());

      // return success
      return true;
    }

    //! @brief increase minimal sequence separation
    //! @param MIN_SEQ_SEPARATION
    //! @return true, if increment was successful, flase if the new seq separtion is smaller than the current
    bool AANeighborList::IncreaseMinimalSequenceSeparation( const size_t MIN_SEQ_SEPARATION)
    {
      // sequence separation can only be increased, if the current one is smaller
      if( m_MinimalSequenceSeparation > MIN_SEQ_SEPARATION)
      {
        return false;
      }

      // nothing to do
      if( m_MinimalSequenceSeparation == MIN_SEQ_SEPARATION)
      {
        return true;
      }

      storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> > new_neigh;
      new_neigh.AllocateMemory( m_Neighbors.GetSize());

      // iterate over all neighbors
      for( auto itr( m_Neighbors.Begin()), itr_end( m_Neighbors.End()); itr != itr_end; ++itr)
      {
        const size_t seq_separation( biol::SequenceSeparation( *m_CenterAminoAcid, *itr->First()));

        // if chain id does not match, remove it
        if( !util::IsDefined( seq_separation) || seq_separation >= MIN_SEQ_SEPARATION)
        {
          new_neigh.PushBack( *itr);
        }
      }

      // update
      m_MinimalSequenceSeparation = MIN_SEQ_SEPARATION;

      // swap vectors
      new_neigh.InternalData().swap( m_Neighbors.InternalData());

      // return success
      return true;
    }

    //! @brief Remove residues located on the same SSE
    //! @param MODEL the protein model containing the SSEs
    void AANeighborList::PruneResiduesSameSSE( const ProteinModel &MODEL)
    {
      util::SiPtr< const SSE> sse( MODEL.GetSSE( *m_CenterAminoAcid));
      storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> > new_neigh;
      new_neigh.AllocateMemory( m_Neighbors.GetSize());

      // iterate over all neighbors
      for( auto itr( m_Neighbors.Begin()), itr_end( m_Neighbors.End()); itr != itr_end; ++itr)
      {
        // if chain id does not match, remove it
        if( sse != MODEL.GetSSE( *( itr->First())))
        {
          new_neigh.PushBack( *itr);
        }
      }

      // swap vectors
      new_neigh.InternalData().swap( m_Neighbors.InternalData());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AANeighborList::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_DistanceCutoff           , ISTREAM);
      io::Serialize::Read( m_MinimalSequenceSeparation, ISTREAM);
      io::Serialize::Read( m_ConsiderDifferentChain   , ISTREAM);
//      io::Serialize::Read( m_CenterAminoAcid          , ISTREAM);
      io::Serialize::Read( m_Neighbors              , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AANeighborList::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_DistanceCutoff           , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_MinimalSequenceSeparation, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ConsiderDifferentChain   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CenterAminoAcid          , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Neighbors                , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
