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
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AANeighborListContainer::s_Instance
    (
      GetObjectInstances().AddInstance
      (
        new AANeighborListContainer
        (
          util::SiPtrVector< const biol::AABase>(), util::GetUndefined< double>(), util::GetUndefined< size_t>(), false
        )
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from distance cutoff and minimal sequence separation
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainer::AANeighborListContainer
    (
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
    }

    //! @brief construct from a list of amino acids, distance cutoff and minimal sequence separation
    //! @param AMINO_ACIDS the list of amino acids in question
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainer::AANeighborListContainer
    (
      const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
      iterator itr_guess( m_NeighborLists.Begin());
      // iterate through the residues
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr( AMINO_ACIDS.Begin()), aa_itr_end( AMINO_ACIDS.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // insert AANeighborList into the set
        itr_guess = m_NeighborLists.InsertElement
        (
          itr_guess,
          std::make_pair
          (
            *aa_itr,
            AANeighborList( **aa_itr, util::SiPtrVector< const biol::AABase>(), m_DistanceCutoff, m_MinimalSequenceSeparation, m_ConsiderDifferentChain)
          )
        );
      }
      if( AMINO_ACIDS.GetSize() > size_t( 1))
      {
        VoxelGridAA voxel_grid( m_DistanceCutoff);
        voxel_grid.SetObjects( AMINO_ACIDS);
        auto neighbors( voxel_grid.GetNeighbors( m_DistanceCutoff));
        auto itr_last_found( m_NeighborLists.Begin()), itr_last_found_b( m_NeighborLists.Begin());
        for( auto itr_neigh( neighbors.Begin()), itr_neigh_end( neighbors.End()); itr_neigh != itr_neigh_end; ++itr_neigh)
        {
          if( !m_ConsiderDifferentChain && itr_neigh->First()->GetChainID() != itr_neigh->Second()->GetChainID())
          {
            continue;
          }
          if( biol::SequenceSeparation( *itr_neigh->First(), *itr_neigh->Second()) < m_MinimalSequenceSeparation)
          {
            continue;
          }
          if( itr_neigh->First() != itr_last_found->first)
          {
            itr_last_found = m_NeighborLists.InternalData().find( *itr_neigh->First());
          }
          itr_last_found->second.m_Neighbors.PushBack( storage::Pair< util::SiPtr< const biol::AABase>, double>( itr_neigh->Second(), itr_neigh->Third()));
          m_NeighborLists[ *itr_neigh->Second()].m_Neighbors.PushBack( storage::Pair< util::SiPtr< const biol::AABase>, double>( itr_neigh->First(), itr_neigh->Third()));
        }
      }
    }

    //! @brief construct from two lists of amino acids, distance cutoff and minimal sequence separation
    //! @details only inter pairwise distances and neighbors are calculated
    //! @param AMINO_ACIDS_A, AMINO_ACIDS_B the list of amino acids in question
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainer::AANeighborListContainer
    (
      const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS_A,
      const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS_B,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
      iterator itr_guess( m_NeighborLists.Begin());
      // iterate through the residues a
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr( AMINO_ACIDS_A.Begin()), aa_itr_end( AMINO_ACIDS_A.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // insert AANeighborList into the set
        itr_guess = m_NeighborLists.InsertElement
        (
          itr_guess,
          std::make_pair
          (
            *aa_itr,
            AANeighborList
            (
              **aa_itr,
              util::SiPtrVector< const biol::AABase>(),
              m_DistanceCutoff,
              m_MinimalSequenceSeparation,
              m_ConsiderDifferentChain
            )
          )
        );
      }

      // iterate through the residues b
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr( AMINO_ACIDS_B.Begin()), aa_itr_end( AMINO_ACIDS_B.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // insert AANeighborList into the set
        itr_guess = m_NeighborLists.InsertElement
        (
          itr_guess,
          std::make_pair
          (
            *aa_itr,
            AANeighborList
            (
              **aa_itr,
              util::SiPtrVector< const biol::AABase>(),
              m_DistanceCutoff,
              m_MinimalSequenceSeparation,
              m_ConsiderDifferentChain
            )
          )
        );
      }

      VoxelGridAA voxel_grid( m_DistanceCutoff);
      voxel_grid.SetObjects( AMINO_ACIDS_B);
      for( auto itr_a( AMINO_ACIDS_A.Begin()), itr_a_end( AMINO_ACIDS_A.End()); itr_a != itr_a_end; ++itr_a)
      {
        auto neighbors( voxel_grid.GetNeighbors( **itr_a, m_DistanceCutoff));
        AANeighborList &neigh_a( m_NeighborLists[ **itr_a]);
        const char chain_id( ( *itr_a)->GetChainID());
        for( auto itr_neigh( neighbors.Begin()), itr_neigh_end( neighbors.End()); itr_neigh != itr_neigh_end; ++itr_neigh)
        {
          if( !m_ConsiderDifferentChain && chain_id != itr_neigh->First()->GetChainID())
          {
            continue;
          }
          if( biol::SequenceSeparation( **itr_a, *itr_neigh->First()) < m_MinimalSequenceSeparation)
          {
            continue;
          }
          neigh_a.m_Neighbors.PushBack( storage::Pair< util::SiPtr< const biol::AABase>, double>( itr_neigh->First(), itr_neigh->Second()));
          m_NeighborLists[ *itr_neigh->First()].m_Neighbors.PushBack( storage::Pair< util::SiPtr< const biol::AABase>, double>( *itr_a, itr_neigh->Second()));
        }
      }
    }

    //! @brief Copy constructor
    //! @param AMINO_ACIDS the list of amino acids in question
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    AANeighborListContainer::AANeighborListContainer
    (
      const AANeighborListContainer &AMINO_ACIDS,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
      for( auto itr( AMINO_ACIDS.Begin()), itr_end( AMINO_ACIDS.End()); itr != itr_end; ++itr)
      {
        m_NeighborLists.InternalData().insert
        (
          m_NeighborLists.End(),
          std::pair< util::SiPtr< const biol::AABase>, AANeighborList>
          (
            itr->first,
            AANeighborList( itr->second, m_DistanceCutoff, m_MinimalSequenceSeparation, m_ConsiderDifferentChain)
          )
        );
      }
    }

    //! @brief construct from already-found amino acids (e.g. using assemble::VoxelGridAA.GetSSEConnections())
    //! @param NEIGHBORS List of neighbors
    AANeighborListContainer::AANeighborListContainer
    (
      const storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> > &NEIGHBORS,
      const util::SiPtrVector< const biol::AABase> &AAS,
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    ) :
      m_DistanceCutoff( DISTANCE_CUTOFF),
      m_MinimalSequenceSeparation( MIN_SEQ_SEPARATION),
      m_ConsiderDifferentChain( CONSIDER_DIFFERENT_CHAIN)
    {
      for( auto itr( AAS.Begin()), itr_end( AAS.End()); itr != itr_end; ++itr)
      {
        m_NeighborLists.InternalData().insert
        (
          m_NeighborLists.End(),
          std::pair< util::SiPtr< const biol::AABase>, AANeighborList>
          (
            *itr,
            AANeighborList( **itr, AANeighborList::NeighborContainerType(), m_DistanceCutoff, m_MinimalSequenceSeparation, m_ConsiderDifferentChain)
          )
        );
      }
      for( auto itr( NEIGHBORS.Begin()), itr_end( NEIGHBORS.End()); itr != itr_end; ++itr)
      {
        m_NeighborLists[ itr->First()].PushBack( itr->Second(), itr->Third());
        m_NeighborLists[ itr->Second()].PushBack( itr->First(), itr->Third());
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the total number of neighbors
    //! @return total number of neighbors
    size_t AANeighborListContainer::GetNumberNeighbors() const
    {
      // initialize count
      size_t count( 0);

      // iterate over all the lists
      for( const_iterator itr( m_NeighborLists.Begin()), itr_end( m_NeighborLists.End()); itr != itr_end; ++itr)
      {
        // sum up the number of members in each list
        count += itr->second.GetSize();
      }

      // end
      return count;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief prune the AANeighborListContainer
    //! @param DISTANCE_CUTOFF above which neighbors are not stored
    //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
    //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
    void AANeighborListContainer::Prune
    (
      const double DISTANCE_CUTOFF,
      const size_t MIN_SEQ_SEPARATION,
      const bool CONSIDER_DIFFERENT_CHAIN
    )
    {
      if
      (
        DISTANCE_CUTOFF == m_DistanceCutoff
        && MIN_SEQ_SEPARATION == m_MinimalSequenceSeparation
        && ( CONSIDER_DIFFERENT_CHAIN || !m_ConsiderDifferentChain)
      )
      {
        return;
      }
      BCL_Assert
      (
        m_ConsiderDifferentChain >= CONSIDER_DIFFERENT_CHAIN,
        "cannot consider different chain id, since this list, does not contain this information"
      );
      // optimization : if there is only one chain in the neighborlist, there will be no need for removal
      // of other chains
      if( !CONSIDER_DIFFERENT_CHAIN && m_ConsiderDifferentChain && !m_NeighborLists.IsEmpty())
      {
        const char chain_id_first( m_NeighborLists.Begin()->first->GetChainID());
        bool all_same_chain_id( true);
        for( const_iterator itr( m_NeighborLists.Begin()), itr_end( m_NeighborLists.End()); itr != itr_end; ++itr)
        {
          if( itr->first->GetChainID() != chain_id_first)
          {
            all_same_chain_id = false;
            break;
          }
        }
        if( all_same_chain_id)
        {
          m_ConsiderDifferentChain = CONSIDER_DIFFERENT_CHAIN;
          for( iterator itr( m_NeighborLists.Begin()), itr_end( m_NeighborLists.End()); itr != itr_end; ++itr)
          {
            itr->second.m_ConsiderDifferentChain = false;
          }
        }
      }
      size_t old_size( 0), new_size( 0);
      // iterate over all this neighbor list
      for( iterator itr( m_NeighborLists.Begin()), itr_end( m_NeighborLists.End()); itr != itr_end; ++itr)
      {
        AANeighborList &this_list( itr->second);
        old_size += this_list.GetSize();
        BCL_Assert
        (
          this_list.DecreaseDistanceCutoff( DISTANCE_CUTOFF), "could not reduce distance cutoff"
        );
        BCL_Assert
        (
          this_list.IncreaseMinimalSequenceSeparation( MIN_SEQ_SEPARATION),
          "could not increase minimal sequence separation"
        );
        if( !CONSIDER_DIFFERENT_CHAIN)
        {
          this_list.RemoveNeighborsWithDifferentChainID();
        }
        new_size += this_list.GetSize();
      }
      BCL_MessageDbg
      (
        "$$$: Pruning neighbor list container from " + util::Format()( m_DistanceCutoff)
        + " to " + util::Format()( DISTANCE_CUTOFF) + " A. min seq sep "
        + util::Format()( m_MinimalSequenceSeparation) + " -> " + util::Format()( MIN_SEQ_SEPARATION)
        + " diff chain: " + util::Format()( m_ConsiderDifferentChain) + " -> " + util::Format()( CONSIDER_DIFFERENT_CHAIN)
        + " neighbors decrease: " + util::Format()( old_size) + " -> " + util::Format()( new_size)
      );
      // update parameters
      m_DistanceCutoff = DISTANCE_CUTOFF;
      m_MinimalSequenceSeparation = MIN_SEQ_SEPARATION;
      m_ConsiderDifferentChain = CONSIDER_DIFFERENT_CHAIN;
    }

    //! @brief Remove residues located on the same SSE
    //! @param MODEL the protein model containing the SSEs
    void AANeighborListContainer::PruneResiduesSameSSE( const ProteinModel &MODEL)
    {
      for( iterator itr( m_NeighborLists.Begin()), itr_end( m_NeighborLists.End()); itr != itr_end; ++itr)
      {
        AANeighborList &this_list( itr->second);
        this_list.PruneResiduesSameSSE( MODEL);
      }
    }

    //! @brief return the intersect size with the given AANeighborListContainer
    //! @param CONTAINER The container with which the intersect is going to be calculated
    //! @return number of overlapping edges
    size_t AANeighborListContainer::IntersectionSize( const AANeighborListContainer &CONTAINER) const
    {
      size_t n_overlap( 0);

      // if the parameters are different
      if
      (
        m_DistanceCutoff != CONTAINER.m_DistanceCutoff ||
        m_MinimalSequenceSeparation != CONTAINER.m_MinimalSequenceSeparation ||
        m_ConsiderDifferentChain != CONTAINER.m_ConsiderDifferentChain
      )
      {
        // return an empty container
        BCL_MessageCrt
        (
          "Can not calculate intersection for given containers because they have different parameters"
        );
        return 0;
      }

      // iterate over the list of members in this container
      for( const_iterator itr( m_NeighborLists.Begin()), itr_end( m_NeighborLists.End()); itr != itr_end; ++itr)
      {
        // search for the list with the same center amino acid in the given container
        const const_iterator itr_other( CONTAINER.m_NeighborLists.Find( itr->first));

        // if not found continue
        if( itr_other != CONTAINER.End())
        {
          // calculate the intersection between two neighbor lists
          n_overlap += itr->second.IntersectionSize( itr_other->second);
        }

      }

      // end
      return n_overlap;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AANeighborListContainer::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_DistanceCutoff           , ISTREAM);
      io::Serialize::Read( m_MinimalSequenceSeparation, ISTREAM);
      io::Serialize::Read( m_ConsiderDifferentChain   , ISTREAM);
      io::Serialize::Read( m_NeighborLists            , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AANeighborListContainer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_DistanceCutoff           , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_MinimalSequenceSeparation, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ConsiderDifferentChain   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NeighborLists            , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
