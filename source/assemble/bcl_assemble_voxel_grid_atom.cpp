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
#include "assemble/bcl_assemble_voxel_grid_atom.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> VoxelGridAtom::s_Instance
    (
      GetObjectInstances().AddInstance( new VoxelGridAtom())
    );

    //! @brief set the group ids for the atoms (indexed by atom PDB_ID). Atoms in the same group will not be reported
    //!        as in-contact.
    //! @param GROUP_IDS a vector containing group ids for all atoms of interest.
    void VoxelGridAtom::SetGroupIDs( const storage::Vector< size_t> &GROUP_IDS)
    {
      BCL_Assert( m_MaxPDBID <= GROUP_IDS.GetSize(), "Inadequate # of group IDs provided!");
      m_PdbIdToSeqID = GROUP_IDS;
      m_PdbIdToGroupIdPtr = util::ToSiPtr( m_PdbIdToSeqID);
    }

    //! @brief Extract an atoms coordinates (Ha2 for Glycine, Cb-Atom otherwise)
    //! @param ATOM pointer to the atom object
    //! @return Reference to the coordinates. null-ptr if coordinates are undefined
    util::SiPtr< const linal::Vector3D> VoxelGridAtom::ExtractPosition( const util::SiPtr< const biol::Atom> &ATOM) const
    {
      // Check if it's coordinates are defined
      if( ATOM->GetCoordinates().IsDefined() && ATOM->GetPdbID() >= 0 && ATOM->GetPdbID() < 1000000)
      { // if so return them
        m_MaxPDBID = std::max( size_t( ATOM->GetPdbID()), m_MaxPDBID);
        return util::SiPtr< const linal::Vector3D>( ATOM->GetCoordinates());
      }
      return util::SiPtr< const linal::Vector3D>();
    }

    //! @brief Update the VoxelGrid with new data
    //! @param NEW_DATA SiPtrVector of the new data
    void VoxelGridAtom::SetObjects( const util::SiPtrVector< const biol::Atom> &NEW_DATA)
    {
      m_MaxPDBID = 0;
      util::VoxelGrid< biol::Atom>::SetObjects( NEW_DATA);
      // reset group ids
      m_PdbIdToGroupIdPtr = util::SiPtr< const storage::Vector< size_t> >();
      m_PdbIdToSSEID.Reset();
      m_PdbIdToSeqID.Reset();
      m_NumberSSEs = 0;
    }

    //! @brief Exclude atoms from the same amino acid from being reported
    //! @param AAs aa vector
    //! @param SEQ_EXCLUSION minimum sequence separation between reported matches
    void VoxelGridAtom::ExcludeResiduesSameAA
    (
      const util::SiPtrVector< const biol::AABase> &AAS,
      const size_t &SEQ_EXCLUSION
    )
    {
      m_SeqExclusion = SEQ_EXCLUSION;
      size_t logical_aa_seq_id( 1);
      char prev_chain( '\0');
      int prev_seq_id( -1000);
      m_PdbIdToSeqID.Reset();
      m_PdbIdToSeqID.Resize( m_MaxPDBID + 1, 0);
      for( auto itr_aa( AAS.Begin()), itr_aa_end( AAS.End()); itr_aa != itr_aa_end; ++itr_aa, ++logical_aa_seq_id)
      {
        if( prev_chain != ( *itr_aa)->GetChainID())
        {
          prev_chain = ( *itr_aa)->GetChainID();
          logical_aa_seq_id += m_SeqExclusion + 1;
        }
        else if( prev_seq_id + 1 < ( *itr_aa)->GetSeqID())
        {
          logical_aa_seq_id += ( *itr_aa)->GetSeqID() - prev_seq_id - 1;
        }
        prev_seq_id = ( *itr_aa)->GetSeqID();
        for
        (
          auto itr_atom( ( *itr_aa)->GetAtoms().Begin()), itr_atom_end( ( *itr_aa)->GetAtoms().End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          if( ( *itr_atom)->GetPdbID() >= 0 && ( *itr_atom)->GetPdbID() < m_MaxPDBID)
          {
            m_PdbIdToSeqID( ( *itr_atom)->GetPdbID()) = logical_aa_seq_id;
          }
        }
      }
      m_PdbIdToGroupIdPtr = util::ToSiPtr( m_PdbIdToSeqID);
    }

    //! @brief Exclude residues from the same sse acid from being reported
    //! @note ExcludeResiduesSameAA should not be called if this function is used
    //! @param AAs aa vector
    void VoxelGridAtom::ExcludeResiduesSameSSE
    (
      const util::SiPtrVector< const SSE> &SSES,
      const util::SiPtrVector< const biol::AABase> &AAS,
      const size_t &SEQ_EXCLUSION
    )
    {
      ExcludeResiduesSameAA( AAS, SEQ_EXCLUSION);
      m_PdbIdToSSEID.Reset();
      m_PdbIdToSSEID.Resize( m_MaxPDBID + 1, util::GetUndefined< size_t>());
      m_NumberSSEs = SSES.GetSize();
      // hash sses for all amino acid ids
      size_t current_sse_id( 0);
      for( auto itr_sses( SSES.Begin()), itr_sse_end( SSES.End()); itr_sses != itr_sse_end; ++itr_sses, ++current_sse_id)
      {
        if( !( *itr_sses)->GetType()->IsStructured())
        {
          continue;
        }
        for( auto itr_aa( ( *itr_sses)->Begin()), itr_aa_end( ( *itr_sses)->End()); itr_aa != itr_aa_end; ++itr_aa)
        {
          for
          (
            auto itr_atom( ( *itr_aa)->GetAtoms().Begin()), itr_atom_end( ( *itr_aa)->GetAtoms().End());
            itr_atom != itr_atom_end;
            ++itr_atom
          )
          {
            if( ( *itr_atom)->GetPdbID() >= 0 && ( *itr_atom)->GetPdbID() < m_MaxPDBID)
            {
              m_PdbIdToSSEID( ( *itr_atom)->GetPdbID()) = current_sse_id;
            }
          }
        }
      }
      m_PdbIdToGroupIdPtr = util::ToSiPtr( m_PdbIdToSSEID);
    }

    namespace
    {
      //! Return true if the distance of A is less than B
      bool LessThanDistance
      (
        const storage::Triplet
        <
          util::SiPtr< const biol::Atom>,
          util::SiPtr< const biol::Atom>,
          double
        > &A,
        const storage::Triplet
        <
          util::SiPtr< const biol::Atom>,
          util::SiPtr< const biol::Atom>,
          double
        > &B
      )
      {
        return A.Third() < B.Third();
      }
    }
    //! @brief Get a matrix with counts of AAs interacting at particular distance between each SSE
    linal::Matrix< float> VoxelGridAtom::GetSSEInteractionMatrix( const double &RESOLUTION) const
    {
      BCL_Assert( m_PdbIdToGroupIdPtr.IsDefined() && &*m_PdbIdToGroupIdPtr == &m_PdbIdToSSEID, "SSEs not provided!");
      linal::Matrix< float> interactions( m_NumberSSEs, m_NumberSSEs, float( 0));
      auto interactions_vector( GetNeighbors( RESOLUTION));
      interactions_vector.Sort( LessThanDistance);
      storage::Vector< double> closest_distance( m_PdbIdToSSEID.GetSize(), util::GetUndefined< double>());
      storage::Vector< size_t> residue_counts( m_PdbIdToSSEID.GetSize(), size_t( 0));
      util::SiPtrVector
      <
        storage::Triplet
        <
          util::SiPtr< const biol::Atom>,
          util::SiPtr< const biol::Atom>,
          double
        >
      > true_closest_distances;
      for( auto itr( interactions_vector.Begin()), itr_end( interactions_vector.End()); itr != itr_end; ++itr)
      {
        if
        (
          m_PdbIdToSSEID( itr->First()->GetPdbID()) == util::GetUndefined< size_t>()
          || m_PdbIdToSSEID( itr->Second()->GetPdbID()) == util::GetUndefined< size_t>()
          || std::abs( int( m_PdbIdToSeqID( itr->First()->GetPdbID())) - int( m_PdbIdToSeqID( itr->Second()->GetPdbID()))) <= m_SeqExclusion
          || m_PdbIdToSSEID( itr->First()->GetPdbID()) == m_PdbIdToSSEID( itr->Second()->GetPdbID())
        )
        {
          continue;
        }
        if( !util::IsDefined( closest_distance( itr->First()->GetPdbID())))
        {
          closest_distance( itr->First()->GetPdbID()) = itr->Third();
        }
        if( !util::IsDefined( closest_distance( itr->Second()->GetPdbID())))
        {
          closest_distance( itr->Second()->GetPdbID()) = itr->Third();
        }
        if( closest_distance( itr->First()->GetPdbID()) + 1.0 < itr->Third())
        {
          continue;
        }
        if( closest_distance( itr->Second()->GetPdbID()) + 1.0 < itr->Third())
        {
          continue;
        }
        true_closest_distances.PushBack( *itr);
        ++residue_counts( itr->First()->GetPdbID());
        ++residue_counts( itr->Second()->GetPdbID());
      }
      for( auto itr( true_closest_distances.Begin()), itr_end( true_closest_distances.End()); itr != itr_end; ++itr)
      {
        const auto &triplet( **itr);
        const double interaction_weight
        (
          4.0
          /
          std::max( 4.0, double( residue_counts( triplet.First()->GetPdbID()) + residue_counts( triplet.Second()->GetPdbID())))
        );
        interactions( m_PdbIdToSSEID( triplet.First()->GetPdbID()), m_PdbIdToSSEID( triplet.Second()->GetPdbID())) += interaction_weight;
        interactions( m_PdbIdToSSEID( triplet.Second()->GetPdbID()), m_PdbIdToSSEID( triplet.First()->GetPdbID())) += interaction_weight;
      }
      return interactions;
    }

  } // namespace assemble
} // namespace bcl
