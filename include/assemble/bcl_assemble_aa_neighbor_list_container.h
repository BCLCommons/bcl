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

#ifndef BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_H_
#define BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_aa_neighbor_list.h"
#include "biol/bcl_biol_aa_compare.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborListContainer
    //! @brief Class for storing a list of AANeighborList instances to represent all amino acid distances in a protein
    //! @details This class provides storage for AANeighborList instances for each amino acid in a given vector of
    //! amino acids. It uses a set to ensure there is only one AANeighborList for each amino acid. The constructor
    //! initializes the lists according to the given distance cutoff (distances above the cutoff are omitted) and
    //! the minimal sequence separation so that neighbor residues on the sequence are not on each others' lists.
    //!
    //! @see @link example_assemble_aa_neighbor_list_container.cpp @endlink
    //! @author woetzen
    //! @date Jun 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborListContainer :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      typedef storage::Map< util::SiPtr< const biol::AABase>, AANeighborList, biol::AALessThanSeqID> NeighborListsType;

      //! set of neighbor list
      NeighborListsType m_NeighborLists;

      //! distance cutoff above which no neighbors are stored
      double m_DistanceCutoff;

      //! minimal sequence separation below which, no neighbors are stored
      size_t m_MinimalSequenceSeparation;

      //! bool indicating whether or not amino acids of different chains are considered neighbors
      bool m_ConsiderDifferentChain;

      //! @brief iterator on container of neighbor list
      typedef NeighborListsType::iterator iterator;

    public:

      //! @brief iterator on container of neighbor list
      typedef NeighborListsType::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from distance cutoff and minimal sequence separation
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborListContainer
      (
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief construct from a list of amino acids, distance cutoff and minimal sequence separation
      //! @param AMINO_ACIDS the list of amino acids in question
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborListContainer
      (
        const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief construct from two lists of amino acids, distance cutoff and minimal sequence separation
      //! @details only inter pairwise distances and neighbors are calculated
      //! @param AMINO_ACIDS_A, AMINO_ACIDS_B the list of amino acids in question
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborListContainer
      (
        const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS_A,
        const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS_B,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief Copy constructor
      //! @param AMINO_ACIDS the list of amino acids in question
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborListContainer
      (
        const AANeighborListContainer &AMINO_ACIDS,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief construct from already-found amino acids (e.g. using assemble::VoxelGridAA.GetSSEConnections())
      //! @param NEIGHBORS List of neighbors
      AANeighborListContainer
      (
        const storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> > &NEIGHBORS,
        const util::SiPtrVector< const biol::AABase> &AAS,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief Clone function
      //! @return pointer to new AANeighborListContainer
      AANeighborListContainer *Clone() const
      {
        return new AANeighborListContainer( *this);
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

      //! @brief get the size
      //! @return number of center amino acids
      size_t GetSize() const
      {
        return m_NeighborLists.GetSize();
      }

      //! @brief return whether is empty
      //! @return whether is empty
      bool IsEmpty() const
      {
        return m_NeighborLists.IsEmpty();
      }

      //! @brief get the total number of neighbors
      //! @return total number of neighbors
      size_t GetNumberNeighbors() const;

      //! @brief access to distance cutoff
      //! @return the distance above which no neighbors are stored
      double GetDistanceCutoff() const
      {
        return m_DistanceCutoff;
      }

      //! @brief access to minimal sequence separation
      //! @return minimal sequence separation, under which no neighbors are stored, if amino acid is in the same chain
      size_t GetMinimalSequenceSeparation() const
      {
        return m_MinimalSequenceSeparation;
      }

      //! @brief are different chain amino acids considered neighbors
      //! @return bool indicating whether or not amino acids of different chains are considered neighbors
      bool GetConsiderDifferentChain() const
      {
        return m_ConsiderDifferentChain;
      }

      //! @brief iterator on first neighbor
      //! @return const_iterator on first neighbor with distance
      const_iterator Begin() const
      {
        return m_NeighborLists.Begin();
      }

      //! @brief iterator on end neighbor
      //! @return const_iterator of the end neighbor
      const_iterator End() const
      {
        return m_NeighborLists.End();
      }

      //! @brief find a neighbor iterator
      //! @param AMINO_ACID the amino acid in consideration
      //! @return const_iterator to amino acid and its distance, end if not found
      const_iterator Find( const biol::AABase &AMINO_ACID) const
      {
        return m_NeighborLists.Find( util::ToSiPtr( AMINO_ACID));
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief prune the AANeighborListContainer
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      void Prune
      (
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief Remove residues located on the same SSE
      //! @param MODEL the protein model containing the SSEs
      void PruneResiduesSameSSE( const ProteinModel &MODEL);

      //! @brief return the intersect size with the given AANeighborListContainer
      //! @param CONTAINER The container with which the intersect is going to be calculated
      //! @return number of overlapping edges
      size_t IntersectionSize( const AANeighborListContainer &CONTAINER) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AANeighborListContainer

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_AA_NEIGHBOR_LIST_CONTAINER_H_
