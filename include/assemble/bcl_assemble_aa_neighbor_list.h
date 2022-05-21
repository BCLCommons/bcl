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

#ifndef BCL_ASSEMBLE_AA_NEIGHBOR_LIST_H_
#define BCL_ASSEMBLE_AA_NEIGHBOR_LIST_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_compare.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AANeighborList
    //! @brief Class for storing for an amino acid, distances to all other amino acids
    //! @details This is a storage class for storing amino acid - amino acid distances. An instance of this class stores
    //! for a specific amino acid, distances to all other amino acids. The class also contains information as to the
    //! distance cutoff used, ( distances above the cutoff are not stored in the member map), as well as the minimal
    //! sequence separation.
    //!
    //! @see @link example_assemble_aa_neighbor_list.cpp @endlink
    //! @author woetzen
    //! @date Jun 26, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AANeighborList :
      public util::ObjectInterface
    {
    /////////////
    // friends //
    /////////////

      //! AANeighborListContainer so that it can access constructor for easier finding in set
      friend class AANeighborListContainer;

    public:

      //! Typedef for internal data structure type
      typedef storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> > NeighborContainerType;

    private:

    //////////
    // data //
    //////////

      //! pointer to center amino acid, for which this neighbor list is generated for
      util::SiPtr< const biol::AABase> m_CenterAminoAcid;

      //! map of neighbor amino acids, with their distance
      NeighborContainerType m_Neighbors;

      //! distance cutoff above which no neighbors are stored
      double m_DistanceCutoff;

      //! minimal sequence separation below which, no neighbors are stored
      size_t m_MinimalSequenceSeparation;

      //! bool indicating whether or not amino acids of different chains are considered neighbors
      bool m_ConsiderDifferentChain;

    public:

      //! @brief iterator on neighbor map
      typedef storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> >::const_iterator const_iterator;
      typedef storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> >::iterator       iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief access to default center aa for undefined list
      static const biol::AABase &GetDefaultCenterAA();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief helper constructor to find a neighbor list in a AANeihborListCotainer from given amino acid
      //! @brief CENTER_AMINO_ACID amino acid for which the neighbor list is generated
      AANeighborList( const biol::AABase &CENTER_AMINO_ACID = GetDefaultCenterAA());

      //! @brief constructor from central AA, neighboring AAs and their distances
      //! @param CENTER_AMINO_ACID amino acid for which the neighbor list is generated
      //! @param AMINO_ACIDS_DISTANCES the list of amino acids in question along w/ their distances
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborList
      (
        const biol::AABase &CENTER_AMINO_ACID,
        const NeighborContainerType &AMINO_ACIDS_DISTANCES,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief construct from center amino acids, a list of amino acids, distance cutoff and minimal sequence separation
      //! @param CENTER_AMINO_ACID amino acid for which the neighbor list is generated
      //! @param AMINO_ACIDS the list of amino acids in question
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborList
      (
        const biol::AABase &CENTER_AMINO_ACID,
        const util::SiPtrVector< const biol::AABase> &AMINO_ACIDS,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief construct from a neighbor list setup with parameters that'd make it a superset of the desired list
      //! @param PARENT the original neighborlist
      //! @param DISTANCE_CUTOFF above which neighbors are not stored
      //! @param MIN_SEQ_SEPARATION minimal sequence distance (if within same chain) under which neighbors are not stored
      //! @param CONSIDER_DIFFERENT_CHAIN bool indicating whether or not amino acids of different chains are considered neighbors
      AANeighborList
      (
        const AANeighborList &PARENT,
        const double DISTANCE_CUTOFF,
        const size_t MIN_SEQ_SEPARATION,
        const bool CONSIDER_DIFFERENT_CHAIN
      );

      //! @brief Clone function
      //! @return pointer to new AANeighborList
      AANeighborList *Clone() const
      {
        return new AANeighborList( *this);
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

      //! @brief access to the center amino acid
      //! @return SiPtr to the center amino acid, that neighbor list is for
      const util::SiPtr< const biol::AABase> &GetCenterAminoAcid() const
      {
        return m_CenterAminoAcid;
      }

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

      //! @brief iterator on first neighbor
      //! @return const_iterator on first neighbor with distance
      const_iterator Begin() const
      {
        return m_Neighbors.Begin();
      }

      //! @brief iterator on end neighbor
      //! @return const_iterator of the end neighbor
      const_iterator End() const
      {
        return m_Neighbors.End();
      }

      //! @brief return the size
      //! @return size
      size_t GetSize() const
      {
        return m_Neighbors.GetSize();
      }

      //! @brief return whether the list is empty
      //! @return whether the list is empty
      bool IsEmpty() const
      {
        return m_Neighbors.IsEmpty();
      }

    ////////////////
    // operations //
    ////////////////

    private:

      //! @brief push back another amino acid into the neighbor. Checks that the added AA/dist are appropriate for this list
      //!        should be done prior to this
      void PushBack( const util::SiPtr< const biol::AABase> &AA, const double &DIST);

      //! @brief return the intersect with the given neighbor list, the distances will indicate the differences from
      //! this list to the given list for the common residues
      //! @param AA_NEIGHBOR_LIST the neighbor list to calculate the intersection with
      //! @return Number of elements in the intersection of the lists
      size_t IntersectionSize( const AANeighborList &AA_NEIGHBOR_LIST) const;

    public:

      //! @brief Remove neighbors with different chain id
      void RemoveNeighborsWithDifferentChainID();

      //! @brief decrease the distance cutoff
      //! @param DISTANCE_CUTOFF
      //! @return true, on success, false if DISTANCE_CUTOFF is larger than current distance cutoff
      bool DecreaseDistanceCutoff( const double DISTANCE_CUTOFF);

      //! @brief increase minimal sequence separation
      //! @param MIN_SEQ_SEPARATION
      //! @return true, if increment was successful, flase if the new seq separtion is smaller than the current
      bool IncreaseMinimalSequenceSeparation( const size_t MIN_SEQ_SEPARATION);

      //! @brief Remove residues located on the same SSE
      //! @param MODEL the protein model containing the SSEs
      void PruneResiduesSameSSE( const ProteinModel &MODEL);

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

    }; // class AANeighborList

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_AA_NEIGHBOR_LIST_H_
