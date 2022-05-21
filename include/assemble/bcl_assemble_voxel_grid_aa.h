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

#ifndef BCL_ASSEMBLE_VOXEL_GRID_AA_H_
#define BCL_ASSEMBLE_VOXEL_GRID_AA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_complete.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_voxel_grid.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class VoxelGridAA
    //! @brief Manages the Slicelists for aminoacid-objects
    //! @details Once constructed with a collection of objects which are derived from AABase, an instance of this
    //!          class provides the relevant aminoacids for one given aminoacid depending on their positions in 3D.
    //!          For this calculation the firstSideChainAtom's position is taken into account for all aminoacids
    //!          excepting Glycine, for which the Ha2-atom's position is used.
    //!
    //! @see @link example_assemble_voxel_grid_aa.cpp @endlink
    //! @author mendenjl
    //! @date Nov 22, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API VoxelGridAA :
      public util::VoxelGrid< biol::AABase>
    {

    private:

    //////////
    // data //
    //////////

      bool m_PreferCA; //!< if true, use CA coordinates (rather than first sidechain atom) whenever they are available

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param RESOLUTION resolution of the grid
      //! @param CACHE_EDGES whether to cache edges between adjacent voxels. Useful when reusing the same grid repeatedly
      VoxelGridAA( const double &RESOLUTION = 4.0, const bool &CACHE_EDGES = false, const bool &PREFER_CA = false) :
        util::VoxelGrid< biol::AABase>( RESOLUTION, CACHE_EDGES),
        m_PreferCA( PREFER_CA)
      {
      }

    //////////
    // data //
    //////////

      //! @brief copy constructor
      //! @return pointer to a copy of the actual object
      VoxelGridAA *Clone() const
      {
        return new VoxelGridAA( *this);
      }

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns string with the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief Extract an AA-instance's coordinates (Ha2 for Glycine, Cb-Atom otherwise)
      //! @param AA pointer to the amino acid object
      //! @return Reference to the coordinates. null-ptr if coordinates are undefined
      util::SiPtr< const linal::Vector3D> ExtractPosition( const util::SiPtr< const biol::AABase> &AA) const
      {
        // get firstSideChainAtom
        const biol::Atom &atom( m_PreferCA ? AA->GetCA() : AA->GetFirstSidechainAtom());

        // Check if it's coordinates are defined
        if( atom.AllCoordinatesDefined())
        { // if so return them
          return util::SiPtr< const linal::Vector3D>( atom.GetCoordinates());
        }
        return util::SiPtr< const linal::Vector3D>();
      }

      //! @brief check if two AA objects are the same, their equalness being defined by SequenceID and ChainID
      //! @param ITEM_1 first AA to compare
      //! @param ITEM_2 second AA to compare
      //! @return true if AAs are the same, false otherwise
      bool IsSameItem( const biol::AABase &ITEM_1, const biol::AABase &ITEM_2) const
      {
        return ( ITEM_1.GetSeqID() == ITEM_2.GetSeqID()) && ( ITEM_1.GetChainID() == ITEM_2.GetChainID());
      }

      //! @brief calculates AA pair clash score for a vector of AAs via Slicelist usage
      //! @note NOTE: This might have to be moved into a scoring class derived from protein model
      //! @param AAS vector of pointers to the AAs
      //! @return score of the AA pair clash
      double GetAAClashScore( const util::SiPtrVector< const biol::AABase> &AAS);

      //! @brief calculates AA pair clash score for a vector of AAs via Slicelist usage
      //! @note NOTE: This might have to be moved into a scoring class derived from protein model
      //! @param AAS vector of pointers to the AAs
      //! @return score of the AA pair clash
      double GetAAClashScore();

      //! @brief Get a matrix with counts of AAs interacting at particular distance between each SSE
      linal::Matrix< float> GetSSEInteractionMatrix
      (
        const util::SiPtrVector< const SSE> &SSES,
        const util::SiPtrVector< const biol::AABase> &AAS,
        const size_t &SEQ_EXCLUSION,
        const double &RESOLUTION,
        const bool   &CONSIDER_LOOPS,
        const double &MIN_CONTACT_P,
        const bool   &CONSIDER_POINT_CONTACTS,
        const bool   &MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY = false
      );

      //! @brief Get the SSE connections
      //! @param MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY if set, only return connections between SSEs that are mutually
      //!        closer than any other connection with each of the residues involved to any external sse
      storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
      GetSSEConnections
      (
        const util::SiPtrVector< const SSE> &SSES,
        const util::SiPtrVector< const biol::AABase> &AAS,
        const size_t &SEQ_EXCLUSION,
        const double &RESOLUTION,
        const bool   &CONSIDER_LOOPS,
        const bool   &MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY = false,
        const bool   &CONSIDER_INTRA_LOOP_CLASHES = true
      );

      //! @brief Get the SSE connections
      size_t GetMinSSEMovesToRemoveClashes
      (
        const util::SiPtrVector< const SSE> &SSES,
        const util::SiPtrVector< const biol::AABase> &AAS,
        const bool   &CONSIDER_LOOPS
      );

      //! @brief Get the SSE connections
      storage::Vector< storage::Pair< util::SiPtr< const SSE>, linal::Vector3D> >
      GetMinSSEMoveIDsToRemoveClashes
      (
        const util::SiPtrVector< const SSE> &SSES,
        const util::SiPtrVector< const biol::AABase> &AAS,
        const bool   &CONSIDER_LOOPS
      );

    }; //class VoxelGridAA

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_VOXEL_GRID_AA_H_
