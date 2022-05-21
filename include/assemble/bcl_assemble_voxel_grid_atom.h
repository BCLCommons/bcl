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

#ifndef BCL_ASSEMBLE_VOXEL_GRID_ATOM_H_
#define BCL_ASSEMBLE_VOXEL_GRID_ATOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_matrix.h"
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
    //! @class VoxelGridAtom
    //! @brief Manages the Slicelists for aminoacid-objects
    //! @details Once constructed with a collection of objects which are derived from AABase, an instance of this
    //!          class provides the relevant aminoacids for one given aminoacid depending on their positions in 3D.
    //!          For this calculation the firstSideChainAtom's position is taken into account for all aminoacids
    //!          excepting Glycine, for which the Ha2-atom's position is used.
    //!
    //! @see @link example_assemble_voxel_grid_atom.cpp @endlink
    //! @author mendenjl
    //! @date Mar 10, 2017
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API VoxelGridAtom :
      public util::VoxelGrid< biol::Atom>
    {

    private:

      storage::Vector< size_t> m_PdbIdToSSEID;    //!< PDB - SSE ID mapping. if this vector is provided
      storage::Vector< size_t> m_PdbIdToSeqID;    //!< PDB - AA ID mapping
      size_t                   m_SeqExclusion;    //!< Minimum sequence separation
      mutable size_t           m_MaxPDBID;        //!< Maximum pdb id observed
      size_t                   m_NumberSSEs;      //!< Number of atoms interacting between SSEs

      //<! Points to the appropriate vector to use for excluding contacts
      util::SiPtr< const storage::Vector< size_t> > m_PdbIdToGroupIdPtr;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param RESOLUTION resolution of the grid
      //! @param SEQ_EXCLUSION minimum sequence distance between residues (<=)
      //! @param CACHE_EDGES whether to cache edges between adjacent voxels. Useful when reusing the same grid repeatedly
      VoxelGridAtom
      (
        const double &RESOLUTION = 4.0,
        const bool &CACHE_EDGES = false
      ) :
        util::VoxelGrid< biol::Atom>( RESOLUTION, CACHE_EDGES),
        m_SeqExclusion( 0),
        m_MaxPDBID( 0)
      {
      }

    //////////
    // data //
    //////////

      //! @brief copy constructor
      //! @return pointer to a copy of the actual object
      VoxelGridAtom *Clone() const
      {
        return new VoxelGridAtom( *this);
      }

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns string with the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief set the group ids for the atoms (indexed by atom PDB_ID). Atoms in the same group will not be reported
      //!        as in-contact.
      //! @param GROUP_IDS a vector containing group ids for all atoms of interest.
      void SetGroupIDs( const storage::Vector< size_t> &GROUP_IDS);

      //! @brief Exclude atoms from the same amino acid from being reported
      //! @param AAs aa vector
      //! @param SEQ_EXCLUSION minimum sequence separation between reported matches
      void ExcludeResiduesSameAA
      (
        const util::SiPtrVector< const biol::AABase> &AAS,
        const size_t &SEQ_EXCLUSION = 0
      );

      //! @brief Exclude residues from the same sse acid from being reported
      //! @note ExcludeResiduesSameAA should not be called if this function is used
      //! @param AAs aa vector
      void ExcludeResiduesSameSSE
      (
        const util::SiPtrVector< const SSE> &SSES,
        const util::SiPtrVector< const biol::AABase> &AAS,
        const size_t &SEQ_EXCLUSION = 0
      );

      //! @brief Extract an atoms coordinates (Ha2 for Glycine, Cb-Atom otherwise)
      //! @param ATOM pointer to the atom object
      //! @return Reference to the coordinates. null-ptr if coordinates are undefined
      util::SiPtr< const linal::Vector3D> ExtractPosition( const util::SiPtr< const biol::Atom> &ATOM) const;

      //! @brief Update the VoxelGrid with new data
      //! @param NEW_DATA SiPtrVector of the new data
      void SetObjects( const util::SiPtrVector< const biol::Atom> &NEW_DATA);

      //! @brief check if two AA objects are the same, their equalness being defined by SequenceID and ChainID
      //! @param ITEM_1 first AA to compare
      //! @param ITEM_2 second AA to compare
      //! @return true if AAs are the same, false otherwise
      bool IsSameItem( const biol::Atom &ITEM_1, const biol::Atom &ITEM_2) const
      {
        return ITEM_1.GetPdbID() == ITEM_2.GetPdbID();
      }

      //! @brief Get a matrix with counts of AAs interacting at particular distance between each SSE
      linal::Matrix< float> GetSSEInteractionMatrix( const double &RESOLUTION) const;

    }; //class VoxelGridAtom

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_VOXEL_GRID_ATOM_H_
