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

#ifndef BCL_ASSEMBLE_VOXEL_GRID_MUTATION_H_
#define BCL_ASSEMBLE_VOXEL_GRID_MUTATION_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_complete.h"
#include "biol/bcl_biol_atom.h"
#include "biol/bcl_biol_mutation.h"
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
    //! @class VoxelGridMutation
    //! @brief Manages the Slicelists for aminoacid-objects
    //! @details Once constructed with a collection of objects which are derived from AABase, an instance of this
    //!          class provides the relevant aminoacids for one given aminoacid depending on their positions in 3D.
    //!          For this calculation the firstSideChainAtom's position is taken into account for all aminoacids
    //!          excepting Glycine, for which the Ha2-atom's position is used.
    //!
    //! @see @link example_assemble_voxel_grid_mutation.cpp @endlink
    //! @author mendenjl
    //! @date Feb 01, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API VoxelGridMutation :
      public util::VoxelGrid< biol::Mutation>
    {

    private:

    //////////
    // data //
    //////////

      bool m_PreferCA; //!< if true, use CA coordinates (rather than first sidechain atom) whenever they are available
      bool m_PreferCenter;
      mutable storage::List< linal::Vector3D> m_Centers;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param RESOLUTION resolution of the grid
      //! @param CACHE_EDGES whether to cache edges between adjacent voxels. Useful when reusing the same grid repeatedly
      VoxelGridMutation( const double &RESOLUTION = 4.0, const bool &CACHE_EDGES = false, const bool &PREFER_CA = false, const bool &PREFER_CENTER = true) :
        util::VoxelGrid< biol::Mutation>( RESOLUTION, CACHE_EDGES),
        m_PreferCA( PREFER_CA),
        m_PreferCenter( PREFER_CENTER)
      {
      }

    //////////
    // data //
    //////////

      //! @brief copy constructor
      //! @return pointer to a copy of the actual object
      VoxelGridMutation *Clone() const
      {
        return new VoxelGridMutation( *this);
      }

      //! @brief GetClassIdentifier returns class name of the object
      //! @return returns string with the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief Update the VoxelGrid with new data
      //! @param NEW_DATA SiPtrVector of the new data
      void SetObjects( const util::SiPtrVector< const biol::Mutation> &NEW_DATA)
      {
        m_Centers.Reset();
        util::VoxelGrid< biol::Mutation>::SetObjects( NEW_DATA);
      }

      //! @brief Extract an AA-instance's coordinates (Ha2 for Glycine, Cb-Atom otherwise)
      //! @param AA pointer to the amino acid object
      //! @return Reference to the coordinates. null-ptr if coordinates are undefined
      util::SiPtr< const linal::Vector3D> ExtractPosition( const util::SiPtr< const biol::Mutation> &MUTATION) const
      {
        BCL_Assert( MUTATION->GetChainIDs().size() == size_t( 1), "Cannot extract one position from multi-chain mutation");
        return ExtractPositionFromChain( MUTATION, size_t( 0));
      }

      //! @brief extract the 3D coordinates of a given t_DataType input. TO BE IMPLEMENTED BY DERIVED CLASSES if one input can have multiple positions
      //! @param INPUT pointer to the t_DataType
      //! @return 3D Vector of the input's coordinates
      util::SiPtrVector< const linal::Vector3D> ExtractPositions( const util::SiPtr< const biol::Mutation> &INPUT) const;

      //! @brief check if two AA objects are the same, their equalness being defined by SequenceID and ChainID
      //! @param ITEM_1 first AA to compare
      //! @param ITEM_2 second AA to compare
      //! @return true if AAs are the same, false otherwise
      bool IsSameItem( const biol::Mutation &ITEM_1, const biol::Mutation &ITEM_2) const
      {
        return ITEM_1 == ITEM_2;
      }

    private:

      //! @brief Extract an AA-instance's coordinates (Ha2 for Glycine, Cb-Atom otherwise)
      //! @param AA pointer to the amino acid object
      //! @return Reference to the coordinates. null-ptr if coordinates are undefined
      util::SiPtr< const linal::Vector3D> ExtractPositionFromChain
      (
        const util::SiPtr< const biol::Mutation> &MUTATION,
        const size_t &CHAIN_INDEX
      ) const;

    }; //class VoxelGridMutation

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_VOXEL_GRID_MUTATION_H_
