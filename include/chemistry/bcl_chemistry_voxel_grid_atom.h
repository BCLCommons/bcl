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

#ifndef BCL_CHEMISTRY_VOXEL_GRID_ATOM_H_
#define BCL_CHEMISTRY_VOXEL_GRID_ATOM_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_voxel_grid.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
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
    //! @see @link example_chemistry_voxel_grid_atom.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Jan 15, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API VoxelGridAtom :
      public util::VoxelGrid< AtomConformationalInterface>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param RESOLUTION resolution of the grid
      //! @param CACHE_EDGES whether to cache edges between adjacent voxels. Useful when reusing the same grid repeatedly
      VoxelGridAtom( const double &RESOLUTION = 4.0, const bool &CACHE_EDGES = false) :
        util::VoxelGrid< AtomConformationalInterface>( RESOLUTION, CACHE_EDGES)
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

      //! @brief Extract an AA-instance's coordinates (Ha2 for Glycine, Cb-Atom otherwise)
      //! @param AA pointer to the amino acid object
      //! @return Reference to the coordinates. null-ptr if coordinates are undefined
      util::SiPtr< const linal::Vector3D> ExtractPosition( const util::SiPtr< const AtomConformationalInterface> &ATOM) const
      {
        if( !ATOM->GetPosition().IsDefined())
        {
          return util::SiPtr< const linal::Vector3D>();
        }
        return util::SiPtr< const linal::Vector3D>( ATOM->GetPosition());
      }

      //! @brief check if two AA objects are the same, their equalness being defined by SequenceID and ChainID
      //! @param ITEM_1 first AA to compare
      //! @param ITEM_2 second AA to compare
      //! @return true if AAs are the same, false otherwise
      bool IsSameItem( const AtomConformationalInterface &ITEM_1, const AtomConformationalInterface &ITEM_2) const
      {
        return &ITEM_1 == &ITEM_2;
      }

    }; //class VoxelGridAtom

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_VOXEL_GRID_ATOM_H_
