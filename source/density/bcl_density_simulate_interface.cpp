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
#include "density/bcl_density_simulate_interface.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_mask_3d.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief measure the maximal map extend
    //! @param ATOMS list of atoms
    //! @return min coord xyz and max coord xyz
    storage::VectorND< 2, linal::Vector3D> SimulateInterface::DetermineGridCorners
    (
      const util::SiPtrVector< const biol::Atom> &ATOMS
    )
    {
      util::SiPtrVector< const linal::Vector3D> coordinates;
      coordinates.AllocateMemory( ATOMS.GetSize());

      // iterate over all atoms to acquire coordinates
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator atom_itr( ATOMS.Begin()), atom_itr_end( ATOMS.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        coordinates.PushBack( util::ToSiPtr( ( *atom_itr)->GetCoordinates()));
      }

      // determine grid corners for the coordinates
      return Mask3d::DetermineGridCorners( coordinates);
    }

  } // namespace density
} // namespace bcl
