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

#ifndef BCL_DENSITY_SIMULATE_INTERFACE_H_
#define BCL_DENSITY_SIMULATE_INTERFACE_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SimulateInterface
    //! @brief an interface class for simulating density maps from Atoms
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Jun 21, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SimulateInterface :
      public math::FunctionInterfaceSerializable< util::SiPtrVector< const biol::Atom>, Map>
    {

    private:

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SimulateInterface
      virtual SimulateInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief set the resolution
      //! @param RESOLUTION the resolution for the density map to be generated
      virtual void SetResolution( const double RESOLUTION) = 0;

      //! @brief set the resolution
      virtual double GetResolution() const = 0;

      //! @brief set the grid spacing
      //! @param GRID_SPACING the width of a grid element in x, y and z
      virtual void SetGridSpacing( const linal::Vector3D &GRID_SPACING) = 0;

      //! @brief set the margin
      //! @param MARGIN number of additional cells next to last atom occupied cells
      virtual void SetMargin( const size_t MARGIN) = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief generate simulated density from given list of atoms
      //! @param ATOMS siptrvector of atoms
      //! @return a simulated density map
      virtual Map operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief measure the maximal map extend
      //! @param ATOMS list of atoms
      //! @return min coord xyz and max coord xyz
      static storage::VectorND< 2, linal::Vector3D> DetermineGridCorners
      (
        const util::SiPtrVector< const biol::Atom> &ATOMS
      );

    }; // class SimulateInterface

  } // namespace density
} // namespace bcl

#endif // BCL_DENSITY_SIMULATE_INTERFACE_H_ 
