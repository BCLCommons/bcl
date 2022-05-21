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

#ifndef BCL_DENSITY_SIMULATORS_H_
#define BCL_DENSITY_SIMULATORS_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_density_simulate_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Simulators
    //! @brief this class enumerates implementations of SimulateInterface
    //! @details It provides convenience functions, to easily construct an Simulate object with the members set correctly
    //!
    //! @see @link example_density_simulators.cpp @endlink
    //! @author woetzen
    //! @date Jul 13, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Simulators :
      public util::Enumerate< util::ShPtr< SimulateInterface>, Simulators>
    {
      friend class util::Enumerate< util::ShPtr< SimulateInterface>, Simulators>;
    public:

    //////////
    // data //
    //////////

      const Simulator e_GaussianSphere;   //!< atom gaussian sphere, with atomic weight intensity
      const Simulator e_NoSmoothing;      //!< trilinear interpolation without smoothing
      const Simulator e_Gaussian;         //!< trilinear interpolation with gaussian smoothing
      const Simulator e_Triangular;       //!< trilinear interpolation with triangular smoothing
      const Simulator e_SemiEpanechnikov; //!< trilinear interpolation with semi epanechnikov smoothing
      const Simulator e_Epanechnikov;     //!< trilinear interpolation with epanechnikov smoothing
      const Simulator e_HardSphere;       //!< atom as hard sphere

      //! @brief default grid spacing
      static const linal::Vector3D &GetDefaultGridSpacing();

      //! @brief default resolution for density simulation
      static double GetDefaultResolution();

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Simulators();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a simulator from given grid spacing and resolution
      //! @param SIMULATOR the simluator enum
      //! @param GRID_SPACING desired spacing for density map
      //! @param RESOLUTION desired resolution for density map
      //! @return ShPtr to SimulatInterface for that simluator with grid spacing and resolution set
      util::ShPtr< SimulateInterface> CreateSimulator
      (
        const Simulator &SIMULATOR,
        const linal::Vector3D &GRID_SPACING,
        const double RESOLUTION
      ) const;

    }; // class Simulators

    //! @brief construct on access function for all Simulators
    //! @return reference to only instances of Simulators
    BCL_API
    Simulators &GetSimulators();

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< density::SimulateInterface>, density::Simulators>;

  } // namespace util
} // namespace bcl

#endif // BCL_DENSITY_SIMULATORS_H_ 
