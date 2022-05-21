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
#include "density/bcl_density_simulators.h"

// includes from bcl - sorted alphabetically
#include "density/bcl_density_simulate_default.h"
#include "density/bcl_density_simulate_gaussian_sphere.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    //! @brief default grid spacing
    const linal::Vector3D &Simulators::GetDefaultGridSpacing()
    {
      static const linal::Vector3D s_default_grid_spacing( 2.2, 2.2, 2.2);
      return s_default_grid_spacing;
    }

    //! @brief default resolution for density simulation
    double Simulators::GetDefaultResolution()
    {
      return 6.6;
    }

    //! @brief default constructor
    Simulators::Simulators() :
      e_GaussianSphere  ( AddEnum( "GaussianSphere"                        , util::ShPtr< SimulateInterface>( new SimulateGaussianSphere( GetDefaultGridSpacing(), GetDefaultResolution())))),
      e_NoSmoothing     ( AddEnum( "TrilinearInterpolation"                , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_HardSphere)))),
      e_Gaussian        ( AddEnum( "TrilinearInterpolationGaussian"        , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_Gaussian)))),
      e_Triangular      ( AddEnum( "TrilinearInterpolationTriangular"      , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_Triangular)))),
      e_SemiEpanechnikov( AddEnum( "TrilinearInterpolationSemiEpanechnikov", util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_SemiEpanechnikov)))),
      e_Epanechnikov    ( AddEnum( "TrilinearInterpolationEpanechnikov"    , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_Epanechnikov)))),
      e_HardSphere      ( AddEnum( "TrilinearInterpolationHardSphere"      , util::ShPtr< SimulateInterface>( new SimulateDefault( GetDefaultGridSpacing(), GetDefaultResolution(), SimulateDefault::e_HardSphere))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Simulators::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct a simulator from given grid spacing and resolution
    //! @param SIMULATOR the simluator enum
    //! @param RESOLUTION desired resolution for density map
    //! @brief GRID_SPACING desired spacing for density map
    //! @return ShPtr to SimulatInterface for that simluator with grid spacing and resolution set
    util::ShPtr< SimulateInterface> Simulators::CreateSimulator
    (
      const Simulator &SIMULATOR,
      const linal::Vector3D &GRID_SPACING,
      const double RESOLUTION
    ) const
    {
      // copy the SimulateInterface derived class for that enumerator
      util::ShPtr< SimulateInterface> simulator( SIMULATOR->HardCopy());

      // set the grid spacing and resolution
      simulator->SetGridSpacing( GRID_SPACING);
      simulator->SetResolution( RESOLUTION);

      // end
      return simulator;
    }

    //! @brief construct on access function for all Simulators
    //! @return reference to only instances of Simulators
    Simulators &GetSimulators()
    {
      return Simulators::GetEnums();
    }

  } // namespace density

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< density::SimulateInterface>, density::Simulators>;

  } // namespace util
} // namespace bcl
