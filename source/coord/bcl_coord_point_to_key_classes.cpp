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
#include "coord/bcl_coord_point_to_key_classes.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_point_to_key_spherical_radius.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PointToKeyClasses::PointToKeyClasses() :
      e_SphericalRadius( AddEnum( PointToKeySphericalRadius::GetCoordinateSystemDescriptor(), util::ShPtr< PointToKeyInterface>( new PointToKeySphericalRadius( PointToKeyInterface::GetDefaultAngularResolution()))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeyClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get point to key function for given name
    //! @param NAME
    //! @param ANGULAR_RESOLUTION
    //! @param DISTANCE_RESOLUTION
    //! @return ShPtr to PointToKey function
    util::ShPtr< PointToKeyInterface>
    PointToKeyClasses::GetPointToKeyFunction
    (
      const std::string &NAME,
      const double ANGULAR_RESOLUTION,
      const double DISTANCE_RESOLUTION
    ) const
    {
      // get enum with that name and make clone of point to key function
      util::ShPtr< PointToKeyInterface> function( ( *EnumType( NAME))->Clone());

      // set resolution
      function->SetAngularResolution( ANGULAR_RESOLUTION);
      function->SetDistanceResolution( DISTANCE_RESOLUTION);

      // end
      return function;
    }

    PointToKeyClasses &GetPointToKeyClasses()
    {
      return PointToKeyClasses::GetEnums();
    }

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< coord::PointToKeyInterface>, coord::PointToKeyClasses>;

  } // namespace util
} // namespace bcl
