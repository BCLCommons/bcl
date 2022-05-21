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
#include "coord/bcl_coord_geometric_hash_storage_classes.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "coord/bcl_coord_geometric_hash_storage_hash_map.h"
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
    GeometricHashStorageClasses::GeometricHashStorageClasses() :
      e_HashMap( AddEnum( "HashMap", util::ShPtr< GeometricHashStorageInterface>( new GeometricHashStorageHashMap())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GeometricHashStorageClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief flag for choosing the geometris hash storage over the commandline
    const util::ShPtr< command::FlagInterface> &
    GeometricHashStorageClasses::GetGeometricHashStorageClassFlag() const
    {
      static const util::ShPtr< command::FlagInterface> s_geometric_hash_storage_class_flag
      (
        new command::FlagStatic
        (
          "hash_storage",
          "choice of storage class",
          command::Parameter
          (
            "hash_storage_class",
            "name for one of the available hash storage classes - depending on size and speed requirements",
            command::ParameterCheckEnumerate< GeometricHashStorageClasses>(),
            e_HashMap
          )
        )
      );

      // end
      return s_geometric_hash_storage_class_flag;
    }

    //! @brief construct a geometric hash storage class form the commandline argument
    util::ShPtr< GeometricHashStorageInterface>
    GeometricHashStorageClasses::ConstructFromCommandline
    (
      const std::string &MRC_NAME,
      const double MRC_RESOLUTION,
      const double MRC_VOXELSIZE,
      const storage::VectorND< 3, size_t> &MRC_EXTENSION,
      const size_t NUMBER_POINTS,
      const double FEATURE_DISTANCE,
      const double RATIO_INTENSITY_GRADIENT,
      const storage::VectorND< 4, double> &THRESHOLD,
      const double FEATURE_RADIUS,
      const PointToKeyInterface &POINT_TO_KEY,
      const size_t MIN_NUMBER_NEIGHBORS,
      const double MIN_NEIGHBOR_DISTANCE
    ) const
    {
      return util::ShPtr< GeometricHashStorageInterface>
        (
          (
            *GeometricHashStorageClass( GetGeometricHashStorageClassFlag()->GetFirstParameter()->GetValue())
          )->Construct
          (
            MRC_NAME,
            MRC_RESOLUTION,
            MRC_VOXELSIZE,
            MRC_EXTENSION,
            NUMBER_POINTS,
            FEATURE_DISTANCE,
            RATIO_INTENSITY_GRADIENT,
            THRESHOLD,
            FEATURE_RADIUS,
            POINT_TO_KEY,
            MIN_NUMBER_NEIGHBORS,
            MIN_NEIGHBOR_DISTANCE
          )
        );
    }

    //! @brief construct on access function for all GeometricHashStorageClasses
    //! @return reference to only instances of GeometricHashStorageClasses
    GeometricHashStorageClasses &GetGeometricHashStorageClasses()
    {
      return GeometricHashStorageClasses::GetEnums();
    }

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< coord::GeometricHashStorageInterface>, coord::GeometricHashStorageClasses>;

  } // namespace util
} // namespace bcl
