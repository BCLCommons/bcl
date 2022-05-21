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

#ifndef BCL_COORD_GEOMETRIC_HASH_STORAGE_CLASSES_H_
#define BCL_COORD_GEOMETRIC_HASH_STORAGE_CLASSES_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GeometricHashStorageClasses
    //! @brief enumerates implementations of hash storage classes
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_coord_geometric_hash_storage_classes.cpp @endlink
    //! @author woetzen
    //! @date May 29, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GeometricHashStorageClasses :
      public util::Enumerate< util::ShPtr< GeometricHashStorageInterface>, GeometricHashStorageClasses>
    {
      friend class util::Enumerate< util::ShPtr< GeometricHashStorageInterface>, GeometricHashStorageClasses>;

    private:

    //////////
    // data //
    //////////

    public:

      //! hashmap storage
      const GeometricHashStorageClass e_HashMap;

    //////////
    // data //
    //////////

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GeometricHashStorageClasses();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief flag for choosing the geometris hash storage over the commandline
      const util::ShPtr< command::FlagInterface> &GetGeometricHashStorageClassFlag() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a geometric hash storage class form the commandline argument
      util::ShPtr< GeometricHashStorageInterface> ConstructFromCommandline
      (
        const std::string &MRC_NAME,
        const double MRC_RESOLUTION,
        const double MRC_VOXELSIZE,
        const storage::VectorND< 3, size_t> &MRC_EXTENSION,
        const size_t NUMBER_POINTS,
        const double FEATURE_DISTANCE,
        const double RATIO_INTENSITY_GRADIENT,
        const storage::VectorND< 4, double> &THRESHOLD,
        const double RADIUS,
        const PointToKeyInterface &POINT_TO_KEY,
        const size_t MIN_NUMBER_NEIGHBORS,
        const double MIN_NEIGHBOR_DISTANCE
      ) const;

    }; // class GeometricHashStorageClasses

    //! @brief construct on access function for all GeometricHashStorageClasses
    //! @return reference to only instances of GeometricHashStorageClasses
    BCL_API
    GeometricHashStorageClasses &GetGeometricHashStorageClasses();

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< coord::GeometricHashStorageInterface>, coord::GeometricHashStorageClasses>;

  } // namespace util
} // namespace bcl

#endif // BCL_COORD_GEOMETRIC_HASH_STORAGE_CLASSES_H_
