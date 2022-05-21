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

#ifndef BCL_COORD_POINT_TO_KEY_CLASSES_H_
#define BCL_COORD_POINT_TO_KEY_CLASSES_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_coord_point_to_key_interface.h"
#include "util/bcl_util_enumerate.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeyClasses
    //! @brief enumerates all point to key functions
    //!
    //! @see @link example_coord_point_to_key_classes.cpp @endlink
    //! @author woetzen
    //! @date 17-08-2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PointToKeyClasses :
      public util::Enumerate< util::ShPtr< PointToKeyInterface>, PointToKeyClasses>
    {
      friend class util::Enumerate< util::ShPtr< PointToKeyInterface>, PointToKeyClasses>;

    public:

    //////////
    // data //
    //////////

      PointToKey e_SphericalRadius;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PointToKeyClasses();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get point to key function for given name
      //! @param NAME
      //! @param ANGULAR_RESOLUTION
      //! @param DISTANCE_RESOLUTION
      //! @return ShPtr to PointToKey function
      util::ShPtr< PointToKeyInterface>
      GetPointToKeyFunction
      (
        const std::string &NAME,
        const double ANGULAR_RESOLUTION = PointToKeyInterface::GetDefaultAngularResolution(),
        const double DISTANCE_RESOLUTION = PointToKeyInterface::GetDefaultDistanceResolution()
      ) const;

    }; //class PointToKeyClasses

    BCL_API
    PointToKeyClasses &GetPointToKeyClasses();

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< coord::PointToKeyInterface>, coord::PointToKeyClasses>;

  } // namespace util
} // namespace bcl

#endif // BCL_COORD_POINT_TO_KEY_CLASSES_H_ 
