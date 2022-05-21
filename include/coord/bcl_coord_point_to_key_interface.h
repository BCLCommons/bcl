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

#ifndef BCL_COORD_POINT_TO_KEY_INTERFACE_H_
#define BCL_COORD_POINT_TO_KEY_INTERFACE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeyInterface
    //! @brief class, converts a Vector3D to a Triplet of int
    //! @details this can be used to quantize cartesian points onto a quantized grid (assign the point to a particular space bin)
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Nov 5, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PointToKeyInterface :
      public util::FunctionInterface< linal::Vector3D, storage::Triplet< int, int, int> >
    {

    public:

    //////////
    // data //
    //////////

      //! @brief default angular resolution
      static double GetDefaultAngularResolution();

      //! @brief default distance resolution
      static double GetDefaultDistanceResolution();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! virtual copy constructor
      virtual PointToKeyInterface *Clone() const = 0;

      //! virtual destructor
      virtual ~PointToKeyInterface();

    /////////////////
    // data access //
    /////////////////

      //! @brief return the CoordinateSystem as string
      //! @return string describing the coordinate system
      virtual const std::string &GetCoordinateSystem() const = 0;

      //! @brief get angular resolution
      //! @return angular resolution
      virtual double GetAngularResolution() const = 0;

      //! @brief set angular resolution
      //! @param ANGULAR_RESOLUTION for this convert function
      virtual void SetAngularResolution( const double ANGULAR_RESOLUTION) = 0;

      //! @brief get distance resolution
      //! @return distance resolution
      virtual double GetDistanceResolution() const = 0;

      //! @brief set distance resolution
      //! @param DISTANCE_RESOLUTION for this convert function
      virtual void SetDistanceResolution( const double DISTANCE_RESOLUTION) = 0;

    ///////////////
    // operators //
    ///////////////

      //! operator converting point to a triplet of ints
      virtual storage::Triplet< int, int, int> operator()( const linal::Vector3D &POINT) const = 0;

      //! operator converting Triplet to a point
      virtual linal::Vector3D operator()( const storage::Triplet< int, int, int> &TRIPLET) const = 0;

    }; // class PointToKeyInterface

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_POINT_TO_KEY_INTERFACE_H_ 
