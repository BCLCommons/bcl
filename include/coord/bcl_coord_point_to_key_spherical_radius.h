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

#ifndef BCL_COORD_POINT_TO_KEY_SPHERICAL_RADIUS_H_
#define BCL_COORD_POINT_TO_KEY_SPHERICAL_RADIUS_H_

// include the namespace header
#include "bcl_coord.h"

// includes from bcl - sorted alphabetically
#include "bcl_coord_point_to_key_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeySphericalRadius
    //! @brief PointToKeySphericalRadius class, operator will quantize point to three ints according to a Spherical coordinate frame
    //!
    //! @see @link example_coord_point_to_key_spherical_radius.cpp @endlink
    //! @author woetzen
    //! @date Aug 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PointToKeySphericalRadius :
      public PointToKeyInterface
    {
    public:

    //////////
    // data //
    //////////

      //! @brief the descriptor for this coordinate system
      //! @return the name for this coordinate system that can be used for the enumerator
      static const std::string &GetCoordinateSystemDescriptor();

    private:

      //! @brief angular resolution
      double m_AngularResolution;

      //! @brief distance resolution
      double m_DistanceResolution;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from Resolution
      PointToKeySphericalRadius( const double ANGULAR_RESOLUTION);

      //! copy constructor
      PointToKeySphericalRadius *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the CoordinateSystem as string
      //! @return string describing the coordinate system
      const std::string &GetCoordinateSystem() const;

      //! @brief get angular resolution
      //! @return angular resolution
      double GetAngularResolution() const;

      //! @brief set angular resolution
      //! @param ANGULAR_RESOLUTION for this convert function
      void SetAngularResolution( const double ANGULAR_RESOLUTION);

      //! @brief get distance resolution
      //! @return distance resolution
      double GetDistanceResolution() const;

      //! @brief set distance resolution
      //! @param DISTANCE_RESOLUTION for this convert function
      void SetDistanceResolution( const double DISTANCE_RESOLUTION);

    ///////////////
    // operators //
    ///////////////

      //! operator will quantize point to three ints within a coordinate frame
      storage::Triplet< int, int, int> operator()( const linal::Vector3D &POINT) const;

      //! operator converting Triplet to a point
      linal::Vector3D operator()( const storage::Triplet< int, int, int> &TRIPLET) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write to OSTREAM
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read from ISTREAM
      std::istream &Read( std::istream &ISTREAM);

    }; // class PointToKeySphericalRadius

  } // namespace coord
} // namespace bcl

#endif // BCL_COORD_POINT_TO_KEY_SPHERICAL_RADIUS_H_
