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

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "coord/bcl_coord_axes.h"
#include "coord/bcl_coord_point_to_key_classes.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeySpherical
    //! @brief perator will quantize point to three ints according to a Spherical coordinate frame
    //!
    //! @see @link example_coord_point_to_key_spherical.cpp @endlink
    //! @author woetzen
    //! @date Aug 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API PointToKeySpherical :
      public PointToKeyInterface
    {
    private:

    //////////
    // data //
    //////////

      static const std::string &GetCoordinateSystemDescriptor();

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

      //! instance of the PointToKey enum
      static PointToKey e_Spherical;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from Resolution
      PointToKeySpherical( const double ANGULAR_RESOLUTION);

      //! copy constructor
      PointToKeySpherical *Clone() const;

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

      //! operator will quantize point to three ints within a coordinate frame
      storage::Triplet< int, int, int> operator()( const linal::Vector3D &POINT) const;

      //! operator converting Triplet to a point
      linal::Vector3D operator()( const storage::Triplet< int, int, int> &TRIPLET) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! write to OSTREAM
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read from ISTREAM
      std::istream &Read( std::istream &ISTREAM);

    }; // class PointToKeySpherical

  //////////
  // data //
  //////////

    const std::string &PointToKeySpherical::GetCoordinateSystemDescriptor()
    {
      static const std::string s_coordinate_system( "Spherical");
      return s_coordinate_system;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointToKeySpherical::s_Instance
    (
      GetObjectInstances().AddInstance( new PointToKeySpherical( PointToKeyInterface::GetDefaultAngularResolution()))
    );

    //! instance of the PointToKey enum
    PointToKey PointToKeySpherical::e_Spherical
    (
      GetPointToKeyClasses().AddEnum( PointToKeySpherical::GetCoordinateSystemDescriptor(), util::ShPtr< PointToKeyInterface>( new PointToKeySpherical( PointToKeyInterface::GetDefaultAngularResolution())))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct from Resolution
    PointToKeySpherical::PointToKeySpherical( const double ANGULAR_RESOLUTION) :
      m_AngularResolution( ANGULAR_RESOLUTION),
      m_DistanceResolution( PointToKeyInterface::GetDefaultDistanceResolution())
    {
    }

    //! virtual copy constructor
    PointToKeySpherical *PointToKeySpherical::Clone() const
    {
      return new PointToKeySpherical( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeySpherical::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the CoordinateSystem as string
    //! @return string describing the coordinate system
    const std::string &PointToKeySpherical::GetCoordinateSystem() const
    {
      // return
      return GetCoordinateSystemDescriptor();
    }

    //! @brief get angular resolution
    //! @return angular resolution
    double PointToKeySpherical::GetAngularResolution() const
    {
      return m_AngularResolution;
    }

    //! @brief set angular resolution
    //! @param ANGULAR_RESOLUTION for this convert function
    void PointToKeySpherical::SetAngularResolution( const double ANGULAR_RESOLUTION)
    {
      m_AngularResolution = ANGULAR_RESOLUTION;
    }

    //! @brief get distance resolution
    //! @return distance resolution
    double PointToKeySpherical::GetDistanceResolution() const
    {
      return m_DistanceResolution;
    }

    //! @brief set distance resolution
    //! @param DISTANCE_RESOLUTION for this convert function
    void PointToKeySpherical::SetDistanceResolution( const double DISTANCE_RESOLUTION)
    {
      m_DistanceResolution = DISTANCE_RESOLUTION;
    }

    //! operator will quantize point to three ints within a coordinate frame
    storage::Triplet< int, int, int> PointToKeySpherical::operator()( const linal::Vector3D &POINT) const
    {
      // derive double values
      const double len( POINT.Norm());                   // length
      const double costheta( POINT.Z() / len);           // costheta
      const double phi( atan2( POINT.Y(), POINT.X()));   // phi
      const double theta( acos( costheta) / math::g_Pi); // theta

      // convert to integers
      const size_t lenint = size_t( len / m_DistanceResolution);
      const size_t thetaint = size_t( theta * m_AngularResolution);
      const size_t points = size_t( sin( math::g_Pi * thetaint / m_AngularResolution) * m_AngularResolution); // number phi points
      const size_t phiint = size_t( ( phi / math::g_Pi) * ( points + 1));

      // convert to hash key
      return storage::Triplet< int, int, int>( int( lenint), int( thetaint), int( phiint));
    }

    //! operator converting Triplet to a point
    linal::Vector3D PointToKeySpherical::operator()( const storage::Triplet< int, int, int> &TRIPLET) const
    {
      // get integer values
      const double lenint( TRIPLET.First());
      const double thetaint( TRIPLET.Second());
      const double phiint( TRIPLET.Third());

      // Calculate angles and distances
      const double len( lenint * m_DistanceResolution);
      const double theta( thetaint / m_AngularResolution);
      const size_t points( size_t( sin( math::g_Pi * thetaint / m_AngularResolution) * m_AngularResolution));
      const double phi( phiint * math::g_Pi / double( points + 1));

      // create the point and set the coordinates
      linal::Vector3D point;
      const double costheta( cos( theta * math::g_Pi));
      const double z( costheta * len);
      const double hypo( math::Sqrt( math::Sqr( len) - math::Sqr( z))); // = r*sin(theta) just faster to calculate that way when z and len is given
      const double x( hypo * cos( phi));
      const double y( hypo * sin( phi));

      //end
      return linal::Vector3D( x, y, z);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write to OSTREAM
    std::ostream &PointToKeySpherical::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AngularResolution, OSTREAM, INDENT);
      io::Serialize::Write( m_DistanceResolution, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read from ISTREAM
    std::istream &PointToKeySpherical::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AngularResolution, ISTREAM);
      io::Serialize::Read( m_DistanceResolution, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace coord
} // namespace bcl
