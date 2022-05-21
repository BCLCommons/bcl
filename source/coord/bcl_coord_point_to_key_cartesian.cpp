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
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeyCartesian
    //! @brief operator will quantize point to three ints according to a Cartesian frame
    //!
    //! @see @link example_coord_point_to_key_cartesian.cpp @endlink
    //! @author woetzen
    //! @date Aug 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API PointToKeyCartesian :
      public PointToKeyInterface
    {
    private:
    //////////
    // data //
    //////////

      static const std::string &GetCoordinateSystemDescriptor();

      //! @brief distance resolution
      double m_DistanceResolution;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! instance of the PointToKey enum
      static PointToKey e_Cartesian;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from Resolution
      PointToKeyCartesian( const double DISTANCE_RESOLUTION);

      //! virtual copy constructor
      PointToKeyCartesian *Clone() const;

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

    }; // class PointToKeyCartesian

  //////////
  // data //
  //////////

    const std::string &PointToKeyCartesian::GetCoordinateSystemDescriptor()
    {
      static const std::string s_coordinate_system( "Cartesian");
      return s_coordinate_system;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointToKeyCartesian::s_Instance
    (
      GetObjectInstances().AddInstance( new PointToKeyCartesian( PointToKeyInterface::GetDefaultDistanceResolution()))
    );

    //! instance of the PointToKey enum
    PointToKey PointToKeyCartesian::e_Cartesian
    (
      GetPointToKeyClasses().AddEnum( PointToKeyCartesian::GetCoordinateSystemDescriptor(), util::ShPtr< PointToKeyInterface>( new PointToKeyCartesian( PointToKeyInterface::GetDefaultDistanceResolution())))
    );

    //! construct from Resolution
    PointToKeyCartesian::PointToKeyCartesian( const double DISTANCE_RESOLUTION) :
      m_DistanceResolution( DISTANCE_RESOLUTION)
    {
    }

    //! virtual copy constructor
    PointToKeyCartesian *PointToKeyCartesian::Clone() const
    {
      return new PointToKeyCartesian( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeyCartesian::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the CoordinateSystem as string
    //! @return string describing the coordinate system
    const std::string &PointToKeyCartesian::GetCoordinateSystem() const
    {
      // return
      return GetCoordinateSystemDescriptor();
    }

    //! @brief get angular resolution
    //! @return angular resolution
    double PointToKeyCartesian::GetAngularResolution() const
    {
      return util::GetUndefinedSize_t();
    }

    //! @brief set angular resolution
    //! @param ANGULAR_RESOLUTION for this convert function
    void PointToKeyCartesian::SetAngularResolution( const double ANGULAR_RESOLUTION)
    {
      return;
    }

    //! @brief get distance resolution
    //! @return distance resolution
    double PointToKeyCartesian::GetDistanceResolution() const
    {
      return m_DistanceResolution;
    }

    //! @brief set distance resolution
    //! @param DISTANCE_RESOLUTION for this convert function
    void PointToKeyCartesian::SetDistanceResolution( const double DISTANCE_RESOLUTION)
    {
      m_DistanceResolution = DISTANCE_RESOLUTION;
    }

    //! operator will quantize point to three ints within a coordinate frame
    storage::Triplet< int, int, int> PointToKeyCartesian::operator()( const linal::Vector3D &POINT) const
    {
      // convert x y z to size_t and convert to hash key
      return storage::Triplet< int, int, int>
             (
               int( POINT.X() / m_DistanceResolution),
               int( POINT.Y() / m_DistanceResolution),
               int( POINT.Z() / m_DistanceResolution)
             );
    }

    //! operator converting Triplet to a point
    linal::Vector3D PointToKeyCartesian::operator()( const storage::Triplet< int, int, int> &TRIPLET) const
    {
      return linal::Vector3D
             (
               double( TRIPLET.First())  * m_DistanceResolution,
               double( TRIPLET.Second()) * m_DistanceResolution,
               double( TRIPLET.Third())  * m_DistanceResolution
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write to OSTREAM
    std::ostream &PointToKeyCartesian::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_DistanceResolution, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read from ISTREAM
    std::istream &PointToKeyCartesian::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_DistanceResolution, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace coord
} // namespace bcl
