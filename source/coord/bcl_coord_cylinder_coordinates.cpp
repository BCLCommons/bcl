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
#include "coord/bcl_coord_cylinder_coordinates.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CylinderCoordinates::s_Instance
    (
      GetObjectInstances().AddInstance( new CylinderCoordinates())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CylinderCoordinates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief the default constructor
    CylinderCoordinates::CylinderCoordinates() :
      m_Height( double( 0.0)),
      m_Radius( double( 0.0)),
      m_Angle( double( 0.0))
    {
    }

    //! @brief construct CylinderCoordinates from sizes
    //! @param HEIGHT height
    //! @param RADIUS radius
    //! @param ANGLE angle
    CylinderCoordinates::CylinderCoordinates( const double HEIGHT, const double RADIUS, const double ANGLE) :
      m_Height( HEIGHT),
      m_Radius( RADIUS),
      m_Angle( ANGLE)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new CylinderCoordinates copied from this one
    CylinderCoordinates *CylinderCoordinates::Clone() const
    {
      return new CylinderCoordinates( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate cartesian coordinates from cylinderical coordinates
    //! @return Cartesian coordinates generated from cylinderical coordinates
    linal::Vector3D CylinderCoordinates::GetCartesianCoordinates() const
    {
      return linal::Vector3D
      (
        GetRadius() * cos( GetAngle()),
        GetRadius() * sin( GetAngle()),
        GetHeight()
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read CylinderCoordinates from std::istream
    std::istream &CylinderCoordinates::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Height, ISTREAM);
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_Angle, ISTREAM);

      // end
      return ISTREAM;
    }

    //! write CylinderCoordinates into std::ostream
    std::ostream &CylinderCoordinates::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Height, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Radius, OSTREAM) << '\t';
      io::Serialize::Write( m_Angle, OSTREAM);

      // end
      return OSTREAM;
     }

  } // namespace coord
} // namespace bcl

