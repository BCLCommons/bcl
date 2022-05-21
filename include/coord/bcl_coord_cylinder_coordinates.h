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

#ifndef BCL_COORD_CYLINDER_COORDINATES_H_
#define BCL_COORD_CYLINDER_COORDINATES_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CylinderCoordinates
    //! @brief Class for representing cylindrical coordinates for a body
    //!
    //! @see @link example_coord_cylinder_coordinates.cpp @endlink
    //! @author karakam, woetzen
    //! @date: 11/11/2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CylinderCoordinates :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      double m_Height;  //!< Height in the cylinder
      double m_Radius;  //!< Distance from main axis
      double m_Angle;   //!< Angle from origin in radians

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief the default constructor
      CylinderCoordinates();

      //! @brief construct CylinderCoordinates from sizes
      //! @param HEIGHT height
      //! @param RADIUS radius
      //! @param ANGLE angle
      CylinderCoordinates( const double HEIGHT, const double RADIUS, const double ANGLE);

      //! @brief copy constructor
      //! @return pointer to a new CylinderCoordinates copied from this one
      CylinderCoordinates *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return height
      //! @return height
      double GetHeight() const
      {
        return m_Height;
      }

      //! @brief return radius
      //! @return radius
      double GetRadius() const
      {
        return m_Radius;
      }

      //! @brief return angle in radians
      //! @return angle in radians
      double GetAngle() const
      {
        return m_Angle;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief generate Cartesian coordinates from cylindrical coordinates
      //! @return Cartesian coordinates generated from cylindrical coordinates
      linal::Vector3D GetCartesianCoordinates() const;

    //////////////
    // operator //
    //////////////

      //! @brief operator ==
      //! @param RHS right hand side Cylinder coordinate to compare to
      //! @return true if height, radius and angle are the same
      bool operator ==( const CylinderCoordinates &RHS) const
      {
        return m_Height == RHS.m_Height && m_Radius == RHS.m_Radius && m_Angle == RHS.m_Angle;
      }

      //! @brief operator !=
      //! @param RHS right hand side Cylinder coordinate to compare to
      //! @return true if any of height, radius and angle are different
      bool operator !=( const CylinderCoordinates &RHS) const
      {
        return !operator ==( RHS);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read CylinderCoordinates from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write CylinderCoordinates into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class CylinderCoordinates

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_CYLINDER_COORDINATES_H_
