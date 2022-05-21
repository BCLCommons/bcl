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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "coord/bcl_coord_cylinder_coordinates.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_cylinder_coordinates.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordCylinderCoordinates :
    public ExampleInterface
  {
  public:

    ExampleCoordCylinderCoordinates *Clone() const
    { return new ExampleCoordCylinderCoordinates( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////////
  // helper function //
  /////////////////////

    //! compares vectors like operator ==, but returns true if difference between both is smaller than EPSILON
    static bool EqualWithinTolerance( const coord::CylinderCoordinates &CYLINDER_A, const coord::CylinderCoordinates &CYLINDER_B, const double &RELATIVE_TOLERANCE = 0.001, const double &ABSOLUTE_TOLERANCE = std::numeric_limits< double>::epsilon())
    {
      if( !math::EqualWithinTolerance( CYLINDER_A.GetHeight(), CYLINDER_B.GetHeight(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE))
      {
        return false;
      }
      if( !math::EqualWithinTolerance( CYLINDER_A.GetRadius(), CYLINDER_B.GetRadius(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE))
      {
        return false;
      }
      if( !math::EqualWithinTolerance( CYLINDER_A.GetAngle(), CYLINDER_B.GetAngle(), RELATIVE_TOLERANCE, ABSOLUTE_TOLERANCE))
      {
        return false;
      }

      return true;
    }

    int Run() const
    {
    //////////////////
    // constructors //
    //////////////////

      // default constructor - initializes everything to zero
      coord::CylinderCoordinates cc_default;

      BCL_Example_Check
      (
        cc_default.GetAngle() == double( 0.0)
        && cc_default.GetHeight() == double( 0.0)
        && cc_default.GetRadius() == double( 0.0),
        "Default constructor should initialize everything to zero"
      );

      // construct from height, radius and angle
      coord::CylinderCoordinates cc_hra( double( 1.1), double( 2.2), math::g_Pi);
      BCL_Example_Check
      (
        cc_hra.GetAngle() == math::g_Pi
        && cc_hra.GetHeight() == double( 1.1)
        && cc_hra.GetRadius() == double( 2.2),
        "constructor from height radius and angle does not set members correctly: " + util::Format()( cc_hra)
      );

      // copy constructor
      coord::CylinderCoordinates cc_copy( cc_hra);
      BCL_Example_Check
      (
        cc_hra.GetAngle() == math::g_Pi
        && cc_hra.GetHeight() == double( 1.1)
        && cc_hra.GetRadius() == double( 2.2),
        "copy constructor does not set members correctly: " + util::Format()( cc_copy)
      );

      // clone
      util::ObjectInterface *ptr( cc_hra.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd( "Cylindrical height is " + util::Format()( cc_hra.GetHeight()));
      BCL_MessageStd( "Cylindrical radius is " + util::Format()( cc_hra.GetRadius()));
      BCL_MessageStd( "Cylindrical angle is " + util::Format()( cc_hra.GetAngle()));

      BCL_MessageStd( "static name of class is is " + GetStaticClassName< coord::CylinderCoordinates>());
      BCL_MessageStd( "class identifier is " + ptr->GetClassIdentifier());

      BCL_Example_Check
      (
        GetStaticClassName< coord::CylinderCoordinates>() == ptr->GetClassIdentifier(),
        "incorrect class identifier"
      );

    ////////////////
    // operations //
    ////////////////

      // convert into Cartesian coordinates
      BCL_MessageStd
      (
        "converted into Cartesian coordinates: " + util::Format()( cc_hra.GetCartesianCoordinates())
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          cc_hra.GetCartesianCoordinates(),
          linal::Vector3D( double( -2.2), double( 0.0), double( 1.1)),
          double( 0.0), double( 0.000001)
        ),
        "conversion into Cartesian coordinates failed: " + util::Format()( cc_hra) +
        " -> " + util::Format()( cc_hra.GetCartesianCoordinates())
      );

    ///////////////
    // operators //
    ///////////////

      BCL_Example_Check
      (
        cc_default != cc_hra,
        "!= operator yields true for: " + util::Format()( cc_default) + " != " + util::Format()( cc_hra)
      );

      BCL_Example_Check
      (
        cc_copy == cc_hra,
        "!= operator yields true for: " + util::Format()( cc_copy) + " != " + util::Format()( cc_hra)
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( cc_hra);
      BCL_MessageVrb( "read object");
      coord::CylinderCoordinates cc_read;
      ReadBCLObject( cc_read);

      BCL_Example_Check
      (
        EqualWithinTolerance( cc_read, cc_hra),
        "written cylinder coordinates are not equal to read cylinder coordinates " +
        util::Format()( cc_hra) + " != " + util::Format()( cc_read)
      );

      delete ptr;

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordCylinderCoordinates

  const ExampleClass::EnumType ExampleCoordCylinderCoordinates::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordCylinderCoordinates())
  );

} // namespace bcl
