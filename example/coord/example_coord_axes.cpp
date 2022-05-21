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
#include "coord/bcl_coord_axes.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_axes.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordAxes :
    public ExampleInterface
  {
  public:

    ExampleCoordAxes *Clone() const
    {
      return new ExampleCoordAxes( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // three enums for x, y and z axis
      coord::Axis x_axis( coord::GetAxes().e_X);
      coord::Axis y_axis( coord::GetAxes().e_Y);
      coord::Axis z_axis( coord::GetAxes().e_Z);

      // check that they have a 1.0 at the correct place
      BCL_Example_Check
      (
        x_axis->operator()( 0) == 1.0, "pos 0 for x axis should be 1.0"
      );
      BCL_Example_Check
      (
        y_axis->operator()( 1) == 1.0, "pos 1 for y axis should be 1.0"
      );
      BCL_Example_Check
      (
        z_axis->operator()( 2) == 1.0, "pos 2 for z axis should be 1.0"
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // undefined
      BCL_Example_Check
      (
        coord::GetAxes().e_Undefined->IsDefined(),
        "undefined axes should be undefined"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write axis to file
      BCL_MessageStd( "testing read and write functionalities");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( y_axis);
      BCL_MessageVrb( "read object");
      coord::Axis read_axis;
      ReadBCLObject( read_axis);

      BCL_Example_Check
      (
        y_axis == read_axis,
        "written and read axis is different"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordAxes

  const ExampleClass::EnumType ExampleCoordAxes::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordAxes())
  );

} // namespace bcl
