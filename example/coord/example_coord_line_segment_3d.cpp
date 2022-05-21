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
#include "coord/bcl_coord_line_segment_3d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_line_segment_3d.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordLineSegment3D :
    public ExampleInterface
  {
  public:

    ExampleCoordLineSegment3D *Clone() const
    {
      return new ExampleCoordLineSegment3D( *this);
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
      // TODO Check global functions: CalculateDistancePointFromLineSegment and ShortestConnectionBetweenLineSegments3D

      // initialize some helper variables
      const linal::Vector3D start_point( 0.0, 0.0, 0.0);
      const linal::Vector3D end_point( 5.0, 0.0, 0.0);
      const linal::Vector3D end_point_6( 6.0, 0.0, 0.0);
      const double length_5( 5.0);
      const double length_6( 6.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      coord::LineSegment3D line_segment_default;

      // construct from start and end point
      coord::LineSegment3D line_segment_1( start_point, end_point);

      // copy constructor
      coord::LineSegment3D line_segment_copy( line_segment_1);

      // clone
      coord::LineSegment3D *ptr( line_segment_1.Clone());

    /////////////////
    // data access //
    /////////////////

      // class name and identifier
      BCL_MessageStd( "class name: " + ptr->GetClassIdentifier());
      BCL_Example_Check
      (
        GetStaticClassName< coord::LineSegment3D>() == "bcl::coord::LineSegment3D"
        && ptr->GetClassIdentifier() == GetStaticClassName< coord::LineSegment3D>(),
          "incorrect class name: static class name: " + GetStaticClassName< coord::LineSegment3D>()
        + " class identifier: " + ptr->GetClassIdentifier()
      );

      // start, end point and direction
      BCL_Example_Check
      (
        line_segment_1.GetStartPoint() == start_point,
          "start point should be " + util::Format()( start_point) + " but is "
        + util::Format()( line_segment_1.GetStartPoint())
      );

      BCL_Example_Check
      (
        line_segment_1.GetEndPoint() == end_point,
          "end point should be " + util::Format()( end_point) + " but is "
        + util::Format()( line_segment_1.GetEndPoint())
      );

      BCL_Example_Check
      (
        line_segment_1.GetDirection() == end_point,
          "direction should be " + util::Format()( end_point) + " but is "
        + util::Format()( line_segment_1.GetDirection())
      );

      // set the start point
      line_segment_default.SetStartPoint( start_point);

      // set the end point
      line_segment_default.SetEndPoint( end_point_6);

      // start, end point and direction
      BCL_Example_Check
      (
        line_segment_default.GetStartPoint() == start_point,
          "start point should be " + util::Format()( start_point) + " but is "
        + util::Format()( line_segment_default.GetStartPoint())
      );

      BCL_Example_Check
      (
        line_segment_default.GetEndPoint() == end_point_6,
          "end point should be " + util::Format()( end_point_6) + " but is "
        + util::Format()( line_segment_default.GetEndPoint())
      );

      BCL_Example_Check
      (
        line_segment_default.GetDirection() == end_point_6,
          "direction should be " + util::Format()( end_point_6) + " but is "
        + util::Format()( line_segment_default.GetDirection())
      );

      // set direction
      line_segment_default.SetDirection( end_point);

      BCL_Example_Check
      (
        line_segment_default.GetDirection() == end_point,
          "direction should be " + util::Format()( end_point) + " but is "
        + util::Format()( line_segment_default.GetDirection())
      );
      line_segment_default.SetDirection( end_point_6);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // length
      BCL_Example_Check
      (
        line_segment_1.GetLength() == length_5,
        "length should be " + util::Format()( length_5) + " but is " + util::Format()( line_segment_default.GetLength())
      );
      BCL_Example_Check
      (
        line_segment_default.GetLength() == length_6,
        "length should be " + util::Format()( length_6) + " but is " + util::Format()( line_segment_default.GetLength())
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( line_segment_1);
      // read
      coord::LineSegment3D line_segment_read;
      ReadBCLObject( line_segment_read);

      BCL_Example_Check
      (
        line_segment_read.GetStartPoint() == line_segment_1.GetStartPoint()
        && line_segment_read.GetEndPoint() == line_segment_1.GetEndPoint()
        && line_segment_read.GetDirection() == line_segment_1.GetDirection(),
        "written and read line segment is different"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // clean up
      delete ptr;

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordLineSegment3D

  const ExampleClass::EnumType ExampleCoordLineSegment3D::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordLineSegment3D())
  );

} // namespace bcl

