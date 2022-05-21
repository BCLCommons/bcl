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
#include "coord/bcl_coord_line_segment_2d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_2d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_line_segment_2d.cpp
  //!
  //! @author mendenjl
  //! @date Dec 05, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordLineSegment2D :
    public ExampleInterface
  {
  public:

    ExampleCoordLineSegment2D *Clone() const
    {
      return new ExampleCoordLineSegment2D( *this);
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

      // initialize some helper variables
      const linal::Vector2D start_point( 0.0, 0.0);
      const linal::Vector2D end_point( 5.0, 5.0);
      const linal::Vector2D end_point_6( 6.0, 0.0);
      const double length_5( end_point.Norm());
      const double length_6( 6.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      coord::LineSegment2D line_segment_default;

      // construct from start and end point
      coord::LineSegment2D line_segment_1( start_point, end_point);

      // copy constructor
      coord::LineSegment2D line_segment_copy( line_segment_1);

    /////////////////
    // data access //
    /////////////////

      // start, end point and direction
      BCL_ExampleCheck( line_segment_1.GetStartPoint(), start_point);
      BCL_ExampleCheck( line_segment_1.GetEndPoint(), end_point);
      BCL_ExampleCheck( line_segment_1.GetDirection(), end_point);

      // set the start point
      line_segment_default.SetStartPoint( start_point);

      // set the end point
      line_segment_default.SetEndPoint( end_point_6);

      // start, end point and direction
      BCL_ExampleCheck( line_segment_default.GetStartPoint(), start_point);
      BCL_ExampleCheck( line_segment_default.GetEndPoint(), end_point_6);
      BCL_ExampleCheck( line_segment_default.GetDirection(), end_point_6);

      // set direction
      line_segment_default.SetDirection( end_point);
      BCL_ExampleIndirectCheck( line_segment_default.GetDirection(), end_point, "SetDirection");
      line_segment_default.SetDirection( end_point_6);

      // Check intersects
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).DistanceToPoint( end_point),
        0.0
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).DistanceToPoint( start_point),
        0.0
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        coord::LineSegment2D( start_point, end_point).DistanceToPoint( linal::Vector2D( 2.5, 5.0)),
        1.76777,
        0.001
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).Overlaps
        (
          coord::LineSegment2D( start_point, 2.0 * end_point)
        ),
        true
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).Overlaps
        (
          coord::LineSegment2D( start_point, end_point)
        ),
        true
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).Overlaps
        (
          coord::LineSegment2D( start_point, -end_point)
        ),
        false
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).DoesIntersect
        (
          coord::LineSegment2D( start_point, -end_point)
        ),
        true
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( start_point, end_point).ShortestLineBetween
        (
          coord::LineSegment2D( start_point, -end_point)
        ).First(),
        coord::LineSegment2D( linal::Vector2D( 0.0, 0.0), linal::Vector2D( 0.0, 0.0))
      );
      BCL_ExampleCheck
      (
        coord::LineSegment2D( linal::Vector2D( 0.0, 0.0), linal::Vector2D( 5.0, 5.0)).ShortestLineBetween
        (
          coord::LineSegment2D( linal::Vector2D( 1.0, 0.0), linal::Vector2D( 10.0, 5.0))
        ).First(),
        coord::LineSegment2D( linal::Vector2D( 0.5, 0.5), linal::Vector2D( 1.0, 0.0))
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // length
      BCL_ExampleCheck( line_segment_1.GetLength(), length_5);
      BCL_ExampleCheck( line_segment_default.GetLength(), length_6);

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( line_segment_1);
      // read
      coord::LineSegment2D line_segment_read;
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

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordLineSegment2D

  const ExampleClass::EnumType ExampleCoordLineSegment2D::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordLineSegment2D())
  );

} // namespace bcl

