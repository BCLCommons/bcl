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
#include "coord/bcl_coord_polygon.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_2d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_polygon.cpp
  //!
  //! @author mendenjl
  //! @date Dec 06, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordPolygon :
    public ExampleInterface
  {
  public:

    ExampleCoordPolygon *Clone() const
    { return new ExampleCoordPolygon( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      // construct a few example polygons
      coord::Polygon point;
      point.PushBack( linal::Vector2D( 1.0, 1.0));
      coord::Polygon line;
      line.PushBack( linal::Vector2D( -1.0, 2.0));
      line.PushBack( linal::Vector2D( 0.0, 0.0));
      coord::Polygon triangle;
      triangle.PushBack( linal::Vector2D( 1.0, 1.0));
      triangle.PushBack( linal::Vector2D( 3.0, 1.0));
      triangle.PushBack( linal::Vector2D( 2.0, 2.0));
      coord::Polygon irregular;
      irregular.PushBack( linal::Vector2D( 1.0, 1.0));
      irregular.PushBack( linal::Vector2D( 3.0, 1.0));
      irregular.PushBack( linal::Vector2D( 5.0, 2.0));
      irregular.PushBack( linal::Vector2D( 7.0, 1.5));
      irregular.PushBack( linal::Vector2D( 6.0, 3.0));
      irregular.PushBack( linal::Vector2D( -1.0, 2.0));
      irregular.PushBack( linal::Vector2D( -2.0, -2.0));

      // empty
      BCL_ExampleCheck( coord::Polygon().GetArea(), 0.0);
      BCL_ExampleCheck( coord::Polygon().GetPerimeter(), 0.0);
      BCL_ExampleCheck( coord::Polygon().IsCornerOf( linal::Vector2D( 0.0, 0.0)), false);
      BCL_ExampleCheck( coord::Polygon().IsWithin( linal::Vector2D( 0.0, 0.0)), false);

      // point
      BCL_ExampleCheck( point.GetArea(), 0.0);
      BCL_ExampleCheck( point.GetPerimeter(), 0.0);
      BCL_ExampleCheck( point.IsCornerOf( linal::Vector2D( 1.0, 1.0)), true);
      BCL_ExampleCheck( point.IsCornerOf( linal::Vector2D( 0.0, 0.0)), false);
      BCL_ExampleCheck( point.IsWithin( linal::Vector2D( 1.0, 1.0)), true);
      BCL_ExampleCheck( point.IsWithin( linal::Vector2D( 0.0, 0.0)), false);
      BCL_ExampleCheck( point.FindNearestSide( linal::Vector2D( 0.0, 0.0)).second, linal::Vector2D( 1.0, 1.0).Norm());

      // line
      BCL_ExampleCheck( line.GetArea(), 0.0);
      BCL_ExampleCheckWithinAbsTolerance( line.GetPerimeter(), math::Sqrt( 5.0), 0.01);
      BCL_ExampleCheck( line.IsCornerOf( linal::Vector2D( 1.0, 1.0)), false);
      BCL_ExampleCheck( line.IsCornerOf( linal::Vector2D( 0.0, 0.0)), true);
      BCL_ExampleCheck( line.IsWithin( linal::Vector2D( 1.0, 1.0)), false);
      BCL_ExampleCheckWithinAbsTolerance( line.FindNearestSide( linal::Vector2D( -0.5, 1.0)).second, 0.0, 0.0001);

      // triangle
      BCL_ExampleCheckWithinAbsTolerance( triangle.GetArea(), 1.0, 0.01);
      BCL_ExampleCheckWithinAbsTolerance( triangle.GetPerimeter(), 2.0 + 2.0 * math::Sqrt( 2.0), 0.01);
      BCL_ExampleCheck( triangle.IsCornerOf( linal::Vector2D( 3.0, 1.0)), true);
      BCL_ExampleCheck( triangle.IsCornerOf( linal::Vector2D( 0.0, 0.0)), false);
      BCL_ExampleCheck( triangle.IsWithin( linal::Vector2D( 3.0, 1.0)), true);
      BCL_ExampleCheck( triangle.IsWithin( linal::Vector2D( 3.0, 1.9)), false);
      BCL_ExampleCheck( triangle.IsWithin( linal::Vector2D( 2.0, 1.5)), true);
      BCL_ExampleCheckWithinAbsTolerance( triangle.FindNearestSide( linal::Vector2D( 2.0, 1.5)).second, 0.353553, 0.0001);

      // irregular
      BCL_ExampleCheckWithinAbsTolerance( irregular.GetArea(), 12.75, 0.01);
      BCL_ExampleCheckWithinAbsTolerance( irregular.GetPerimeter(), 23.5372, 0.01);
      BCL_ExampleCheck( irregular.IsCornerOf( linal::Vector2D( 3.0, 1.0)), true);
      BCL_ExampleCheck( irregular.IsCornerOf( linal::Vector2D( 0.0, 0.0)), false);
      BCL_ExampleCheck( irregular.IsWithin( linal::Vector2D( 3.0, 1.0)), true);
      BCL_ExampleCheck( irregular.IsWithin( linal::Vector2D( 3.0, 0.5)), false);
      BCL_ExampleCheck( irregular.IsWithin( linal::Vector2D( 2.0, 1.5)), true);
      BCL_ExampleCheck( irregular.IsWithin( linal::Vector2D( -1.0, -1.0)), true);
      BCL_ExampleCheckWithinAbsTolerance( irregular.FindNearestSide( linal::Vector2D( 2.0, 1.5)).second, 0.5, 0.0001);

      // generate 100 points within a circle, along with 100 points inside the circle, to check the convex hull
      // algorithm
      storage::Vector< linal::Vector2D> points_in_circle_radius_one;
      coord::Polygon actual_convex_hull;
      storage::Vector< double> phis;
      for( size_t i( 0); i < 100; ++i)
      {
        //random angle between 0 and 2pi
        const double phi( 2 * math::g_Pi * random::GetGlobalRandom().Double());
        phis.PushBack( phi);

        // get the point on the circle's circumference
        linal::Vector2D point_in_circle( std::cos( phi), std::sin( phi));

        // add the vector within the circle to the set of points in the circle
        points_in_circle_radius_one.PushBack( point_in_circle);

        // get a random distance
        const double distance( math::Sqrt( random::GetGlobalRandom().Random( 1.0)));

        // add the vector somewhere randomly within the circle to the set of points
        points_in_circle_radius_one.PushBack( point_in_circle * distance);
      }
      phis.Sort( std::less< double>());

      // compute the actual convex hull using the phis to just get the ordered points on the surface of the circle
      for( size_t i( 0); i < 100; ++i)
      {
        actual_convex_hull.PushBack( linal::Vector2D( std::cos( phis( i)), std::sin( phis( i))));
      }

      // compute the convex hull using the internal algorithm
      coord::Polygon computed_convex_hull( coord::Polygon::ConvexHull( points_in_circle_radius_one));

      // check perimeters and areas, which should be identical
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        actual_convex_hull.GetArea(),
        computed_convex_hull.GetArea(),
        0.01,
        "ConvexHull"
      );
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        actual_convex_hull.GetPerimeter(),
        computed_convex_hull.GetPerimeter(),
        0.01,
        "ConvexHull"
      );
      BCL_ExampleIndirectCheck
      (
        actual_convex_hull.GetNumberOfSides(),
        computed_convex_hull.GetNumberOfSides(),
        "ConvexHull"
      );

      // check that all points in the phis vector were within the convex hull
      bool all_points_within_convex_hull( true);
      for
      (
        storage::Vector< linal::Vector2D>::const_iterator
          itr( points_in_circle_radius_one.Begin()), itr_end( points_in_circle_radius_one.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !computed_convex_hull.IsWithin( *itr))
        {
          all_points_within_convex_hull = false;
          break;
        }
      }
      BCL_ExampleIndirectCheck
      (
        all_points_within_convex_hull,
        true,
        "computed convex hull should contain all the points"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleCoordPolygon

  const ExampleClass::EnumType ExampleCoordPolygon::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordPolygon())
  );

} // namespace bcl
