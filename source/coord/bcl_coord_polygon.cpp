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
#include "coord/bcl_coord_polygon.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_path.h"
#include "linal/bcl_linal_vector_2d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_limits.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Polygon::s_Instance
    (
      GetObjectInstances().AddInstance( new Polygon())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    Polygon::Polygon()
    {
    }

    //! construct from points
    Polygon::Polygon( const storage::Vector< linal::Vector2D> &POINTS)
    {
      for
      (
        storage::Vector< linal::Vector2D>::const_iterator itr( POINTS.Begin()), itr_end( POINTS.End());
        itr != itr_end;
        ++itr
      )
      {
        PushBack( *itr);
      }
    }

    //! copy constructor
    Polygon *Polygon::Clone() const
    {
      return new Polygon( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Polygon::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns number of sides of the polygon
    //! @return number of sides of the polygon
    size_t Polygon::GetNumberOfSides() const
    {
      return
        m_Sides.GetSize() != size_t( 1)
        ? m_Sides.GetSize()
        : m_Sides.FirstElement().GetStartPoint() != m_Sides.FirstElement().GetEndPoint();
    }

    //! @brief returns number of sides of the polygon
    //! @return perimeter of the polygon
    double Polygon::GetPerimeter() const
    {
      double perimeter( 0.0);
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        perimeter += itr->GetDirection().Norm();
      }
      return perimeter;
    }

    //! @brief get the area of the polygon
    //! @return area of the polygon
    double Polygon::GetArea() const
    {
      // handle trivial cases 1st
      if( m_Sides.GetSize() < size_t( 2))
      {
        return 0.0;
      }

      double area( 0.0);
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        area += linal::CrossProduct( itr->GetStartPoint(), itr->GetEndPoint());
      }
      return math::Absolute( area * 0.5);
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    Polygon::iterator Polygon::Begin()
    {
      return m_Sides.Begin();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    Polygon::iterator Polygon::End()
    {
      return m_Sides.End();
    }

    //! @brief return iterator on begin
    //! @return iterator pointing to the beginning of the container, i.e. the first element
    Polygon::const_iterator Polygon::Begin() const
    {
      return m_Sides.Begin();
    }

    //! @brief return iterator on end
    //! @return iterator pointing to the end of the container, i.e. behind the last element
    Polygon::const_iterator Polygon::End() const
    {
      return m_Sides.End();
    }

    //! @brief add a line between the last point in the polygon and the first point
    //! @param POINT the point to add
    void Polygon::PushBack( const linal::Vector2D &POINT)
    {
      UpdateMinMaxBounds( POINT);
      if( m_Sides.IsEmpty())
      {
        // empty polygon
        m_Sides.PushBack( LineSegment2D( POINT, POINT));
      }
      else if( m_Sides.GetSize() == size_t( 1))
      {
        if( m_Sides.FirstElement().GetStartPoint() == m_Sides.FirstElement().GetEndPoint())
        {
          // starting from single vertex
          m_Sides( 0).SetEndPoint( POINT);
        }
        else
        {
          // starting from a single line, need to draw two lines, one from the end of the current line to point
          // and another from point back to the start
          m_Sides.PushBack( LineSegment2D( m_Sides.FirstElement().GetEndPoint(), POINT));
          m_Sides.PushBack( LineSegment2D( POINT, m_Sides.FirstElement().GetStartPoint()));
        }
      }
      else
      {
        m_Sides.PushBack( LineSegment2D( POINT, m_Sides.FirstElement().GetStartPoint()));
        m_Sides( m_Sides.GetSize() - 2).SetEndPoint( POINT);
      }
    }

    //! @brief conversion to a vector of line segments
    Polygon::operator const storage::Vector< LineSegment2D> &() const
    {
      return m_Sides;
    }

    //! @brief conversion to a vector of line segments
    Polygon::operator storage::Vector< LineSegment2D> &()
    {
      return m_Sides;
    }

    //! @brief Get the underlying sides of the polygon
    //! @return the sides of the polygon
    const storage::Vector< LineSegment2D> &Polygon::GetSides() const
    {
      return m_Sides;
    }

    //! @brief Get the underlying sides of the polygon
    //! @return the sides of the polygon
    storage::Vector< LineSegment2D> &Polygon::GetSides()
    {
      return m_Sides;
    }

    //! returns the geometric center of the object
    linal::Vector2D Polygon::GetCenter() const
    {
      if( m_Sides.IsEmpty())
      {
        return linal::Vector2D();
      }
      if( m_Sides.GetSize() == size_t( 1))
      {
        // only one sided, this is the only case where the average of the start points does not yield the center,
        // because the end point is not represented as a start point of any line
        return GetBarycenter();
      }
      math::RunningAverage< linal::Vector2D> center;
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        center += itr->GetStartPoint();
      }
      return center.GetAverage();
    }

    //! returns the barycenter of the object
    linal::Vector2D Polygon::GetBarycenter() const
    {
      return m_Bounds.GetStartPoint() + 0.5 * m_Bounds.GetDirection();
    }

    //! @brief test whether a particular point is a corner (intersection of two sides) of this polygon
    //! @param POINT the point to test for being a corner in this polygon
    bool Polygon::IsCornerOf( const linal::Vector2D &POINT) const
    {
      if( m_Sides.GetSize() == size_t( 1))
      {
        // only one sided, this is the only case where the average of the start points does not yield the center,
        // because the end point is not represented as a start point of any line
        return POINT == m_Sides( 0).GetStartPoint() || POINT == m_Sides( 0).GetEndPoint();
      }
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( POINT == itr->GetStartPoint())
        {
          return true;
        }
      }
      return false;
    }

    //! @brief test whether a particular point is within the polygon
    //! @param POINT the point to test for inclusion in the polygon
    bool Polygon::IsWithin( const linal::Vector2D &POINT) const
    {
      // handle the trivial case that POINT lies outside the bounds
      if
      (
        POINT.X() > m_Bounds.GetEndPoint().X()
        || POINT.Y() > m_Bounds.GetEndPoint().Y()
        || POINT.X() < m_Bounds.GetStartPoint().X()
        || POINT.Y() < m_Bounds.GetStartPoint().Y()
      )
      {
        return false;
      }

      // handle the trivial case that POINT is at one of the vertices of the polygon
      if( IsCornerOf( POINT))
      {
        return true;
      }

      // This uses the ray-casting algorithm
      // The idea is to draw a line from the point to anywhere outside the polygon and count the number of lines that
      // the line crosses.  If the number is even, the point lies outside the polygon.  If it is odd, it lies within the
      // polygon.  In the case that the polygon is convex, the number of crossings will always be 0, 1, or 2, but for
      // other polygon types, any number is possible
      size_t number_crossings( 0);

      // construct the end point of the line as slightly larger than the outer bounds of the box
      linal::Vector2D outside_polygon( m_Bounds.GetEndPoint());
      outside_polygon.X() += math::Absolute( outside_polygon.X()) * 0.1 + 0.1;
      outside_polygon.Y() += math::Absolute( outside_polygon.Y()) * 0.1 + 0.1;

      LineSegment2D from_point_to_outside_polygon( POINT, outside_polygon);

      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        if( from_point_to_outside_polygon.DoesIntersect( *itr))
        {
          if( from_point_to_outside_polygon.Overlaps( *itr))
          {
            // POINT lies along this side of the polygon; and the line we constructed happens to be collinear
            return true;
          }
          ++number_crossings;
        }
      }
      return ( number_crossings % size_t( 2)) == size_t( 1);
    }

    //! @brief Find the side nearest to a point, which may be inside or outside of the
    //! @param POINT the point to test for inclusion in the polygon
    //! @return an iterator to the nearest side and the actual distance
    std::pair< Polygon::const_iterator, double> Polygon::FindNearestSide( const linal::Vector2D &POINT) const
    {
      if( m_Sides.GetSize() == size_t( 1) && m_Sides.FirstElement().GetStartPoint() == m_Sides.FirstElement().GetEndPoint())
      {
        return std::make_pair( Begin(), linal::Distance( m_Sides( 0).GetStartPoint(), POINT));
      }
      std::pair< Polygon::const_iterator, double> nearest_side_distance( End(), math::GetHighestBoundedValue< double>());
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        const double current_line_distance( itr->DistanceToPoint( POINT));
        if( current_line_distance < nearest_side_distance.second)
        {
          nearest_side_distance.first = itr;
          nearest_side_distance.second = current_line_distance;
        }
      }
      return nearest_side_distance;
    }

    namespace
    {

      //! @brief get the next point on a 2d convex hull
      //! @param PREVIOUS the previous point on the convex hull
      //! @param POINTS the set of points that can be chosen from
      //! @return the next point on the convex hull, subject to MAX_DESIRED_DISTANCE
      linal::Vector2D GetNextConvexHullPoint
      (
        const linal::Vector2D &PREVIOUS,
        const storage::Vector< linal::Vector2D> &POINTS,
        std::string &CHOSEN_POINTS
      )
      {
        linal::Vector2D next_hull_point( PREVIOUS);
        const double previous_x( PREVIOUS.X()), previous_y( PREVIOUS.Y());
        size_t point_index( 0), next_index( 0);

        for
        (
          storage::Vector< linal::Vector2D>::const_iterator
            itr( POINTS.Begin()), itr_end( POINTS.End());
          itr != itr_end;
          ++itr, ++point_index
        )
        {
          // determine whether this point is CCW from the current next x
          const double vector_cross_product
          (
              ( next_hull_point.X() - previous_x) * ( itr->Y() - previous_y)
            - ( itr->X() - previous_x) * ( next_hull_point.Y() - previous_y)
          );

          // skip CW points
          if( vector_cross_product > 0.0)
          {
            continue;
          }

          if
          (
            // Colinear point closer to the previous point on the hull
            vector_cross_product == 0.0
            && linal::SquareDistance( *itr, PREVIOUS) <= linal::SquareDistance( next_hull_point, PREVIOUS)
          )
          {
            continue;
          }

          // accept this candidate point as the next hull point
          next_hull_point = *itr;
          next_index = point_index;
        }
        const bool already_seen( CHOSEN_POINTS[ next_index] == '1');
        CHOSEN_POINTS[ next_index] = '1';
        return already_seen ? PREVIOUS : next_hull_point;
      }
    }

    //! @brief create a polygon from the convex hull, which is the minimal convex polygon that can enclose all given points
    //! @param POINTS the points to enclose with the convex hull
    //! @return the minimal convex hull
    Polygon Polygon::ConvexHull
    (
      const storage::Vector< linal::Vector2D> &POINTS,
      const double &MAX_SIDE_LENGTH,
      const double &RADIUS
    )
    {
      if( POINTS.IsEmpty())
      {
        return Polygon();
      }

      // Find the minimum point in POINTS
      linal::Vector2D minimum( POINTS.FirstElement());
      size_t min_index( 0), index( 1);
      std::string chosen_points( POINTS.GetSize(), '0');
      for
      (
        storage::Vector< linal::Vector2D>::const_iterator
          itr( POINTS.Begin() + 1), itr_end( POINTS.End());
        itr != itr_end;
        ++itr, ++index
      )
      {
        if( *itr < minimum)
        {
          minimum = *itr;
          min_index = index;
        }
      }

      Polygon convex_hull;
      convex_hull.PushBack( minimum);
      chosen_points[ min_index] = '1';

      // Find each point on the hull in CCW order
      linal::Vector2D next_hull_point( minimum);
      do
      {
        convex_hull.PushBack( next_hull_point);
        next_hull_point = GetNextConvexHullPoint( next_hull_point, POINTS, chosen_points);
      } while
        (
          next_hull_point != minimum
          && !( next_hull_point == convex_hull.GetSides().LastElement().GetEndPoint())
          && !( next_hull_point == convex_hull.GetSides().LastElement().GetStartPoint())
        );

      if( !util::IsDefined( MAX_SIDE_LENGTH))
      {
        return convex_hull;
      }
      // determine which lines are longer than the desired distance
      const size_t initial_number_lines( convex_hull.GetSides().GetSize());

      // hash vector to determine which lines are to be split
      std::string line_must_be_split( initial_number_lines, '0');

      // determine which lines must be split
      size_t number_lines_to_split( 0);
      for( size_t line_number( 0); line_number < initial_number_lines; ++line_number)
      {
        if( convex_hull.GetSides()( line_number).GetLength() > MAX_SIDE_LENGTH)
        {
          line_must_be_split[ line_number] = '1';
          ++number_lines_to_split;
        }
      }

      const Polygon &const_convex_hull( convex_hull);
      // if no lines were longer, return the convex hull
      if( !number_lines_to_split || initial_number_lines == size_t( 1))
      {
        return convex_hull;
      }

      // for each point, determine which line it lies nearest
      // and place the point into a container for that line, if the line is to be split
      storage::Vector< storage::Vector< linal::Vector2D> > points_nearest_line( initial_number_lines);

      // get the barycenter of the initial hull
      const linal::Vector2D initial_barycenter( convex_hull.GetBarycenter());
      for
      (
        storage::Vector< linal::Vector2D>::const_iterator
          itr_points( POINTS.Begin()), itr_points_end( POINTS.End());
        itr_points != itr_points_end;
        ++itr_points
      )
      {
        // skip corner points, since they are nearest to two different lines and must be the start and end points of
        // any new paths
        if( convex_hull.IsCornerOf( *itr_points))
        {
          continue;
        }

        // find the closest line and distance from it
        std::pair< Polygon::const_iterator, double> closest_line_and_distance
        (
          convex_hull.FindNearestSide( *itr_points)
        );

        const size_t nearest_line_id( std::distance( const_convex_hull.Begin(), closest_line_and_distance.first));
        if( line_must_be_split[ nearest_line_id] == '1')
        {
          points_nearest_line( nearest_line_id).PushBack( *itr_points);
        }
      }

      // create a new polygon to contain the distance limited hull
      Polygon distance_limited_hull;

      // for each line
      const double undefined_edge_value( math::GetHighestBoundedValue< double>());
      for( size_t line_number( 0); line_number < initial_number_lines; ++line_number)
      {
        const LineSegment2D &current_side( convex_hull.GetSides()( line_number));
        if( line_must_be_split[ line_number] == '0' || points_nearest_line( line_number).IsEmpty())
        {
          distance_limited_hull.PushBack( current_side.GetStartPoint());
          continue;
        }
        // add the starting and ending points to the nearest lines
        points_nearest_line( line_number).PushBack( current_side.GetStartPoint());
        points_nearest_line( line_number).PushBack( current_side.GetEndPoint());
        const storage::Vector< linal::Vector2D> &nearest_points( points_nearest_line( line_number));

        // determine the indices of the two endpoints on this hull line in the points_nearest_line vector
        const size_t start_point_id( nearest_points.GetSize() - 2);
        const size_t end_point_id( nearest_points.GetSize() - 1);

        // get the path with small side lengths that goes from one side of the convex hull on this side to the other
        // this is seeking to minimize the sum of areas given by rectangles drawn around each line
        graph::Path path;
        const size_t n_points( nearest_points.GetSize());
        graph::ConstGraph< size_t, double> distance_graph;
        {
          linal::Matrix< double> distances( n_points, n_points, undefined_edge_value);
          for( size_t point_a( 0); point_a < n_points; ++point_a)
          {
            for( size_t point_b( point_a + 1); point_b < n_points; ++point_b)
            {
              distances( point_a, point_b)
                = linal::SquareDistance( nearest_points( point_a), nearest_points( point_b));
            }
          }

          distance_graph = graph::ConstGraph< size_t, double>
                           (
                             storage::Vector< size_t>( n_points, size_t( 1)),
                             distances,
                             undefined_edge_value,
                             false,
                             true
                           );
        }
        path = graph::Connectivity::FindMinimalPath( distance_graph, start_point_id, end_point_id);

        if( path.GetSize() <= size_t( 2))
        {
          distance_limited_hull.PushBack( current_side.GetStartPoint());
          continue;
        }

        // create a polygon from the given path, skipping the last point which is the start point for the next path
        for( graph::Path::const_iterator itr( path.Begin()), itr_end( path.End() - 1); itr != itr_end; ++itr)
        {
          distance_limited_hull.PushBack( nearest_points( *itr));
        }
      }
      return distance_limited_hull;
    }

    //! @brief update the min and max bounds of the polygon with the current point
    //! @param POINT the point to update the min / max bounds with
    void Polygon::UpdateMinMaxBounds( const linal::Vector2D &NEW_POINT)
    {
      if( m_Sides.IsEmpty())
      {
        m_Bounds = LineSegment2D( NEW_POINT, NEW_POINT);
        return;
      }
      if( NEW_POINT.X() < m_Bounds.GetStartPoint().X())
      {
        m_Bounds.SetStartPoint( linal::Vector2D( NEW_POINT.X(), m_Bounds.GetStartPoint().Y()));
      }
      else if( NEW_POINT.X() > m_Bounds.GetEndPoint().X())
      {
        m_Bounds.SetEndPoint( linal::Vector2D( NEW_POINT.X(), m_Bounds.GetEndPoint().Y()));
      }
      if( NEW_POINT.Y() < m_Bounds.GetStartPoint().Y())
      {
        m_Bounds.SetStartPoint( linal::Vector2D( m_Bounds.GetStartPoint().X(), NEW_POINT.Y()));
      }
      else if( NEW_POINT.Y() > m_Bounds.GetEndPoint().Y())
      {
        m_Bounds.SetEndPoint( linal::Vector2D( m_Bounds.GetEndPoint().X(), NEW_POINT.Y()));
      }
    }

    //! @brief Expand the polygon; each vertex moves out by this much away from the barycenter
    //! @param EXPANSION amount to move each point away from the barycenter
    Polygon &Polygon::Expand( const double &EXPANSION)
    {
      const linal::Vector2D barycenter( GetCenter());
      for( auto itr( m_Sides.Begin()), itr_end( m_Sides.End()); itr != itr_end; ++itr)
      {
        const double prior_length( ( itr->GetStartPoint() - barycenter).Norm());
        const double expansion_ratio_begin( ( prior_length + EXPANSION) / prior_length);
        const double prior_length_end( ( itr->GetEndPoint() - barycenter).Norm());
        const double expansion_ratio_end( ( prior_length_end + EXPANSION) / prior_length_end);
        *itr = LineSegment2D
               (
                 ( itr->GetStartPoint() - barycenter) * expansion_ratio_begin + barycenter,
                 ( itr->GetEndPoint() - barycenter) * expansion_ratio_end + barycenter
               );
        UpdateMinMaxBounds( itr->GetStartPoint());
        UpdateMinMaxBounds( itr->GetEndPoint());
      }
      return *this;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Polygon::Read( std::istream &ISTREAM)
    {
      // read base classes
      io::Serialize::Read( m_Sides, ISTREAM);
      io::Serialize::Read( m_Bounds, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &Polygon::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base classes
      io::Serialize::Write( m_Sides, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bounds, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace coord
} // namespace bcl
