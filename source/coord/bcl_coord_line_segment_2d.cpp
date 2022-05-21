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
#include "coord/bcl_coord_line_segment_2d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_2d_operations.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LineSegment2D::s_Instance
    (
      GetObjectInstances().AddInstance( new LineSegment2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructer
    LineSegment2D::LineSegment2D()
    {
    }

    //! @brief construct a LineSegment from start and end point
    //! @param START_POINT start point of line segment
    //! @param END_POINT end point of line segment
    LineSegment2D::LineSegment2D( const linal::Vector2D &START_POINT, const linal::Vector2D &END_POINT) :
      m_StartPoint( START_POINT),
      m_EndPoint( END_POINT),
      m_Direction( END_POINT - START_POINT)
    {
    }

    //! copy constructor
    LineSegment2D *LineSegment2D::Clone() const
    {
      return new LineSegment2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LineSegment2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! returns the start point of the line segment
    const linal::Vector2D &LineSegment2D::GetStartPoint() const
    {
      return m_StartPoint;
    }

    //! returns the end point of the line segment
    const linal::Vector2D &LineSegment2D::GetEndPoint() const
    {
      return m_EndPoint;
    }

    //! returns the direction line segment - vector from start to end point
    const linal::Vector2D &LineSegment2D::GetDirection() const
    {
      return m_Direction;
    }

    //! set the start point
    void LineSegment2D::SetStartPoint( const linal::Vector2D &START_POINT)
    {
      m_StartPoint = START_POINT;
      RecalculateDirection();
    }

    //! set the end point
    void LineSegment2D::SetEndPoint( const linal::Vector2D &END_POINT)
    {
      m_EndPoint = END_POINT;
      RecalculateDirection();
    }

    //! Set Direction from StartPoint. Direction has to have the same length as the considered linesegment
    void LineSegment2D::SetDirection( const linal::Vector2D &DIRECTION)
    {
      m_Direction = DIRECTION;
      m_EndPoint  = m_StartPoint + m_Direction;
    }

  ////////////////
  // operations //
  ////////////////

    //! returns length of LineSegment2D
    double LineSegment2D::GetLength() const
    {
      return m_Direction.Norm();
    }

    //! @brief get the footprint fraction of a point onto this line
    //! @param POINT the point of interest
    //! @return the footprint: fraction along the line that the point is nearest
    double LineSegment2D::GetFootPointFraction( const linal::Vector2D &POINT) const
    {
      return linal::ScalarProduct( ( POINT - m_StartPoint), m_Direction) / m_Direction.SquareNorm();
    }

    //! @brief Compute the distance from this line to a given point
    //! @param POINT the point of interest
    //! @return distance from this line segment to the point
    double LineSegment2D::DistanceToPoint( const linal::Vector2D &POINT) const
    {
      // calculate footpoint position in terms of fractional length of LINE vector
      const double footpoint_fraction( GetFootPointFraction( POINT));

      // if footpoint_fraction is smaller than 0, then POINT lies beyond LineSegment2D on start point side
      if( footpoint_fraction <= 0.0)
      {
        // return distance between POINT and m_StartPoint
        return linal::Distance( m_StartPoint, POINT);
      }
      // if footpoint_fraction is larger than 1, then POINT lies beyond LineSegment2D on end point side
      if( footpoint_fraction >= 1.0)
      {
        // return distance between POINT and m_EndPoint
        return linal::Distance( m_EndPoint, POINT);
      }

      // return distance between POINT and footpoint
      return linal::Distance( m_StartPoint + footpoint_fraction * m_Direction, POINT);
    }

    //! @brief Compute the distance from this line to a given point
    //! @param LINE line of interest
    //! @return shortest line segment between the two lines (length == 0 if they intersect)
    storage::Pair< LineSegment2D, bool> LineSegment2D::ShortestLineBetween( const LineSegment2D &LINE) const
    {
      // http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
      // this code is identical to the code for 3d vectors
      // static number that is relatively small to be used in finding if two lines are parallel or not
      static const double s_small_number( 0.00001);

      // create references of directions of both linesegments
      const linal::Vector2D &direction_a( m_Direction);
      const linal::Vector2D &direction_b( LINE.m_Direction);

      const linal::Vector2D w( m_StartPoint - LINE.m_StartPoint);
      const double         a( direction_a * direction_a);           // always >= 0
      const double         b( direction_a * direction_b);
      const double         c( direction_b * direction_b);           // always >= 0
      const double         d( direction_a * w);
      const double         e( direction_b * w);
      const double         denominator( a * c - b * b);   // always >= 0
      double               sc, sN, sD = denominator;      // sc = sN / sD, default sD = D >= 0
      double               tc, tN, tD = denominator;      // tc = tN / tD, default tD = D >= 0

      bool footpoints_on_segments( true);

      // compute the line parameters of the two closest points
      if( denominator < s_small_number) // the lines are almost parallel
      {
        sN = 0.0;          // force using point P0 on segment S1
        sD = 1.0;          // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
      }
      else                 // get the closest points on the infinite lines
      {
        sN = ( b * e - c * d);
        tN = ( a * e - b * d);
        if( sN < 0.0)      // sc < 0 => the s=0 edge is visible
        {
          sN = 0.0;
          tN = e;
          tD = c;
          footpoints_on_segments = false;
        }
        else if( sN > sD) // sc > 1 => the s=1 edge is visible
        {
          sN = sD;
          tN = e + b;
          tD = c;
          footpoints_on_segments = false;
        }
      }

      if( tN < 0.0)        // tc < 0 => the t=0 edge is visible
      {
        tN = 0.0;
        // recompute sc for this edge
        if( -d < 0.0)
        {
          sN = 0.0;
        }
        else if( -d > a)
        {
          sN = sD;
        }
        else
        {
          sN = -d;
          sD = a;
        }
        footpoints_on_segments = false;
      }
      else if( tN > tD)    // tc > 1 => the t=1 edge is visible
      {
        tN = tD;
        // recompute sc for this edge
        if( ( -d + b) < 0.0)
        {
          sN = 0;
        }
        else if( ( -d + b) > a)
        {
          sN = sD;
        }
        else
        {
          sN = ( -d + b);
          sD = a;
        }
        footpoints_on_segments = false;
      }

      // finally do the division to get sc and tc
      sc = ( math::Absolute( sN) < s_small_number ? 0.0 : sN / sD);
      tc = ( math::Absolute( tN) < s_small_number ? 0.0 : tN / tD);

      // get the difference of the two closest points
      const linal::Vector2D foot_point_1( m_StartPoint + sc * direction_a);
      const linal::Vector2D foot_point_2( LINE.m_StartPoint + tc * direction_b);

      // return the connection line and whether the connecting linesegment is both sided orthogonal to the given linesegments
      return storage::Pair< LineSegment2D, bool>( LineSegment2D( foot_point_1, foot_point_2), footpoints_on_segments);
    }

    //! @brief test whether two line segments intersect
    //! @param LINE line of interest
    //! @return true if the lines intersect
    bool LineSegment2D::DoesIntersect( const LineSegment2D &LINE) const
    {
      const double denom( linal::CrossProduct( m_Direction, LINE.m_Direction));
      if( denom == 0.0)
      {
        // lines are collinear, test for any overlap.  If overlap is present, the lines definitely intersect
        return
           (( LINE.m_StartPoint.X() < m_StartPoint.X()) != ( LINE.m_StartPoint.X() < m_EndPoint.X()))
        && (( LINE.m_StartPoint.Y() < m_StartPoint.Y()) != ( LINE.m_StartPoint.Y() < m_EndPoint.Y())); // Collinear && overlaps
      }
      const bool denom_positive( denom > 0);

      linal::Vector2D difference( m_StartPoint);
      difference -= LINE.m_StartPoint;
      double s_numer( linal::CrossProduct( m_Direction, difference));
      if( ( s_numer < 0) == denom_positive)
      {
        return false; // No collision
      }

      double t_numer( linal::CrossProduct( LINE.m_Direction, difference));
      return    ( t_numer < 0)     != denom_positive
             && ( s_numer > denom) != denom_positive
             && ( t_numer > denom) != denom_positive;
    }

    //! @brief test whether two line segments overlap, that is, they intersect more than once
    //! @param LINE line of interest
    //! @return true if the lines overlap
    bool LineSegment2D::Overlaps( const LineSegment2D &LINE) const
    {
      if( linal::CrossProduct( m_Direction, LINE.m_Direction) == double( 0.0) && GetLength() > 0.0 && LINE.GetLength() > 0.0)
      {
        // lines are collinear, test for any overlap.  If overlap is present, the lines definitely intersect
        // test for whether X or Y axes overlap
        if
        (
             (( LINE.m_StartPoint.X() < m_StartPoint.X()) != ( LINE.m_StartPoint.X() < m_EndPoint.X()))
          && (( LINE.m_StartPoint.Y() < m_StartPoint.Y()) != ( LINE.m_StartPoint.Y() < m_EndPoint.Y()))
        )
        {
          if
          (
            m_StartPoint != LINE.m_StartPoint
            && m_StartPoint != LINE.m_EndPoint
            && m_EndPoint != LINE.m_StartPoint
            && m_EndPoint != LINE.m_EndPoint
          )
          {
            // no boundary points are identical, so there is non-trivial overlap
            return true;
          }
          // test for identical lines
          if( operator ==( LINE))
          {
            return true;
          }
          if( m_StartPoint == LINE.m_EndPoint && m_EndPoint == LINE.m_StartPoint)
          {
            // reversed line
            return true;
          }
          // exactly one boundary point in common; overlap exists provided that the directions are the same
          return m_Direction * LINE.m_Direction > 0;
        }
      }
      return false;
    }

    //! @brief calculates and returns the reverse LineSegment2D
    //! @return the reverse LineSegment2D
    LineSegment2D LineSegment2D::GetReverse() const
    {
      return LineSegment2D( m_EndPoint, m_StartPoint);
    }

    //! @brief shortens the line segment by the amount given equally from both ends
    //! @param LENGTH length to shorten the line segment by
    void LineSegment2D::Shorten( const double LENGTH)
    {
      // normalize a copy of the direction vector
      linal::Vector2D normalized_direction( m_Direction);
      normalized_direction.Normalize();

      // calculate the vector from the center
      const linal::Vector2D vector_from_center
      (
        std::max( 0.0, ( GetLength() - LENGTH) / 2.0) * normalized_direction
      );

      // calculate the center
      const linal::Vector2D center( m_StartPoint + m_Direction / 2.0);

      // update the members
      m_StartPoint = center - vector_from_center;
      m_EndPoint = center + vector_from_center;
      RecalculateDirection();
    }

    //! @brief equality test
    //! @param LINE other line to test
    bool LineSegment2D::operator ==( const LineSegment2D &LINE) const
    {
      return m_StartPoint == LINE.m_StartPoint && m_EndPoint == LINE.m_EndPoint;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write LineSegment2D to std::ostream
    std::ostream &LineSegment2D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write data
      io::Serialize::Write( m_StartPoint, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EndPoint, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read AA from io::IFStream
    std::istream &LineSegment2D::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_StartPoint, ISTREAM);
      io::Serialize::Read( m_EndPoint, ISTREAM);

      //recalculate direction
      RecalculateDirection();

      //return
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! recalculate m_Direction = m_EndPoint - m_StartPoint
    void LineSegment2D::RecalculateDirection()
    {
      m_Direction = m_EndPoint - m_StartPoint;
    }

  } // namespace coord
} // namespace bcl
