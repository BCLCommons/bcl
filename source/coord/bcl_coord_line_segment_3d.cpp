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
#include "coord/bcl_coord_line_segment_3d.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
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
    const util::SiPtr< const util::ObjectInterface> LineSegment3D::s_Instance
    (
      GetObjectInstances().AddInstance( new LineSegment3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructer
    LineSegment3D::LineSegment3D()
    {
    }

    //! @brief construct a LineSegment from start and end point
    //! @param START_POINT start point of line segment
    //! @param END_POINT end point of line segment
    LineSegment3D::LineSegment3D( const linal::Vector3D &START_POINT, const linal::Vector3D &END_POINT) :
      m_StartPoint( START_POINT),
      m_EndPoint( END_POINT),
      m_Direction( END_POINT - START_POINT)
    {
    }

    //! copy constructor
    LineSegment3D *LineSegment3D::Clone() const
    {
      return new LineSegment3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LineSegment3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! returns the start point of the line segment
    const linal::Vector3D &LineSegment3D::GetStartPoint() const
    {
      return m_StartPoint;
    }

    //! returns the end point of the line segment
    const linal::Vector3D &LineSegment3D::GetEndPoint() const
    {
      return m_EndPoint;
    }

    //! returns the direction line segment - vector from start to end point
    const linal::Vector3D &LineSegment3D::GetDirection() const
    {
      return m_Direction;
    }

    //! set the start point
    void LineSegment3D::SetStartPoint( const linal::Vector3D &START_POINT)
    {
      m_StartPoint = START_POINT;
      RecalculateDirection();
    }

    //! set the end point
    void LineSegment3D::SetEndPoint( const linal::Vector3D &END_POINT)
    {
      m_EndPoint = END_POINT;
      RecalculateDirection();
    }

    //! Set Direction from StartPoint. Direction has to have the same length as the considered linesegment
    void LineSegment3D::SetDirection( const linal::Vector3D &DIRECTION)
    {
      m_Direction = DIRECTION;
      m_EndPoint  = m_StartPoint + m_Direction;
    }

  ////////////////
  // operations //
  ////////////////

    //! returns length of LineSegment3D
    double LineSegment3D::GetLength() const
    {
      return m_Direction.Norm();
    }

    //! @brief get the footprint fraction of a point onto this line
    //! @param POINT the point of interest
    //! @return the footprint: fraction along the line that the point is nearest
    double LineSegment3D::GetFootPointFraction( const linal::Vector3D &POINT) const
    {
      return linal::ScalarProduct( ( POINT - m_StartPoint), m_Direction) / m_Direction.SquareNorm();
    }

    //! @brief calculates and returns the reverse LineSegment3D
    //! @return the reverse LineSegment3D
    LineSegment3D LineSegment3D::GetReverse() const
    {
      return LineSegment3D( m_EndPoint, m_StartPoint);
    }

    //! @brief shortens the line segment by the amount given equally from both ends
    //! @param LENGTH length to shorten the line segment by
    void LineSegment3D::Shorten( const double LENGTH)
    {
      // normalize a copy of the direction vector
      linal::Vector3D normalized_direction( m_Direction);
      normalized_direction.Normalize();

      // calculate the vector from the center
      const linal::Vector3D vector_from_center
      (
        std::max( 0.0, ( GetLength() - LENGTH) / 2.0) * normalized_direction
      );

      // calculate the center
      const linal::Vector3D center( m_StartPoint + m_Direction / 2.0);

      // update the members
      m_StartPoint = center - vector_from_center;
      m_EndPoint = center + vector_from_center;
      RecalculateDirection();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write LineSegment3D to std::ostream
    std::ostream &LineSegment3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write data
      io::Serialize::Write( m_StartPoint, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EndPoint, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read AA from io::IFStream
    std::istream &LineSegment3D::Read( std::istream &ISTREAM)
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
    void LineSegment3D::RecalculateDirection()
    {
      m_Direction = m_EndPoint - m_StartPoint;
    }

    //! returns the shortest distance of two LineSegment3D LINESEGMENT_A and LINESEGMENT_B
    //! http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
    //! the pair contains the LineSegement connecting the two linesegment by the shortest distance. the second value provides you with the\n
    //! information whether the connecting segment is orthogonal to the two linesegments A and B and the linesegments are not parallel
    storage::Pair< LineSegment3D, bool>
    ShortestConnectionBetweenLineSegments3D
    (
      const LineSegment3D &LINESEGMENT_A,
      const LineSegment3D &LINESEGMENT_B
    )
    {
      // This code is slightly adapted from http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm, and hence
      // subject to
      // Copyright 2001 softSurfer, 2012 Dan Sunday
      // This code may be freely used and modified for any purpose
      // providing that this copyright notice is included with it.
      // SoftSurfer makes no warranty for this code, and cannot be held
      // liable for any real or imagined damage resulting from its use.
      // Users of this code must verify correctness for their application.

      // static number that is relatively small to be used in finding if two lines are parallel or not
      static const double s_small_number( 0.00001);

      // create references of directions of both linesegments
      const linal::Vector3D &direction_a( LINESEGMENT_A.GetDirection());
      const linal::Vector3D &direction_b( LINESEGMENT_B.GetDirection());

      const linal::Vector3D w( LINESEGMENT_A.GetStartPoint() - LINESEGMENT_B.GetStartPoint());
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
      const linal::Vector3D foot_point_1( LINESEGMENT_A.GetStartPoint() + sc * direction_a);
      const linal::Vector3D foot_point_2( LINESEGMENT_B.GetStartPoint() + tc * direction_b);

      // return the connection line and whether the connecting linesegment is both sided orthogonal to the given linesegments
      return storage::Pair< LineSegment3D, bool>( LineSegment3D( foot_point_1, foot_point_2), footpoints_on_segments);
    }

    //! @brief returns distance of point POINT from LineSegment3D LINESEGMENT
    //! @param LINESEGMENT from which distance to point is to be calculated
    //! @param POINT point from which distance to LineSegment3D is to be calculated
    //! @return returns the distance as a double, and a bool to indicate, if it is an orthogonal connection
    storage::Pair< double, bool>
    CalculateDistancePointFromLineSegment
    (
      const LineSegment3D &LINESEGMENT,
      const linal::Vector3D &POINT
    )
    {
      // calculate footpoint position in terms of fractional length of LINE vector
      const double footpoint_fraction( LINESEGMENT.GetFootPointFraction( POINT));

      // if footpoint_fraction is smaller than 0, then POINT lies beyond LineSegment3D on start point side
      if( footpoint_fraction <= 0.0)
      {
        // return distance between POINT and m_StartPoint
        return storage::Pair< double, bool>( linal::Distance( LINESEGMENT.GetStartPoint(), POINT), false);
      }
      // if footpoint_fraction is larger than 1, then POINT lies beyond LineSegment3D on end point side
      if( footpoint_fraction >= 1.0)
      {
        // return distance between POINT and m_EndPoint
        return storage::Pair< double, bool>( linal::Distance( LINESEGMENT.GetEndPoint(), POINT), false);
      }

      // in case the projection of POINT onto LineSegment3D falls on the LineSegment3D, calculate actual footpoint
      const linal::Vector3D footpoint( LINESEGMENT.GetStartPoint() + footpoint_fraction * LINESEGMENT.GetDirection());

      // return distance between POINT and footpoint
      return storage::Pair< double, bool>( linal::Distance( footpoint, POINT), true);
    }

  } // namespace coord
} // namespace bcl
