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

#ifndef BCL_COORD_LINE_SEGMENT_3D_H_
#define BCL_COORD_LINE_SEGMENT_3D_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LineSegment3D
    //! @brief A geometry class for a line segment.
    //! @details It has a start point, an end point and redundant the direction stored with it. You can access and change the
    //! members, dependencies are automatically updated. You can ask the line segment for its length.
    //! It's main use is to determine the shortest connection between two line segments.
    //!
    //! @see @link example_coord_line_segment_3d.cpp @endlink
    //! @author woetzen
    //! @date 11.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LineSegment3D :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      linal::Vector3D m_StartPoint; //!< start     point of line segment
      linal::Vector3D m_EndPoint;   //!< end       point of line segment
      linal::Vector3D m_Direction;  //!< direction vector from start to end

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LineSegment3D();

      //! @brief construct a LineSegment from start and end point
      //! @param START_POINT start point of line segment
      //! @param END_POINT end point of line segment
      LineSegment3D( const linal::Vector3D &START_POINT, const linal::Vector3D &END_POINT);

      //! copy constructor
      LineSegment3D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the start point of the line segment
      //! @return start point as Vector3D
      const linal::Vector3D &GetStartPoint() const;

      //! @brief returns the end point of the line segment
      //! @return end point as Vector3D
      const linal::Vector3D &GetEndPoint() const;

      //! @brief returns the direction line segment - vector from start to end point
      //! @return direction as Vector3D
      const linal::Vector3D &GetDirection() const;

      //! @brief set the start point
      //! @param START_POINT new starting point for this line segment
      void SetStartPoint( const linal::Vector3D &START_POINT);

      //! @brief set the end point
      //! @param END_POINT new ending point for this line segment
      void SetEndPoint( const linal::Vector3D &END_POINT);

      //! @brief Set Direction from StartPoint. Direction has to have the same length as the considered line segment
      //! @param DIRECTION the new direction, end point will be recalculated from start point and direction
      void SetDirection( const linal::Vector3D &DIRECTION);

    ////////////////
    // operations //
    ////////////////

      //! @brief returns length of LineSegment3D
      //! @return norm of direction
      double GetLength() const;

      //! @brief get the footprint fraction of a point onto this line
      //! @param POINT the point of interest
      //! @return the footprint: fraction along the line that the point is nearest
      double GetFootPointFraction( const linal::Vector3D &POINT) const;

      //! @brief calculates and returns the reverse LineSegment3D
      //! @return the reverse LineSegment3D
      LineSegment3D GetReverse() const;

      //! @brief shortens the line segment by the amount given equally from both ends
      //! @param LENGTH length to shorten the line segment by
      void Shorten( const double LENGTH);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write LineSegment3D to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read LineSegment3D from util::istream
      std::istream &Read( std::istream &ISTREAM);

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief recalculate direction
      //!        m_Direction = m_EndPoint - m_StartPoint : is called when start or end are set
      void RecalculateDirection();

    }; // class LineSegment3D

    //! @brief returns the shortest distance of two LineSegment3D LINESEGMENT_A and LINESEGMENT_B
    //! the pair contains the LineSegement connecting the two line segment by the shortest distance. the second value provides you with the\n
    //! information whether the connecting segment is orthogonal to the two line segments A and B and the line segments are not parallel
    //! @param LINESEGMENT_A first line segment
    //! @param LINESEGMENT_B second line segment
    //! @return returns the shortest connection as a line segment, and a bool to indicate, if it is an orthogonal connection
    BCL_API
    storage::Pair< LineSegment3D, bool>
    ShortestConnectionBetweenLineSegments3D
    (
      const LineSegment3D &LINESEGMENT_A,
      const LineSegment3D &LINESEGMENT_B
    );

    //! @brief returns distance of point POINT from LineSegment3D LINESEGMENT
    //! @param LINESEGMENT from which distance to point is to be calculated
    //! @param POINT point from which distance to LineSegment3D is to be calculated
    //! @return returns the distance as a double, and a bool to indicate, if it is an orthogonal connection
    BCL_API
    storage::Pair< double, bool>
    CalculateDistancePointFromLineSegment
    (
      const LineSegment3D &LINESEGMENT,
      const linal::Vector3D &POINT
    );

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_LINE_SEGMENT_3D_H_
