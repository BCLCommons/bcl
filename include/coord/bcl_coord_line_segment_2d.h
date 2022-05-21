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

#ifndef BCL_COORD_LINE_SEGMENT_2D_H_
#define BCL_COORD_LINE_SEGMENT_2D_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_2d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LineSegment2D
    //! @brief A geometry class for a line segment.
    //! @details It has a start point, an end point and redundant the direction stored with it. You can access and change the
    //! members, dependencies are automatically updated. You can ask the line segment for its length.
    //! It's main use is to determine the shortest connection between two line segments.
    //!
    //! @see @link example_coord_line_segment_2d.cpp @endlink
    //! @author mendenjl
    //! @date Dec 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LineSegment2D :
      public util::ObjectInterface
    {

    //////////
    // data //
    //////////

      linal::Vector2D m_StartPoint; //!< start     point of line segment
      linal::Vector2D m_EndPoint;   //!< end       point of line segment
      linal::Vector2D m_Direction;  //!< direction vector from start to end

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
      LineSegment2D();

      //! @brief construct a LineSegment from start and end point
      //! @param START_POINT start point of line segment
      //! @param END_POINT end point of line segment
      LineSegment2D( const linal::Vector2D &START_POINT, const linal::Vector2D &END_POINT);

      //! copy constructor
      LineSegment2D *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the start point of the line segment
      //! @return start point as Vector2D
      const linal::Vector2D &GetStartPoint() const;

      //! @brief returns the end point of the line segment
      //! @return end point as Vector2D
      const linal::Vector2D &GetEndPoint() const;

      //! @brief returns the direction line segment - vector from start to end point
      //! @return direction as Vector2D
      const linal::Vector2D &GetDirection() const;

      //! @brief set the start point
      //! @param START_POINT new starting point for this line segment
      void SetStartPoint( const linal::Vector2D &START_POINT);

      //! @brief set the end point
      //! @param END_POINT new ending point for this line segment
      void SetEndPoint( const linal::Vector2D &END_POINT);

      //! @brief Set Direction from StartPoint. Direction has to have the same length as the considered line segment
      //! @param DIRECTION the new direction, end point will be recalculated from start point and direction
      void SetDirection( const linal::Vector2D &DIRECTION);

    ////////////////
    // operations //
    ////////////////

      //! @brief returns length of LineSegment2D
      //! @return norm of direction
      double GetLength() const;

      //! @brief get the footprint fraction of a point onto this line
      //! @param POINT the point of interest
      //! @return the footprint: fraction along the line that the point is nearest
      double GetFootPointFraction( const linal::Vector2D &POINT) const;

      //! @brief Compute the distance from this line to a given point
      //! @param POINT the point of interest
      //! @return distance from this line segment to the point
      double DistanceToPoint( const linal::Vector2D &POINT) const;

      //! @brief Compute the distance from this line to a given point
      //! @param LINE line of interest
      //! @return returns the shortest connection as a line segment, and a bool to indicate, if it is an orthogonal connection
      storage::Pair< LineSegment2D, bool> ShortestLineBetween( const LineSegment2D &LINE) const;

      //! @brief test whether two line segments intersect
      //! @param LINE line of interest
      //! @return true if the lines intersect
      bool DoesIntersect( const LineSegment2D &LINE) const;

      //! @brief test whether two line segments overlap, that is, they intersect more than once
      //! @param LINE line of interest
      //! @return true if the lines overlap
      bool Overlaps( const LineSegment2D &LINE) const;

      //! @brief calculates and returns the reverse LineSegment2D
      //! @return the reverse LineSegment2D
      LineSegment2D GetReverse() const;

      //! @brief shortens the line segment by the amount given equally from both ends
      //! @param LENGTH length to shorten the line segment by
      void Shorten( const double LENGTH);

      //! @brief equality test
      //! @param LINE other line to test
      bool operator ==( const LineSegment2D &LINE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! write LineSegment2D to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read LineSegment2D from util::istream
      std::istream &Read( std::istream &ISTREAM);

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief recalculate direction
      //!        m_Direction = m_EndPoint - m_StartPoint : is called when start or end are set
      void RecalculateDirection();

    }; // class LineSegment2D

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_LINE_SEGMENT_2D_H_
