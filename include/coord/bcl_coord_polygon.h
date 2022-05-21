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

#ifndef BCL_COORD_POLYGON_H_
#define BCL_COORD_POLYGON_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_coord_line_segment_2d.h"
#include "bcl_coord_movable_interface.h"
#include "linal/bcl_linal_vector_2d.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Polygon
    //! @brief models a polygon as a series of connected line segments
    //!
    //! @see @link example_coord_polygon.cpp @endlink
    //! @author mendenjl
    //! @date Dec 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Polygon :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::Vector< LineSegment2D> m_Sides; //!< Actual sides of the polygon
      LineSegment2D m_Bounds; //!< min/max x/y for the current polygon

    public:

      typedef storage::Vector< LineSegment2D>::iterator iterator;
      typedef storage::Vector< LineSegment2D>::const_iterator const_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      Polygon();

      //! construct from points
      Polygon( const storage::Vector< linal::Vector2D> &POINTS);

      //! copy constructor
      Polygon *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns number of sides of the polygon
      //! @return number of sides of the polygon
      size_t GetNumberOfSides() const;

      //! @brief returns number of sides of the polygon
      //! @return number of sides of the polygon
      double GetPerimeter() const;

      //! @brief get the area of the polygon
      //! @return area of the polygon
      double GetArea() const;

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin();

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End();

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const;

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const;

      //! @brief add a line between the last point in the polygon and the first point
      //! @param POINT the point to add
      void PushBack( const linal::Vector2D &POINT);

      //! @brief conversion to a vector of line segments
      operator const storage::Vector< LineSegment2D> &() const;

      //! @brief conversion to a vector of line segments
      operator storage::Vector< LineSegment2D> &();

      //! @brief Get the underlying sides of the polygon
      //! @return the sides of the polygon
      const storage::Vector< LineSegment2D> &GetSides() const;

      //! @brief Get the underlying sides of the polygon
      //! @return the sides of the polygon
      storage::Vector< LineSegment2D> &GetSides();

    ////////////////
    // operations //
    ////////////////

      //! returns the geometric center of the object
      linal::Vector2D GetCenter() const;

      //! returns the barycenter of the object
      linal::Vector2D GetBarycenter() const;

      //! @brief test whether a particular point is a corner (intersection of two sides) of this polygon
      //! @param POINT the point to test for being a corner in this polygon
      bool IsCornerOf( const linal::Vector2D &POINT) const;

      //! @brief test whether a particular point is within the polygon
      //! @param POINT the point to test for inclusion in the polygon
      bool IsWithin( const linal::Vector2D &POINT) const;

      //! @brief Find the side nearest to a point, which need not be inside the polygon
      //! @param POINT the point to test for inclusion in the polygon
      //! @return an iterator to the nearest side and the actual distance
      std::pair< Polygon::const_iterator, double> FindNearestSide( const linal::Vector2D &POINT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief create a polygon from the convex hull, which is the minimal convex polygon that can enclose all given points
      //! @param POINTS the points to enclose with the convex hull
      //! @param MAX_SIDE_LENGTH maximum length of a given side, if the length would be longer, a concavity is introduced into the polygon
      //! @param POINT_RADIUS radius of each of the points
      //! @return the minimal convex hull
      static Polygon ConvexHull
      (
        const storage::Vector< linal::Vector2D> &POINTS,
        const double &MAX_SIDE_LENGTH = std::numeric_limits< double>::infinity(),
        const double &RADIUS = 0.0
      );

      //! @brief Expand the polygon; each vertex moves out by this much away from the barycenter
      //! @param EXPANSION amount to move each point away from the barycenter
      Polygon &Expand( const double &EXPANSION);

    protected:

      //! @brief update the min and max bounds of the polygon with the current point
      //! @param POINT the point to update the min / max bounds with
      void UpdateMinMaxBounds( const linal::Vector2D &NEW_POINT);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Polygon

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_POLYGON_H_
