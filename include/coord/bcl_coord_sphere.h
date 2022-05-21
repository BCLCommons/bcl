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

#ifndef BCL_COORD_SPHERE_H_
#define BCL_COORD_SPHERE_H_

// include the namespace header
#include "bcl_coord.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d.h"
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Sphere
    //! @brief is a class representing a 3d sphere
    //! @details It defines the position in space and the radius of the Sphere.
    //! It provides basic sphere arithmetic - volume, surface area, free surface area in respect to a set of other
    //! spheres, surface interpolation by point representation, algebraic solutions for volume and surface overlap
    //!
    //! @see @link example_coord_sphere.cpp @endlink
    //! @author staritrd, woetzen
    //! @date 26.09.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Sphere :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector3D m_Position; //!< Position
      double          m_Radius;   //!< Radius

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! The default constructor, initializes all values as zero.
      Sphere();

      //! Single value constructor.
      Sphere( const linal::Vector3D &POSITION, const double RADIUS);

      //! Virtual copy constructor
      Sphere *Clone() const
      {
        return new Sphere( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! Access to radius.
      void SetRadius( const double RADIUS)
      {
        BCL_Assert( RADIUS >= 0.0, "radius needs to be >= 0.0");
        m_Radius = RADIUS;
      }

      //! Read only access to radius.
      double GetRadius() const
      {
        return m_Radius;
      }

      //! Access to position
      void SetPosition( const linal::Vector3D &POSITION)
      {
        BCL_Assert( POSITION.IsDefined(), "position is not defined");
        m_Position = POSITION;
      }

      //! Read only access to position.
      const linal::Vector3D &GetPosition() const
      {
        return m_Position;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Calculates the total surface of the sphere.
      //! @return surface = 4 * pi * r^2
      double GetSurface() const
      {
        return SurfaceArea( m_Radius);
      }

      //! @brief Calculates the volume of the sphere.
      //! @return volume = 4/3 * pi * r^3
      double GetVolume() const
      {
        return Volume( m_Radius);
      }

      //! @brief check if the given sphere does overlap with this sphere
      //! @param SPHERE the sphere to be considered for the overlap
      //! @return true, if the sum of the radii is larger than the euclidean distance between the positions
      bool DoesOverlap( const Sphere &SPHERE) const;

      //! @brief check if given point is contained in the sphere
      //! @param POINT a point to be considered
      //! @return true, if the distance point->position is larger than the radius
      bool DoesContain( const linal::Vector3D &POINT) const;

      //! @brief distance and volume overlap
      //! @param SPHERE the sphere that overlaps with this sphere
      //! @return pair of distance and bool - true if volume overlap > 0 exists, false otherwise
      storage::Pair< double, bool> DistanceAndVolumeOverlap( const Sphere &SPHERE) const;

      //! @brief determine the surface overlap of two spheres analytically
      //! http://mathworld.wolfram.com/Sphere-SphereIntersection.html
      //! @param SPHERE the sphere that overlaps with this sphere
      //! @return the surface overlap of two spheres - 0 if they do not overlap, or one sphere is smaller than the other
      //!         and the surfaces do not touch - first is surface of this sphere, second of argument SPHERE
      storage::VectorND< 2, double> SurfaceOverlap( const Sphere &SPHERE) const;

      //! @brief Volume overlap of two spheres
      //! http://mathworld.wolfram.com/Sphere-SphereIntersection.html
      //! @param SPHERE the sphere that overlaps with this sphere
      //! @return the volume of the overlapping lens of this sphere with argument SPHERE - 0 if no overlap
      double VolumeOverlap( const Sphere &SPHERE) const;

      //! @brief represent surface by points using latitude sampling
      //! @param NUMBER number of points - the actual number of points can slightly deviate
      //! @return vector of points, representing the surface
      storage::Vector< linal::Vector3D> PointsOnSurfaceByLatiudes( const size_t NUMBER) const;

      //! @brief represent surface by points randomly place on surface
      //! @param NUMBER number of points
      //! @return vector of points, representing the surface
      storage::Vector< linal::Vector3D> PointsOnSurfaceRandom( const size_t NUMBER) const;

      //! @brief represents surface by points using a spiral sampling
      //! http://sitemason.vanderbilt.edu/page/hmbADS#spiral
      //! @param NUMBER number of points - the actual number of points can slightly deviate
      //! @return vector of points, reprensenting the surface
      storage::Vector< linal::Vector3D> PointsOnSurfaceSpiral( const size_t NUMBER) const;

      //! @brief calculate the free surface area fraction, that is not included by any given overlapping sphere
      //! @param SPHERES a list of spheres to be considered
      //! @param NUMBER the number of points that are used to interpolate the surface of the sphere - the higher, the more accurate
      //! @return fraction [0,1] where 0 is no free surface area, 1 is no overlap with any sphere
      double FreeSurfaceAreaFraction( const util::SiPtrList< const Sphere> &SPHERES, const size_t NUMBER) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! Writes Sphere to ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read sphere from istream
      std::istream &Read( std::istream &ISTREAM);

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief derive cap heights of overlapping spheres
      //! @param DISTANCE the distance of the two caps
      //! @param SPHERE the other sphere
      //! @return pair of heights of caps - first for this cap, second for the other
      storage::VectorND< 2, double> HeightOfCaps( const double DISTANCE, const Sphere &SPHERE) const;

    public:

      //! @brief Surface area of a sphere given a radius
      //! @param RADIUS radius of the sphere
      //! @return the surface area of the sphere; surface = 4 * pi * r^2
      static double SurfaceArea( const double RADIUS)
      {
        return double( 4) * math::g_Pi * RADIUS * RADIUS;
      }

      //! @brief Volume of sphere given a radius
      //! @param RADIUS radius of the sphere
      //! @return the volume of the sphere
      static double Volume( const double RADIUS)
      {
        return double( 4) / double( 3) * math::g_Pi * RADIUS * RADIUS * RADIUS;
      }

      //! @brief radius of sphere from volume
      //! @param VOLUME the volume of the sphere
      //! @return the radius fo the sphere
      static double Radius( const double VOLUME)
      {
        return math::Pow( VOLUME * double( 3) / ( 4 * math::g_Pi), double( 1) / double( 3));
      }

    }; // class Sphere

    //! @brief list overlapping spheres for each sphere
    //! @param SPHERES all spheres that are considered
    //! @return a map, that has a list of overlapping spheres (data) for each sphere (key)
    BCL_API
    storage::Map< util::SiPtr< const Sphere>, util::SiPtrList< const Sphere> >
    ListOverlaps( const util::SiPtrVector< const Sphere> &SPHERES);

    //! @brief calculate the solvent accessible surface area for each given sphere
    //! this is achieved by interpolating the area of a sphere by equally distributing NUMBER of points on the surface
    //! of the sphere and calculating the ration of points that are within the radius of any other sphere to NUMBER
    //! @see Shrake A, Rupley JA. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 79(2):351-71.
    //! @param SPHERES the spheres to be considered
    //! @param NUMBER number of points to interpolate the surface of the sphere - the larger the more accurate, but slower
    //! @return a map that has an absolute SASA for each sphere
    BCL_API
    storage::Map< util::SiPtr< const Sphere>, double>
    FreeSurface
    (
      const util::SiPtrVector< const Sphere> &SPHERES,
      const size_t NUMBER
    );

    //! @brief calculate the solvent accessible surface area ratio for each given sphere
    //! this is achieved by interpolating the area of a sphere by equally distributing NUMBER of points on the surface
    //! of the sphere and calculating the ration of points that are within the radius of any other sphere to NUMBER
    //! @param SPHERES the spheres to be considered
    //! @param NUMBER number of points to interpolate the surface of the sphere - the larger the more accurate, but slower
    //! @return a map that has a relative area for each sphere
    BCL_API
    storage::Map< util::SiPtr< const Sphere>, double>
    FreeSurfaceRatio
    (
      const util::SiPtrVector< const Sphere> &SPHERES,
      const size_t NUMBER
    );

  } // namespace coord
} // namespace bcl

#endif //BCL_COORD_SPHERE_H_
