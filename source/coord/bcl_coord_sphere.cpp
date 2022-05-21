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
#include "coord/bcl_coord_sphere.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Sphere::s_Instance
    (
      GetObjectInstances().AddInstance( new Sphere())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! The default constructor, initializes all values as zero.
    Sphere::Sphere() :
      m_Position( 0, 0, 0),
      m_Radius( 0)
    {
    }

    //! Single value constructor.
    Sphere::Sphere( const linal::Vector3D &POSITION, const double RADIUS) :
      m_Position( POSITION),
      m_Radius( RADIUS)
    {
      BCL_Assert( RADIUS >= 0, "radius needs to be >= 0.0");
      BCL_Assert( POSITION.IsDefined(), "position is not defined");
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Sphere::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief check if the given sphere does overlap with this sphere
    //! @param SPHERE the sphere to be considered for the overlap
    //! @return true, if the sum of the radii is larger than the euclidean distance between the positions
    bool Sphere::DoesOverlap( const Sphere &SPHERE) const
    {
      const linal::Vector3D connection( m_Position - SPHERE.m_Position);
      const double max_distance( m_Radius + SPHERE.m_Radius);

      // if the length of the largest of x, y or z distance is larger that the max distance, the spheres do not overlap
      if
      (
           math::Absolute( connection.X()) >= max_distance
        || math::Absolute( connection.Y()) >= max_distance
        || math::Absolute( connection.Z()) >= max_distance)
      {
        return false;
      }

      // compare the euclidean distance
      if( connection.SquareNorm() >= math::Sqr( max_distance))
      {
        return false;
      }

      // both spheres overlap
      return true;
    }

    //! @brief check if given point is contained in the sphere
    //! @param POINT a point to be considered
    //! @return true, if the distance point->position is larger than the radius
    bool Sphere::DoesContain( const linal::Vector3D &POINT) const
    {
      const linal::Vector3D connection( m_Position - POINT);

      // if the length of the largest of x, y or z distance is larger that the RADIUS, the point is not within
      if
      (
           math::Absolute( connection.X()) >= m_Radius
        || math::Absolute( connection.Y()) >= m_Radius
        || math::Absolute( connection.Z()) >= m_Radius
      )
      {
        return false;
      }

      // compare the euclidean distance
      if( connection.SquareNorm() >= math::Sqr( m_Radius))
      {
        return false;
      }

      // point is within sphere
      return true;
    }

    //! @brief distance and volume overlap
    //! @param SPHERE the sphere that overlaps with this sphere
    //! @return pair of distance and bool - true if volume overlap > 0 exists, false otherwise
    storage::Pair< double, bool> Sphere::DistanceAndVolumeOverlap( const Sphere &SPHERE) const
    {
      // calculat the distance between the center of the spheres, init overlap to true
      storage::Pair< double, bool> distance_overlap( linal::Distance( m_Position, SPHERE.m_Position), true);

      // if distance is larger than sum of radii, then there is no overlap
      if( distance_overlap.First() >= m_Radius + SPHERE.m_Radius)
      {
        distance_overlap.Second() = false;
        return distance_overlap;
      }

      // else there is volume overlap
      return distance_overlap;
    }

    //! @brief determine the surface overlap of two spheres analytically
    //! http://mathworld.wolfram.com/Sphere-SphereIntersection.html
    //! @param SPHERE the sphere that overlaps with this sphere
    //! @return the surface overlap of two spheres - 0 is they do not overlap, or one sphere is smaller than the other
    //!         and the surfaces do not touch
    storage::VectorND< 2, double> Sphere::SurfaceOverlap( const Sphere &SPHERE) const
    {
      // calculat the distance between the center of the spheres
      const storage::Pair< double, bool> distance_overlap( DistanceAndVolumeOverlap( SPHERE));

      // if there is no volume overlap, there is no surface overlap
      if( !distance_overlap.Second())
      {
        return 0.0;
      }

      // if distance + radius of either sphere is smaller equal than radius of the other spher, the smaller sphere is
      // within the other sphere, and there is no overlap
      if( distance_overlap.First() + m_Radius <= SPHERE.m_Radius || distance_overlap.First() + SPHERE.m_Radius <= m_Radius)
      {
        return 0.0;
      }

      // calculate height of both caps
      const storage::VectorND< 2, double> cap_heights( HeightOfCaps( distance_overlap.First(), SPHERE));

      // calculate surface of two caps and return
      return storage::VectorND< 2, double>
      (
        2 * math::g_Pi * m_Radius * cap_heights.First(),
        2 * math::g_Pi * SPHERE.m_Radius * cap_heights.Second()
      );
    }

    //! @brief Volume overlap of two spheres
    //! http://mathworld.wolfram.com/Sphere-SphereIntersection.html
    //! @param SPHERE the sphere that overlaps with this sphere
    //! @return the volume of the overlapping lens of this sphere with argument SPHERE - 0 if no overlap
    double Sphere::VolumeOverlap( const Sphere &SPHERE) const
    {
      // calculat the distance between the center of the spheres
      const storage::Pair< double, bool> distance_overlap( DistanceAndVolumeOverlap( SPHERE));

      // if there is no volume overlap, there is no surface overlap
      if( !distance_overlap.Second())
      {
        return 0.0;
      }

      // if one radius is smaller than the other, and the sphere is completely internal
      if( distance_overlap.First() + m_Radius <= SPHERE.m_Radius)
      {
        return GetVolume();
      }
      if( distance_overlap.First() + SPHERE.m_Radius <= m_Radius)
      {
        return SPHERE.GetVolume();
      }

      // actual partial overlap
      // calculate height of both caps
      const storage::VectorND< 2, double> cap_heights( HeightOfCaps( distance_overlap.First(), SPHERE));

      // sum of both volumes
      return math::g_Pi / ( 12 * distance_overlap.First()) *
             math::Sqr( m_Radius + SPHERE.m_Radius - distance_overlap.First()) *
             (
               math::Sqr( distance_overlap.First())
               + 2 * distance_overlap.First() * SPHERE.m_Radius
               - 3 * math::Sqrt( SPHERE.m_Radius)
               + 2 * distance_overlap.First() * m_Radius
               + 6 * SPHERE.m_Radius * m_Radius
               - 3 * math::Sqr( m_Radius)
            );
    }

    //! @brief represent surface by points using latitude sampling
    // http://www.mpip-mainz.mpg.de/~deserno/science_notes/sphere_equi/sphere_equi.ps
    //! @param NUMBER number of points - the actual number of points can slightly deviate
    //! @return vector of points, representing the surface
    storage::Vector< linal::Vector3D> Sphere::PointsOnSurfaceByLatiudes( const size_t NUMBER) const
    {
      // allocate space for points
      storage::Vector< linal::Vector3D> points;
      points.AllocateMemory( NUMBER);

      // surface occupied per point
      const double surface( GetSurface() / NUMBER);
      // distance
      const double distance( math::Sqrt( surface));

      // number of latitudes
      const size_t n_theta( size_t( math::g_Pi / distance));

      // latitudal distance
      const double d_theta( math::g_Pi / n_theta);

      // distance within latidue
      const double d_phi( surface / d_theta);

      const double pi_over_n_theta( math::g_Pi / n_theta);
      const double two_pi_over_d_phi( 2 * math::g_Pi / d_phi);

      // iterate over latitudes
      for( size_t m( 0); m < n_theta; ++m)
      {
        // angle of current latitude
        const double theta( pi_over_n_theta * ( m + 0.5));
        // number of points on that latitude
        const size_t n_phi( size_t( two_pi_over_d_phi * std::sin( theta)));

        const double two_pi_over_n_phi( 2 * math::g_Pi / n_phi);
        // iterate over all points on that longitude
        for( size_t n( 0); n < n_phi; ++n)
        {
          // position on latitude
          const double phi( two_pi_over_n_phi * n);

          // insert the point
          points.PushBack( SphericalCoordinates( m_Radius, phi, theta).Translate( m_Position));
        }
      }

      // end
      return points;
    }

    //! @brief represent surface by points randomly place on surface
    // http://www.mpip-mainz.mpg.de/~deserno/science_notes/sphere_equi/sphere_equi.ps
    //! @param NUMBER number of points
    //! @return vector of points, representing the surface
    storage::Vector< linal::Vector3D> Sphere::PointsOnSurfaceRandom( const size_t NUMBER) const
    {
      // allocate space for points
      storage::Vector< linal::Vector3D> points;
      points.AllocateMemory( NUMBER);

      const double r_sqr( math::Sqr( m_Radius));

      for( size_t n( 0); n < NUMBER; ++n)
      {
        // random z
        const double z( random::GetGlobalRandom().Random< double>( -m_Radius, m_Radius));

        // some math
        const double var( math::Sqrt( r_sqr - math::Sqr( z)));
        // random angle
        const double phi( random::GetGlobalRandom().Random< double>( 2 * math::g_Pi));
        points.PushBack
        (
          linal::Vector3D
          (
            var * std::cos( phi),
            var * std::sin( phi),
            z
          ).Translate( m_Position)
        );
      }

      // end
      return points;
    }

    //! @brief represents surface by points using a spiral sampling
    //! http://sitemason.vanderbilt.edu/page/hmbADS#spiral
    //! @param NUMBER number of points - the actual number of points can slightly deviate
    //! @return vector of points, representing the surface
    storage::Vector< linal::Vector3D> Sphere::PointsOnSurfaceSpiral( const size_t NUMBER) const
    {
      // allocate space for points
      storage::Vector< linal::Vector3D> points;
      points.AllocateMemory( NUMBER);

      const double p( 0.5);
      const double a( 1.0 - 2.0 * p / ( NUMBER - 3));
      const double b( p * double( NUMBER + 1) / double( NUMBER - 3));

      const double empiric_factor_sqrt_n( 7.2 / math::Sqrt( double( NUMBER)));

      double r_prev( 0);
      double phi_prev( 0);

      // insert first point on pole
      points.PushBack( SphericalCoordinates( m_Radius, phi_prev, math::g_Pi).Translate( m_Position));

      for( size_t k( 2); k < NUMBER; ++k)
      {
        const double k_prime( a * k + b);
        const double h( -1 + 2 * ( k_prime - 1) / ( NUMBER - 1));
        const double r( math::Sqrt( 1.0 - math::Sqr( h)));
        const double theta( std::acos( h));
        double phi( phi_prev + empiric_factor_sqrt_n / ( r_prev + r));

        // reduce phi by 2 * pi
        if( phi > 2 * math::g_Pi)
        {
          phi -= 2 * math::g_Pi;
        }

        // insert new point
        points.PushBack( SphericalCoordinates( m_Radius, phi, theta).Translate( m_Position));

        // update previous
        r_prev = r;
        phi_prev = phi;
      }

      // insert last point on pole
      points.PushBack( SphericalCoordinates( m_Radius, 0, 0).Translate( m_Position));

      // end
      return points;
    }

    //! @brief calculate the free surface area fraction, that is not included by any given overlapping sphere
    //! @param SPHERES a list of spheres to be considered
    //! @param NUMBER the number of points that are used to interpolate the surface of the sphere - the higher, the more accurate
    //! @return fraction [0,1] where 0 is no free surface area, 1 is no overlap with any sphere
    double Sphere::FreeSurfaceAreaFraction( const util::SiPtrList< const Sphere> &SPHERES, const size_t NUMBER) const
    {
      BCL_Assert( NUMBER > 0, "need at least one point to interpolate surface of a sphere!")

      // are there any overlapping spheres
      if( SPHERES.IsEmpty())
      {
        return 1.0;
      }

      // there are overlapping spheres, interpolate surface with points
      const storage::Vector< linal::Vector3D> interpolate_surface_points( PointsOnSurfaceSpiral( NUMBER));
//      const storage::Vector< linal::Vector3D> interpolate_surface_points( PointsOnSurfaceByLatiudes( NUMBER));
//      const storage::Vector< linal::Vector3D> interpolate_surface_points( PointsOnSurfaceRandom( NUMBER));

      // counter for overlapping points
      size_t number_overlapping_points( 0);

      // for each point interpolating the surface, check if there is a sphere, that overlaps that point
      for
      (
        storage::Vector< linal::Vector3D>::const_iterator
          isp_itr( interpolate_surface_points.Begin()), isp_itr_end( interpolate_surface_points.End());
        isp_itr != isp_itr_end;
        ++isp_itr
      )
      {
        // iterate over all overlapping spheres
        for
        (
          util::SiPtrList< const Sphere>::const_iterator os_itr( SPHERES.Begin()), os_itr_end( SPHERES.End());
          os_itr != os_itr_end;
          ++os_itr
        )
        {
          // does this sphere contain the point
          if( ( *os_itr)->DoesContain( *isp_itr))
          {
            // increase the points that are overlapping with any sphere
            ++number_overlapping_points;
            // no need to check if point overlaps with any other sphere
            break;
          }
        }
      }

      // insert the fraction of the sphere sizes
      return double( NUMBER - number_overlapping_points) / double( NUMBER);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! writes to ostream
    std::ostream &Sphere::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Position, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Radius  , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::ostream
    //! @param ISTREAM input stream
    //! @return ostream which was read from
    std::istream &Sphere::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Position, ISTREAM);
      io::Serialize::Read( m_Radius  , ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief derive cap heights of overlapping spheres
    //! @param DISTANCE the distance of the two caps
    //! @param SPHERE the other sphere
    //! @return pair of heights of caps - first for this cap, second for the other
    storage::VectorND< 2, double> Sphere::HeightOfCaps( const double DISTANCE, const Sphere &SPHERE) const
    {
      // calculate the height of the two caps
      return storage::VectorND< 2, double>
      (
        ( SPHERE.m_Radius - m_Radius + DISTANCE) * ( SPHERE.m_Radius + m_Radius - DISTANCE) /
        ( 2 * DISTANCE),
        ( m_Radius - SPHERE.m_Radius + DISTANCE) * ( m_Radius + SPHERE.m_Radius - DISTANCE) /
        ( 2 * DISTANCE)
      );
    }

    //! @brief list overlapping spheres for each sphere
    //! @param SPHERES all spheres that are considered
    //! @return a map, that has a list of overlapping spheres (data) for each sphere (key)
    storage::Map< util::SiPtr< const Sphere>, util::SiPtrList< const Sphere> >
    ListOverlaps( const util::SiPtrVector< const Sphere> &SPHERES)
    {
      storage::Map< util::SiPtr< const Sphere>, util::SiPtrList< const Sphere> > overlaps;

      // need at least two spheres
      if( SPHERES.GetSize() <= 1)
      {
        return overlaps;
      }

      // iterate over all spheres
      for( util::SiPtrVector< const Sphere>::const_iterator itr1( SPHERES.Begin()), itr_end( SPHERES.End()); itr1 != itr_end; ++itr1)
      {
        // iterate over all other spheres
        for( util::SiPtrVector< const Sphere>::const_iterator itr2( itr1 + 1); itr2 != itr_end; ++itr2)
        {
          if( ( *itr1)->DoesOverlap( **itr2))
          {
            overlaps[ *itr1].PushFront( *itr2);
            overlaps[ *itr2].PushFront( *itr1);
          }
        }
      }

      // end
      return overlaps;
    }

    //! @brief calculate the solvent accessible surface area for each given sphere
    //! this is achieved by interpolating the area of a sphere by equally distributing NUMBER of points on the surface
    //! of the sphere and calculating the ration of points that are within the radius of any other sphere to NUMBER
    //! @see Shrake A, Rupley JA. (1973). Environment and exposure to solvent of protein atoms. Lysozyme and insulin. J Mol Biol 79(2):351-71.
    //! @param SPHERES the spheres to be considered
    //! @param NUMBER number of points to interpolate the surface of the sphere - the larger the more accurate, but slower
    //! @return a map that has an absolute SASA for each sphere
    storage::Map< util::SiPtr< const Sphere>, double>
    FreeSurface
    (
      const util::SiPtrVector< const Sphere> &SPHERES,
      const size_t NUMBER
    )
    {
      // initialize the relative free surface areas
      storage::Map< util::SiPtr< const Sphere>, double> free_surface_area( FreeSurfaceRatio( SPHERES, NUMBER));

      // iteate over all spheres and their relative surface area
      for
      (
        storage::Map< util::SiPtr< const Sphere>, double>::iterator
          itr( free_surface_area.Begin()), itr_end( free_surface_area.End());
        itr != itr_end;
        ++itr
      )
      {
        // multiply with absolute surface
        itr->second *= itr->first->GetSurface();
      }

      // end
      return free_surface_area;
    }

    //! @brief calculate the solvent accessible surface area ratio for each given sphere
    //! this is achieved by interpolating the area of a sphere by equally distributing NUMBER of points on the surface
    //! of the sphere and calculating the ration of points that are within the radius of any other sphere to NUMBER
    //! @param SPHERES the spheres to be considered
    //! @param NUMBER number of points to interpolate the surface of the sphere - the larger the more accurate, but slower
    //! @return a map that has a relative area for each sphere
    storage::Map< util::SiPtr< const Sphere>, double>
    FreeSurfaceRatio
    (
      const util::SiPtrVector< const Sphere> &SPHERES,
      const size_t NUMBER
    )
    {
      BCL_Assert( NUMBER > 0, "need at least one point to interpolate surface of a sphere!")

      // Initialize the free surface areas
      storage::Map< util::SiPtr< const Sphere>, double> free_surface_area;

      // the search of the overlaps (should be done outside for a function RecalculateFreeSurface() ).
      const storage::Map< util::SiPtr< const Sphere>, util::SiPtrList< const Sphere> > overlaps( ListOverlaps( SPHERES));

      // iterate over all spheres
      for
      (
        storage::Map< util::SiPtr< const Sphere>, util::SiPtrList< const Sphere> >::const_iterator
          sphere_itr( overlaps.Begin()), sphere_itr_end( overlaps.End());
        sphere_itr != sphere_itr_end;
        ++sphere_itr
      )
      {
        // current sphere
        const util::SiPtr< const Sphere> &current_sphere( sphere_itr->first);

        // insert the fraction of the sphere sizes
        free_surface_area[ current_sphere] = current_sphere->FreeSurfaceAreaFraction( sphere_itr->second, NUMBER);
      }

      // end
      return free_surface_area;
    }

  } // namespace coord
} // namespace bcl
