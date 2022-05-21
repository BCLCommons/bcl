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
#include "coord/bcl_coord.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_polygon.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_2d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_running_min_max.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! @brief checks that all coordinates given in the vector COORDINATES are defined
    //! @param COORDINATES vector of coordinates of interest
    //! @return whether all coordinates given in the vector COORDINATES are defined
    bool AreDefinedCoordinates
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    )
    {
      //iterate over all given coordinates to check each whether they are defined
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        //if current coordinate is not defined return false
        if( !( *coord_itr)->IsDefined())
        {
          return false;
        }
      }
      // otherwise return true
      return true;
    }

    //! @brief calculates the cartesian coordinates from the distance (relative to Zero point) and two angles
    //! @param DISTANCE distance to origin
    //! @param PHI first angle
    //! @param PSI second angle
    //! @return the cartesian coordinates for given distance and two angles
    linal::Vector3D SphericalCoordinates( const double DISTANCE, const double PHI, const double PSI)
    {
      // construct the cartesian coordinates and return it
      return linal::Vector3D
             (
               DISTANCE * cos( PHI) * sin( PSI),
               DISTANCE * sin( PHI) * sin( PSI),
               DISTANCE *             cos( PSI)
             );
    }

    //! @brief transforms the given vector of coordinates by given TransformationMatrix3D
    //! @param COORDINATES Vector of coordinates of interest
    //! @param TRANSFORMATION_MATRIX TransformationMatrix3D to be used
    void TransformCoordinates
    (
      util::SiPtrVector< linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATION_MATRIX
    )
    {
      // iterate over coordinates
      for
      (
        util::SiPtrVector< linal::Vector3D>::iterator ptr( COORDINATES.Begin()), ptr_end( COORDINATES.End());
        ptr != ptr_end;
        ++ptr
      )
      {
        // apply the transformation
        ( *ptr)->Transform( TRANSFORMATION_MATRIX);
      }
    }

    //! @brief calculate a distance matrix from given vector of coordinates
    //! @param COORDINATES Vector of coordinates of interest
    //! @return distance matrix for given coordinates
    linal::Matrix< double> CalculateDistanceMatrix
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    )
    {
      const size_t nr_coordiantes( COORDINATES.GetSize());

      // initialize distance matrix to be returned
      linal::Matrix< double> distance_matrix( nr_coordiantes, nr_coordiantes);

      // iterate over the coordinates
      for( size_t i( 0); i < nr_coordiantes; ++i)
      {
        // iterate over the coordinates that come after
        for( size_t j( i + 1); j < nr_coordiantes; ++j)
        {
          // calculate the distance
          const double distance( ( *COORDINATES( i) - *COORDINATES( j)).Norm());

          // copy distance to corresponding cells in the matrix
          distance_matrix( i, j) = distance;
          distance_matrix( j, i) = distance;
        }
      }
      // end
      return distance_matrix;
    }

    //! @brief calculate difference matrix between two coordinate vectors
    //! @param COORDINATES_A first vector of coordinates of interest
    //! @param COORDINATES_B second vector of coordinates of interest
    //! @return difference matrix between given coordinate vectors
    linal::Matrix< double> CalculateDifferenceDistanceMatrix
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_A,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES_B
    )
    {
      // return the difference between the distance matrices of given coordinates vector pair
      return CalculateDistanceMatrix( COORDINATES_A) - CalculateDistanceMatrix( COORDINATES_B);
    }

    //! computes center of mass of COORDINATES
    linal::Vector3D CenterOfMass( const util::SiPtrVector< const linal::Vector3D> &COORDINATES, const bool SKIP_UNDEFINED)
    {
      //instantiate center
      linal::Vector3D center( double( 0.0));
      size_t count( 0);
      //sum all coordinates up
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()),
          coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        // skip undefined coordinates
        if( SKIP_UNDEFINED && !( *coord_itr)->IsDefined())
        {
          continue;
        }

        center += ( **coord_itr);
        ++count;
      }

      // divide center by the number of coordinates
      if( count > 0)
      {
        center /= double( count);
      }

      return center;
    }

    // compute the square of the radius of gyration for a set of coordinates
    double SquareRadiusOfGyration( const util::SiPtrVector< const linal::Vector3D> &COORDINATES)
    {
      // calculate center of mass
      linal::Vector3D center_of_mass( 0.0, 0.0, 0.0);
      double sum_square_norm( 0.0);

      size_t count( 0);
      // add up square of Euclidean distance to the center of mass of each point
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        // skip undefined coordinates
        if( !( *coord_itr)->IsDefined())
        {
          continue;
        }

        center_of_mass += ( **coord_itr);
        sum_square_norm += ( *coord_itr)->SquareNorm();
        ++count;
      }

      // center of mass as average coordinate
      center_of_mass /= double( count);
      sum_square_norm /= double( count);

      // actual square radius of gyration
      const double square_radius_of_gyration( sum_square_norm - center_of_mass.SquareNorm());

      //end
      return square_radius_of_gyration;
    }

    //! compute the radius of gyration for a set of coordinates
    double RadiusOfGyration( const util::SiPtrVector< const linal::Vector3D> &COORDINATES)
    {
      return math::Sqrt( SquareRadiusOfGyration( COORDINATES));
    }

    //! computes the symmetrie of a set of coordinates
    double SymmetryFactor( const util::ShPtrVector< util::SiPtrVector< const linal::Vector3D> > &SET_OF_COORDINATES)
    {
      //symmetry factor
      double symmetry( 0);

      // number of symmetryelements
      size_t number_elements( SET_OF_COORDINATES.GetSize());
      // number of coordinates per element
      size_t number_coordinates_per_element( SET_OF_COORDINATES( 0)->GetSize());

      //sum over all distances to the center and between the symmetry elements
      double sum_distances( 0);
      //sum of all scalarproducts between center and symmetryelements
      double sum_scalars( 0);

      //differences iterators on all given coordinatevectors
      storage::Vector< std::vector< util::SiPtr< const linal::Vector3D> >::difference_type>
        coordinates_itr_diff( SET_OF_COORDINATES.GetSize());

      //iterator to the first element in the first coordinatevector
      util::SiPtrVector< const linal::Vector3D>::const_iterator
        coordinates_itr( SET_OF_COORDINATES( 0)->Begin()),
        coordinates_itr_end( SET_OF_COORDINATES( 0)->End());
      coordinates_itr_diff( 0) = 0;
      {
        size_t i( 1);
        //assures that all coordinatevectors are similar (same length)
        for
        (
          util::ShPtrVector< util::SiPtrVector< const linal::Vector3D> >::const_iterator
            set_itr( SET_OF_COORDINATES.Begin()), set_itr_end( SET_OF_COORDINATES.End());
          set_itr + 1 != set_itr_end; ++set_itr)
        {
          BCL_Assert
          (
            ( *set_itr)->GetSize() == ( *( set_itr + 1))->GetSize(),
            "symmetry of sets with different number of coordinates"
          );
          //store difference of to iterators of all other sequences to first sequence iterators
          coordinates_itr_diff( i++)
            = util::SiPtrVector< const linal::Vector3D>::const_iterator( ( *( set_itr + 1))->Begin()) - coordinates_itr;
        }
      }

      //calculate the the distance to the center and between all coordinate vectors
      storage::Vector< linal::Vector3D> centers( SET_OF_COORDINATES( 0)->GetSize());
      for
      (
        std::vector< linal::Vector3D>::iterator center_itr( centers.Begin()), center_itr_end( centers.End());
        coordinates_itr != coordinates_itr_end && center_itr != center_itr_end;
        ++center_itr, ++coordinates_itr
      )
      {
        for
        (
          std::vector< std::vector< util::SiPtr< const linal::Vector3D> >::difference_type>::const_iterator
            diff_itr( coordinates_itr_diff.Begin()),
            diff_itr_end( coordinates_itr_diff.End());
          diff_itr != diff_itr_end; ++diff_itr
        )
        {
          ( *center_itr) += ( **( coordinates_itr + ( *diff_itr)));
        }
        ( *center_itr) /= double( number_elements);

        //calculate sum of differences between distances of coordinates in coordinate vector and center
        for
        (
          std::vector< std::vector< util::SiPtr< const linal::Vector3D> >::difference_type>::const_iterator
            diff_itr_a( coordinates_itr_diff.Begin()),
            diff_itr_end( coordinates_itr_diff.End());
          diff_itr_a + 1 != diff_itr_end; ++diff_itr_a
        )
        {
          for
          (
            std::vector< std::vector< util::SiPtr< const linal::Vector3D> >::difference_type>::const_iterator
              diff_itr_b( diff_itr_a + 1);
            diff_itr_b != diff_itr_end; ++diff_itr_b)
          {
            sum_distances += math::Sqr
                             (
                               ( ( **( coordinates_itr + ( *diff_itr_a))) - ( *center_itr)).Norm() -
                               ( ( **( coordinates_itr + ( *diff_itr_b))) - ( *center_itr)).Norm()
                             );
          }
        }
      }

      //reset coordinates_itr to first position
      coordinates_itr = SET_OF_COORDINATES( 0)->Begin();

      //sum up the scalars
      for
      (
        std::vector< linal::Vector3D>::const_iterator
          center_itr_a( centers.Begin()),
          center_itr_end( centers.End());
        center_itr_a != center_itr_end; ++center_itr_a, ++coordinates_itr
      )
        for
        (
          std::vector< linal::Vector3D>::const_iterator
            center_itr_b( centers.Begin());
          center_itr_b != center_itr_end; ++center_itr_b
        )
        {
          if( center_itr_a == center_itr_b) continue;
          for
          (
            std::vector< std::vector< util::SiPtr< const linal::Vector3D> >::difference_type>::const_iterator
              diff_itr( coordinates_itr_diff.Begin()), diff_itr_end( coordinates_itr_diff.End());
            diff_itr != diff_itr_end; ++diff_itr
          )
          {
            linal::Vector3D cr( ( ( **( coordinates_itr + ( *diff_itr))) - ( *center_itr_a))),
                           cc( ( *center_itr_b) - ( *center_itr_a));

            sum_scalars += math::Absolute( cr * cc);
//            sum_scalars += math::Absolute( linal::ProjAngleCosinus( ( *center_itr1), ( **( coordinates_itr + ( *diff_itr))), ( *center_itr1), ( *center_itr2)));

          }
        }

      //add distances and scalars
      symmetry =   sum_distances * 2 / ( number_elements * ( number_elements - 1) * number_coordinates_per_element)
                 + sum_scalars / ( math::Sqr( number_coordinates_per_element) * number_elements);

      //squareroot
      symmetry = math::Sqrt( symmetry);

      //end
      return symmetry;
    } // end SymmetryFactor

    //! NeighborWeight - assign each amino amid a weight (0 <= weight <= 1) based on its distance from the amino acid of interest
    //! if distance <= low_threshold --> weight = 1
    //! if low_threshold < distance < high_threshold --> weight = piece of cos function
    //! if distance >= high_threshold --> weight = 0
    //! @param DISTANCE - the distance between the current AA and the AA of interest
    //! @param LOW_THRESHOLD - the lower bound.  all neighbors within this distance are given a weight of 1
    //! @param HIGH_THRESHOLD - the upper bound.  all neighbors further than this distance are given a weight of 0
    //! @return weight - the weight assigned to the current AA.  0 <= weight <= 1
    double NeighborWeight( const double DISTANCE, const double LOW_THRESHOLD, const double HIGH_THRESHOLD)
    {
      //difference between thresholds
      const double difference( HIGH_THRESHOLD - LOW_THRESHOLD);

      if( DISTANCE <= LOW_THRESHOLD)
      {
        return double( 1.0);
      }
      else if( DISTANCE >= HIGH_THRESHOLD)
      {
        return double( 0.0);
      }
      else
      {
        return ( std::cos( ( DISTANCE - LOW_THRESHOLD) / ( difference) * math::g_Pi) + double( 1.0)) / double( 2.0);
      }

    } // end NeighborWeight

    //! @brief CountNeighbors - calculates the number of amino acid neighbors that the current amino acid has
    //! @param ALL_POINTS - all neighbor candidates.  we examine them to determine
    //!        which ones are actually close enough to be counted as neighbors
    //! @param CURRENT_POINT - the amino acid whose neighbors we are counting
    //! @param LOW_THRESHOLD - the lower bound for neighbor count - all amino acids within this distance are
    //!        counted as a neighbor
    //! @param HIGH_THRESHOLD - the upper bound for neighbor count - all neighbors further than this distance
    //!        are given a weight of 0
    //! all amino acids between low_threshold and high_threshold are given a weight between 0 and 1
    //! @return neighbor_count - how many neighbors current_amino acid has
    double
    CountNeighbors
    (
      const util::SiPtrVector< const linal::Vector3D> &ALL_POINTS,
      const linal::Vector3D &CURRENT_POINT,
      const double LOW_THRESHOLD,
      const double HIGH_THRESHOLD
    )
    {
      // start with the assumption that the current AA has no neighbors
      double neighbor_count( 0.0);

      // loop over all AAs to figure out which ones are neighbors
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          point_itr( ALL_POINTS.Begin()), point_itr_end( ALL_POINTS.End());
        point_itr != point_itr_end; ++point_itr
      )
      {
        if( !( *point_itr)->IsDefined()) continue;

        // calculate the distance between the current amino acid and the current neighbor being considered
        const double dist( ( **point_itr - CURRENT_POINT).Norm());

        neighbor_count += NeighborWeight( dist, LOW_THRESHOLD, HIGH_THRESHOLD);
      }
      return neighbor_count;
    } // end CountNeighbors

    //! @brief NeighborVector - calculates the sum of all vectors between the current AA and its neighboring AAs
    //! @param ALL_POINTS- all neighbor candidates.
    //! @param CURRENT_POINT - the amino acid whose neighbors we are counting
    //! @param LOW_THRESHOLD - the lower bound for neighbor count - all amino acids within this distance are
    //!        counted as a neighbor
    //! @param HIGH_THRESHOLD - the upper bound for neighbor count - all neighbors further than this distance
    //!        are given a weight of 0
    //! @param NORMALIZE_BY_NEIGHBOR_COUNT = true - whether or not the sum of neighbor vectors is normalized
    //!        by the neighbor count
    linal::Vector3D NeighborVector
    (
      const util::SiPtrVector< const linal::Vector3D> &ALL_POINTS,
      const linal::Vector3D &CURRENT_POINT,
      const double LOW_THRESHOLD,
      const double HIGH_THRESHOLD,
      const bool NORMALIZE_BY_NEIGHBOR_COUNT
    )
    {
      // start with the assumption that the current AA has no neighbors
      double neighbor_count( 0.0);
      // vector variable used to hold the sums of all vectors between the current AA and neighboring AAs
      linal::Vector3D vector_sum;

      // if this AA is not defined, get out of this function and return error vector
      if( !CURRENT_POINT.IsDefined())
      {
        return linal::Vector3D( util::GetUndefined< double>());
      }

      //loop over all AAs
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          point_itr( ALL_POINTS.Begin()), point_itr_end( ALL_POINTS.End());
        point_itr != point_itr_end; ++point_itr
      )
      {

        // if this AA is not defined, move along to the next AA
        if( !( *point_itr)->IsDefined()) continue;

        //calculate distance between the current amino acid and the current neighbor being considered
        const double dist( ( ( **point_itr) - CURRENT_POINT).Norm());

        // determine the weight assigned to this neighbor
        const double neighbor_weight = NeighborWeight( dist, LOW_THRESHOLD, HIGH_THRESHOLD);

        // check to ensure it's within the threshhold
        if( neighbor_weight != double( 0.0))
        {
          neighbor_count += neighbor_weight;
          // normalize each neighbor vector such that each gets an initial weight of 1 (before the neighbor_weight)
          vector_sum += ( ( **point_itr) - CURRENT_POINT).Normalize() * neighbor_weight;
        }
      }

      if( NORMALIZE_BY_NEIGHBOR_COUNT)
      {
        vector_sum /= neighbor_count;
      }

      return vector_sum;
    } // end NeighborVector

    //! @brief calculate the density of a given set of points
    //! @param COORDINATES the coordinates to be used
    //! @param RESOLUTION the resolution the density is determined at in x, y and z
    //! @return the density in length_unit^3
    double CalculatePointDensity
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const linal::Vector3D &RESOLUTION
    )
    {
      // store the cubes that contain the points
      const storage::Set< storage::Triplet< int, int, int> > cubes
      (
        QuantizePoints( COORDINATES, RESOLUTION)
      );

      BCL_MessageDbg
      (
        "number coordinates: " + util::Format()( COORDINATES.GetSize()) +
        " number cubes occupied: " + util::Format()( cubes.GetSize())
      );

      // divide the number of coordinates by the number of cubes occupied and their volume
      return double( COORDINATES.GetSize()) /
             ( double( cubes.GetSize()) * RESOLUTION.X() * RESOLUTION.Y() * RESOLUTION.Z());
    }

    //! @brief quantize a given set of points to a certain resolution
    //! @param COORDINATES the coordinates to be used
    //! @param RESOLUTION the resolution the points are quantized to in x,y and z
    //! @return a set of 3d indices to a grid of the size RESOLUTION
    storage::Set< storage::Triplet< int, int, int> > QuantizePoints
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const linal::Vector3D &RESOLUTION
    )
    {
      // store the cubes that contain the points
      storage::Set< storage::Triplet< int, int, int> > cubes;

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        // insert the point as quantized cube
        cubes.Insert
        (
          storage::Triplet< int, int, int>
          (
            int( std::floor( ( *coord_itr)->X() / RESOLUTION.X())),
            int( std::floor( ( *coord_itr)->Y() / RESOLUTION.Y())),
            int( std::floor( ( *coord_itr)->Z() / RESOLUTION.Z()))
          )
        );
      }

      // end
      return cubes;
    }

    //! @brief Estimate volume of points using 2D convex hulls
    //! @param COORDINATES the coordinates to be used
    //! @param SLICE_WIDTH typical distance between adjacent points
    //! @param MAX_DISTANCE maximum distance between points; if above this, induce a concavity into the convex hull
    double EstimateVolume
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const double &SLICE_WIDTH,
      const double &MAX_DISTANCE,
      const double &RADIUS
    )
    {
      const double rad_cubed( RADIUS * RADIUS * RADIUS);
      if( COORDINATES.IsEmpty())
      {
        return 0.0;
      }
      else if( COORDINATES.GetSize() == size_t( 1))
      {
        return 4.0 / 3.0 * math::g_Pi * rad_cubed;
      }
      // most dense possible sphere packing leaves 74% of space empty
      const double min_volume
      (
        4.0 / 3.0 * math::g_Pi * rad_cubed * double( COORDINATES.GetSize()) / 0.74
      );
      math::RunningMinMax< linal::Vector3D> minmax;
      for( auto itr( COORDINATES.Begin()), itr_end( COORDINATES.End()); itr != itr_end; ++itr)
      {
        minmax += **itr;
      }
      const linal::Vector3D &mins( minmax.GetMin());
      const linal::Vector3D &maxs( minmax.GetMax());

      if( !mins.IsDefined() || !maxs.IsDefined() || mins == maxs)
      {
        return min_volume;
      }

      double est_volume( 0.0);
      for( size_t axis( 0), mx_axis( 3); axis < mx_axis; ++axis)
      {
        const size_t n_divisions( ( maxs( axis) - mins( axis)) / SLICE_WIDTH + 1);
        storage::Vector< storage::Vector< linal::Vector2D> > slice_coords( n_divisions);
        const double ma( mins( axis));
        const size_t dim1( axis ? 0 : 1);
        const size_t dim2( axis == 2 ? 1 : 2);

        for( auto itr( COORDINATES.Begin()), itr_end( COORDINATES.End()); itr != itr_end; ++itr)
        {
          const size_t slice( ( ( **itr)( axis) - ma) / SLICE_WIDTH);
          slice_coords( slice).PushBack( linal::Vector2D( ( **itr)( dim1), ( **itr)( dim2)));
        }
        linal::Vector< double> polygon_areas( n_divisions);
        int i( 0);
        for( auto itr( slice_coords.Begin()), itr_end( slice_coords.End()); itr != itr_end; ++itr, ++i)
        {
          polygon_areas( i) = Polygon::ConvexHull( *itr, MAX_DISTANCE, RADIUS).Expand( RADIUS).GetArea();
        }
        double vol_so_far( 0.0);
        if( n_divisions == 1)
        {
          vol_so_far += polygon_areas( 0);
        }
        else if( n_divisions)
        {
          for( size_t i( 1); i < n_divisions; ++i)
          {
            vol_so_far += ( polygon_areas( i) + polygon_areas( i - 1) + math::Sqrt( polygon_areas( i) * polygon_areas( i - 1))) / 3.0;
          }
          vol_so_far += ( polygon_areas( 0) + polygon_areas( polygon_areas.GetSize() - 1)) / 2.0;
        }
        est_volume += vol_so_far * SLICE_WIDTH / 3.0;
      }
      return std::max( double( est_volume), min_volume);
    }

  } // namespace coord
} // namespace bcl

