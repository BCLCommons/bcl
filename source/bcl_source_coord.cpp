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
#include "coord/bcl_coord_axes.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Axes::Axes() :
      e_X( AddEnum( "X_Axis", linal::Vector3D( 1.0, 0.0, 0.0))),
      e_Y( AddEnum( "Y_Axis", linal::Vector3D( 0.0, 1.0, 0.0))),
      e_Z( AddEnum( "Z_Axis", linal::Vector3D( 0.0, 0.0, 1.0)))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Axes::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Axes::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Axes::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief construct on access function for all Axes
    //! @return reference to only instances of Axes
    const Axes &GetAxes()
    {
      return Axes::GetEnums();
    }

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< linal::Vector3D, coord::Axes>;

  } // namespace util
} // namespace bcl
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
#include "coord/bcl_coord_cyclic_coordinate_descent.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_line_segment_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CyclicCoordinateDescent::TargetAndMovingPointPair::s_Instance
    (
      GetObjectInstances().AddInstance( new CyclicCoordinateDescent::TargetAndMovingPointPair())
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CyclicCoordinateDescent::ResultsAndCoefficients::s_Instance
    (
      GetObjectInstances().AddInstance( new CyclicCoordinateDescent::ResultsAndCoefficients())
    );

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CyclicCoordinateDescent::s_Instance
    (
      GetObjectInstances().AddInstance( new CyclicCoordinateDescent())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    CyclicCoordinateDescent::TargetAndMovingPointPair::TargetAndMovingPointPair() :
      m_TargetPoint(),
      m_MovingPoint()
    {
    }

    //! constructor taking a target point and a moving point
    //! @param TARGET_POINT the fixed target point desired to be found
    //! @param MOVING_POINT the point which is being moved in an attempt to reach the target point
    CyclicCoordinateDescent::TargetAndMovingPointPair::TargetAndMovingPointPair
    (
      const linal::Vector3D &TARGET_POINT, const linal::Vector3D &MOVING_POINT
    ) :
      m_TargetPoint( TARGET_POINT),
      m_MovingPoint( MOVING_POINT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new TargetAndMovingPointPair
    CyclicCoordinateDescent::TargetAndMovingPointPair *CyclicCoordinateDescent::TargetAndMovingPointPair::Clone() const
    {
      return new TargetAndMovingPointPair( *this);
    }

    //! @brief default constructor
    CyclicCoordinateDescent::ResultsAndCoefficients::ResultsAndCoefficients() :
      m_OptimalRotation(),
      m_MinimumDistanceSum(),
      m_CoefficientA(),
      m_CoefficientB(),
      m_CoefficientC()
    {
    }

    //! @brief constructor taking all of the member variables
    //! @param OPTIMAL_ROTATION the rotation to minimize the distance sum between target and moving points
    //! @param MINIMUM_DISTANCE_SUM the sum of distances between target and moving points that will result after
    //!        the moving points are rotated by "m_OptimalRotation"
    //! @param COEFFICIENT_A please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_B please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_C please see the reference given above (Canutescu et al. Protein Sci, 2003)
    CyclicCoordinateDescent::ResultsAndCoefficients::ResultsAndCoefficients
    (
      const double OPTIMAL_ROTATION, const double MINIMUM_DISTANCE_SUM, const double COEFFICIENT_A,
      const double COEFFICIENT_B, const double COEFFICIENT_C
    ) :
      m_OptimalRotation( OPTIMAL_ROTATION),
      m_MinimumDistanceSum( MINIMUM_DISTANCE_SUM),
      m_CoefficientA( COEFFICIENT_A),
      m_CoefficientB( COEFFICIENT_B),
      m_CoefficientC( COEFFICIENT_C)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ResultsAndCoefficients
    CyclicCoordinateDescent::ResultsAndCoefficients *CyclicCoordinateDescent::ResultsAndCoefficients::Clone() const
    {
      return new ResultsAndCoefficients( *this);
    }

    //! @brief default constructor
    CyclicCoordinateDescent::CyclicCoordinateDescent()
    {
    }

    //! @brief Clone function
    //! @return pointer to new CyclicCoordinateDescent
    CyclicCoordinateDescent *CyclicCoordinateDescent::Clone() const
    {
      return new CyclicCoordinateDescent( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CyclicCoordinateDescent::TargetAndMovingPointPair::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CyclicCoordinateDescent::ResultsAndCoefficients::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CyclicCoordinateDescent::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief GetTargetPoint gives the fixed target point desired to be found
    //! @return returns linal::Vector3D which is "m_TargetPoint"
    const linal::Vector3D &CyclicCoordinateDescent::TargetAndMovingPointPair::GetTargetPoint() const
    {
      return m_TargetPoint;
    }

    //! @brief GetMovingPoint gives the point which is being moved in an attempt to reach the target point
    //! @return returns linal::Vector3D which is "m_MovingPoint"
    const linal::Vector3D &CyclicCoordinateDescent::TargetAndMovingPointPair::GetMovingPoint() const
    {
      return m_MovingPoint;
    }

    //! @brief GetOptimalRotation gives the rotation to minimize the distance sum between target and moving points
    //! @return returns "m_OptimalRotation"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetOptimalRotation() const
    {
      return m_OptimalRotation;
    }

    //! @brief GetMinimumDistanceSum gives the sum of distances between target and moving points that will result
    //!        after the moving points are rotated by "m_OptimalRotation"
    //! @return returns "m_MinimumDistanceSum"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetMinimumDistanceSum() const
    {
      return m_MinimumDistanceSum;
    }

    //! @brief GetCoefficientA gives "m_CoefficientA"
    //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @return returns "m_CoefficientA"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetCoefficientA() const
    {
      return m_CoefficientA;
    }

    //! @brief GetCoefficientB gives "GetCoefficientB"
    //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @return returns "GetCoefficientB"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetCoefficientB() const
    {
      return m_CoefficientB;
    }

    //! @brief GetCoefficientC gives "GetCoefficientC"
    //!        please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @return returns "GetCoefficientC"
    double CyclicCoordinateDescent::ResultsAndCoefficients::GetCoefficientC() const
    {
      return m_CoefficientC;
    }

    //! @brief GetOptimalRotationandDistance gives rotation to minimize the distance between target and moving points
    //! It calculates the rotation around "ROTATION_AXIS" that needs to be performed on all moving points so that
    //! after the rotation, the sum distance from each of the moving points to corresponding fixed points will be
    //! minimized
    //! @param ROTATION_AXIS this is the axis the rotation will be performed around
    //! @param TARGET_MOVING_POINTS the list of moving points and the corresponding fixed target points the moving
    //!        points are trying to reach
    //! @return returns a ResultsAndCoefficients which holds the results of the calculation
    CyclicCoordinateDescent::ResultsAndCoefficients CyclicCoordinateDescent::GetOptimalRotationandDistance
    (
      const LineSegment3D &ROTATION_AXIS,
      const storage::List< CyclicCoordinateDescent::TargetAndMovingPointPair> &TARGET_MOVING_POINTS
    ) const
    {
      // create variables to hold the coefficients a, b, and c. These will be summed up over the "TARGET_MOVING_POINTS"
      double coefficient_a_sum( 0);
      double coefficient_b_sum( 0);
      double coefficient_c_sum( 0);

      // iterate through "TARGET_MOVING_POINTS" in order to sum up "coefficient_a_sum", "coefficient_b_sum", and
      // "coefficient_c_sum"
      for
      (
        storage::List< CyclicCoordinateDescent::TargetAndMovingPointPair>::const_iterator
          itr( TARGET_MOVING_POINTS.Begin()), itr_end( TARGET_MOVING_POINTS.End());
        itr != itr_end;
        ++itr
      )
      {
        // calculate the coefficient values for the current target and mvoing point pair
        const storage::VectorND< 3, double> current_a_b_c_coefficient_values
        (
          CalculateCurrentCoefficientValues( ROTATION_AXIS, *itr)
        );

        // add the coefficients of "current_a_b_c_coefficient_values" into "coefficient_a_sum", "coefficient_b_sum",
        // and "coefficient_c_sum"
        coefficient_a_sum += current_a_b_c_coefficient_values.First();
        coefficient_b_sum += current_a_b_c_coefficient_values.Second();
        coefficient_c_sum += current_a_b_c_coefficient_values.Third();
      }

      // message the coefficient sums
      BCL_MessageDbg( "coefficient_a_sum " + util::Format()( coefficient_a_sum));
      BCL_MessageDbg( "coefficient_b_sum " + util::Format()( coefficient_b_sum));
      BCL_MessageDbg( "coefficient_c_sum " + util::Format()( coefficient_c_sum));

      // calculate the cosine of alpha
      const double cos_alpha
      (
        coefficient_b_sum /
          math::Sqrt( math::Sqr( coefficient_b_sum) + math::Sqr( coefficient_c_sum))
      );

      // message the cosine of alpha
      BCL_MessageDbg( "cos_alpha " + util::Format()( cos_alpha));

      // calculate the sin of alpha
      const double sin_alpha
      (
        coefficient_c_sum /
          math::Sqrt( math::Sqr( coefficient_b_sum) + math::Sqr( coefficient_c_sum))
      );

      // message the sin of alpha
      BCL_MessageDbg( "sin_alpha " + util::Format()( sin_alpha));

      // calculate theta, this is the optimal rotation which will minimize the sum square distance between the moving
      // and target points
      const double theta( atan2( sin_alpha, cos_alpha));

      // message theta
      BCL_MessageDbg( "theta " + util::Format()( theta));

      // calculate the distance sum that will result after rotating around "ROTATION_AXIS" by "theta"
      const double distance_sum( CalculateDistanceSum( theta, coefficient_a_sum, coefficient_b_sum, coefficient_c_sum));

      // message the resulting distance sum
      BCL_MessageDbg( "distance_sum_no_rot " + util::Format()( CalculateDistanceSum( 0.0, coefficient_a_sum, coefficient_b_sum, coefficient_c_sum)));
      BCL_MessageDbg( "distance_sum " + util::Format()( distance_sum));

      // create a ResultsAndCoefficients
      const CyclicCoordinateDescent::ResultsAndCoefficients optimal_rotation_and_distance_and_coefficients
      (
        theta, distance_sum, coefficient_a_sum, coefficient_b_sum, coefficient_c_sum
      );

      // return the results and coefficients
      return optimal_rotation_and_distance_and_coefficients;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CyclicCoordinateDescent::TargetAndMovingPointPair::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TargetPoint, ISTREAM);
      io::Serialize::Read( m_MovingPoint, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CyclicCoordinateDescent::TargetAndMovingPointPair::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // Write members
      io::Serialize::Write( m_TargetPoint, OSTREAM, INDENT);
      io::Serialize::Write( m_MovingPoint, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CyclicCoordinateDescent::ResultsAndCoefficients::Read( std::istream &ISTREAM)
    {
      // Read members
      io::Serialize::Read( m_OptimalRotation, ISTREAM);
      io::Serialize::Read( m_MinimumDistanceSum, ISTREAM);
      io::Serialize::Read( m_CoefficientA, ISTREAM);
      io::Serialize::Read( m_CoefficientB, ISTREAM);
      io::Serialize::Read( m_CoefficientC, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CyclicCoordinateDescent::ResultsAndCoefficients::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // Write members
      io::Serialize::Write( m_OptimalRotation, OSTREAM, INDENT);
      io::Serialize::Write( m_MinimumDistanceSum, OSTREAM, INDENT);
      io::Serialize::Write( m_CoefficientA, OSTREAM, INDENT);
      io::Serialize::Write( m_CoefficientB, OSTREAM, INDENT);
      io::Serialize::Write( m_CoefficientC, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CyclicCoordinateDescent::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CyclicCoordinateDescent::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief CalculateDistanceSum calculates the sum distance difference between moving and target points
    //! The sum distance difference between the moving points and their corresponding fixed target points will result
    //! if the moving points are rotated by the given rotation angle around the rotation axis used when the
    //! coefficients were calculated.
    //! @param ROTATION the amount of rotation that should be considered
    //! @param COEFFICIENT_A please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_B please see the reference given above (Canutescu et al. Protein Sci, 2003)
    //! @param COEFFICIENT_C please see the reference given above (Canutescu et al. Protein Sci, 2003)
    double CyclicCoordinateDescent::CalculateDistanceSum
    (
      const double ROTATION, const double COEFFICIENT_A, const double COEFFICIENT_B, const double COEFFICIENT_C
    ) const
    {
      // calculate the distance sum as defined by Canutescu et al. Protein Sci, 2003
      const double distance_sum( COEFFICIENT_A - COEFFICIENT_B * cos( ROTATION) - COEFFICIENT_C * sin( ROTATION));

      // return distance sum
      return distance_sum;
    }

    //! @brief CalculateCurrentCoefficientValues calculates the contribution to the total coefficient values for a
    //! single target and moving point pair
    //! @param ROTATION_AXIS the axis around which rotation will occur
    //! @param CURRENT_COORDINATES the pair of target and moving points
    //! @return a VectorND< 3, double> which holds coefficients a, b, and c, respectively
    storage::VectorND< 3, double> CyclicCoordinateDescent::CalculateCurrentCoefficientValues
    (
      const LineSegment3D &ROTATION_AXIS,
      const CyclicCoordinateDescent::TargetAndMovingPointPair &CURRENT_COORDINATES
    ) const
    {
      // O, the footpoint of the moving atoms along the rotation axis
      const linal::Vector3D footpoint
      (
        linal::CalculateFootpoint
        (
          CURRENT_COORDINATES.GetMovingPoint(), ROTATION_AXIS.GetStartPoint(), ROTATION_AXIS.GetDirection()
        )
      );

      // message the footpoint coordinates
      BCL_MessageDbg( "footprint coordinates are " + util::Format()( footpoint));

      // o^, the unit vector of the rotation axis vector
      linal::Vector3D rotation_axis_unit_vector( ROTATION_AXIS.GetDirection());
      rotation_axis_unit_vector.Normalize();

      // message o^
      BCL_MessageDbg( "rotation_axis_unit_vector " + util::Format()( rotation_axis_unit_vector));

      // f, the distance from O to F
      const double footpoint_to_target_distance( linal::Distance( CURRENT_COORDINATES.GetTargetPoint(), footpoint));

      // message f
      BCL_MessageDbg
      (
        "footpoint_to_target_distance " + util::Format()( footpoint_to_target_distance)
      );

      // f->, the vector from O->F
      const linal::Vector3D footpoint_to_target_vector( CURRENT_COORDINATES.GetTargetPoint() - footpoint);

      // message f->
      BCL_MessageDbg( "footpoint_to_target_vector " + util::Format()( footpoint_to_target_vector));

      // f^, the unit vector of f->
      linal::Vector3D footpoint_to_target_unit_vector( footpoint_to_target_vector);
      footpoint_to_target_unit_vector.Normalize();

      // message f^
      BCL_MessageDbg
      (
        "footpoint_to_target_unit_vector " + util::Format()( footpoint_to_target_unit_vector)
      );

      // r, the distance from O to M
      const double footpoint_to_moving_distance( linal::Distance( CURRENT_COORDINATES.GetMovingPoint(), footpoint));

      // message r
      BCL_MessageDbg
      (
        "footpoint_to_moving_distance " + util::Format()( footpoint_to_moving_distance)
      );

      // r->, the vector from O->M
      const linal::Vector3D footpoint_to_moving_vector( CURRENT_COORDINATES.GetMovingPoint() - footpoint);

      // message r->
      BCL_MessageDbg( "footpoint_to_moving_vector " + util::Format()( footpoint_to_moving_vector));

      // r^, the unit vector of r->
      linal::Vector3D footpoint_to_moving_unit_vector( footpoint_to_moving_vector);
      footpoint_to_moving_unit_vector.Normalize();

      // message r^
      BCL_MessageDbg
      (
        "footpoint_to_moving_unit_vector " + util::Format()( footpoint_to_moving_unit_vector)
      );

      // s, the third orthogonal coordinate axis in the coordinate system defined by and o^, r^
      const linal::Vector3D third_coordinate_axis
      (
        linal::CrossProduct( footpoint_to_moving_unit_vector, rotation_axis_unit_vector)
      );

      // message s
      BCL_MessageDbg( "third_coordinate_axis " + util::Format()( third_coordinate_axis));

      // s^, the unit vector of s
      linal::Vector3D third_coordinate_axis_unit_vector( third_coordinate_axis);
      third_coordinate_axis_unit_vector.Normalize();

      // calculate coefficient a
      const double coefficient_a
      (
        math::Sqr( footpoint_to_moving_distance) + math::Sqr( footpoint_to_target_distance)
      );

      // calculate coefficient b
      const double coefficient_b
      (
        2 * footpoint_to_moving_distance * linal::ScalarProduct
        (
          footpoint_to_target_vector, footpoint_to_moving_unit_vector
        )
      );

      // calculate coefficient c
      const double coefficient_c
      (
        2 * footpoint_to_moving_distance * linal::ScalarProduct
        (
          footpoint_to_target_vector, third_coordinate_axis_unit_vector
        )
      );

      // message the three coefficients
      BCL_MessageDbg
      (
        "coefficient a\t" + util::Format()( coefficient_a) + "\tcoefficient b\t"
        + util::Format()( coefficient_b) + "\tcoefficient c\t" + util::Format()( coefficient_c)
      );

      // return the three coefficients
      return storage::VectorND< 3, double>( coefficient_a, coefficient_b, coefficient_c);
    }

    //! @brief operator== for defining the equality between two ResultsAndCoefficients objects
    //! @param RESULT_A the first ResultsAndCoefficients object
    //! @param RESULT_B the second ResultsAndCoefficients object
    //! @return boolean, true if each component of the two ResultsAndCoefficients objects is equal within tolerance
    bool operator==
    (
      const CyclicCoordinateDescent::ResultsAndCoefficients &RESULT_A,
      const CyclicCoordinateDescent::ResultsAndCoefficients &RESULT_B
    )
    {
      // return true if all values are equal within tolerance, false otherwise
      return math::EqualWithinTolerance( RESULT_A.GetOptimalRotation(), RESULT_B.GetOptimalRotation()) &&
        math::EqualWithinTolerance( RESULT_A.GetMinimumDistanceSum(), RESULT_B.GetMinimumDistanceSum()) &&
        math::EqualWithinTolerance( RESULT_A.GetCoefficientA(), RESULT_B.GetCoefficientA()) &&
        math::EqualWithinTolerance( RESULT_A.GetCoefficientB(), RESULT_B.GetCoefficientB()) &&
        math::EqualWithinTolerance( RESULT_A.GetCoefficientC(), RESULT_B.GetCoefficientC());
    }

    //! @brief operator== for defining the equality between two TargetAndMovingPointPair objects
    //! @param PAIR_A the first TargetAndMovingPointPair object
    //! @param PAIR_B the second TargetAndMovingPointPair object
    //! @return boolean, true if the coordinates of the two TargetAndMovingPointPair objects is equal within tolerance
    bool operator==
    (
      const CyclicCoordinateDescent::TargetAndMovingPointPair &PAIR_A,
      const CyclicCoordinateDescent::TargetAndMovingPointPair &PAIR_B
    )
    {
      // return true if the coordinates of the target and moving points in each TargetAndMovingPointPair is equal
      // within tolerance
      return math::EqualWithinTolerance( PAIR_A.GetTargetPoint(), PAIR_B.GetTargetPoint()) &&
        math::EqualWithinTolerance( PAIR_A.GetMovingPoint(), PAIR_B.GetMovingPoint());
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_cylinder_coordinates.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CylinderCoordinates::s_Instance
    (
      GetObjectInstances().AddInstance( new CylinderCoordinates())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CylinderCoordinates::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief the default constructor
    CylinderCoordinates::CylinderCoordinates() :
      m_Height( double( 0.0)),
      m_Radius( double( 0.0)),
      m_Angle( double( 0.0))
    {
    }

    //! @brief construct CylinderCoordinates from sizes
    //! @param HEIGHT height
    //! @param RADIUS radius
    //! @param ANGLE angle
    CylinderCoordinates::CylinderCoordinates( const double HEIGHT, const double RADIUS, const double ANGLE) :
      m_Height( HEIGHT),
      m_Radius( RADIUS),
      m_Angle( ANGLE)
    {
    }

    //! @brief copy constructor
    //! @return pointer to a new CylinderCoordinates copied from this one
    CylinderCoordinates *CylinderCoordinates::Clone() const
    {
      return new CylinderCoordinates( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate cartesian coordinates from cylinderical coordinates
    //! @return Cartesian coordinates generated from cylinderical coordinates
    linal::Vector3D CylinderCoordinates::GetCartesianCoordinates() const
    {
      return linal::Vector3D
      (
        GetRadius() * cos( GetAngle()),
        GetRadius() * sin( GetAngle()),
        GetHeight()
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read CylinderCoordinates from std::istream
    std::istream &CylinderCoordinates::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Height, ISTREAM);
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_Angle, ISTREAM);

      // end
      return ISTREAM;
    }

    //! write CylinderCoordinates into std::ostream
    std::ostream &CylinderCoordinates::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Height, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Radius, OSTREAM) << '\t';
      io::Serialize::Write( m_Angle, OSTREAM);

      // end
      return OSTREAM;
     }

  } // namespace coord
} // namespace bcl

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
#include "coord/bcl_coord_geometric_hashing.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_geometric_hash_storage_hash_map.h"
#include "coord/bcl_coord_point_cloud.h"
#include "density/bcl_density_mask_3d.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_tertiary_function_job_with_data.h"
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
    const util::SiPtr< const util::ObjectInterface> GeometricHashing::s_Instance
    (
      GetObjectInstances().AddInstance( new GeometricHashing())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //!default Constructor
    GeometricHashing::GeometricHashing() :
      m_HashMap( new GeometricHashStorageHashMap())
    {
    }

    //! construct Hashmap from GeometricHashStorage
    GeometricHashing::GeometricHashing( const util::ShPtr< GeometricHashStorageInterface> &STORAGE_SP) :
      m_HashMap( STORAGE_SP)
    {
    }

    //! virtual copy constructor
    GeometricHashing *GeometricHashing::Clone() const
    {
      return new GeometricHashing( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GeometricHashing::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! return pairwise distances and order them in three different sizes
    storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > >
    GeometricHashing::CalculatePointPairs
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const storage::VectorND< 4, double> &THRESHOLD
    )
    {
      //make the vectors as long as the number of given coordinates
      storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > >
      pointpairs;

      storage::Vector< size_t> count( 3);
      // 0: point_neighbors within radius
      // 1: short_point_pairs between lower threshold and lower_threshold + 1/3 * (upper threshold - lower threshold)
      // 2: middle_point_pairs between lower lower_threshold + 1/3 * (upper threshold - lower threshold)  and lower threshold + 2/3 * (upper threshold - lower threshold)
      // 3: long_point_pairs lower threshold + 2/3 * (upper threshold - lower threshold) and upper threshold

//      const storage::VectorND< 4, double> thresholds( CreateEqualOccupiedIntervals( CalculateDistanceDistribution( COORDINATES, THRESHOLD, 30)));
      const double lower_tresh2( math::Sqr( THRESHOLD.First()));
      const double low2(         math::Sqr( THRESHOLD.Second()));
      const double middle2(      math::Sqr( THRESHOLD.Third()));
      const double upper_tresh2( math::Sqr( THRESHOLD.Fourth()));

      const double maxdist( THRESHOLD.Fourth());

      // calculate pair distances and look, whether they are within the threshold and store all pairs of points within this threshold
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr_1( COORDINATES.Begin()),
          coord_itr_end( COORDINATES.End());
        coord_itr_1 != coord_itr_end;
        ++coord_itr_1
      )
      {
        // references to list of shortest, middle and longest connection for easier access
        storage::Set< util::SiPtr< const linal::Vector3D> > &current_short_list( pointpairs.First()[ *coord_itr_1]);
        storage::Set< util::SiPtr< const linal::Vector3D> > &current_middle_list( pointpairs.Second()[ *coord_itr_1]);
        storage::Set< util::SiPtr< const linal::Vector3D> > &current_long_list( pointpairs.Third()[ *coord_itr_1]);
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator
            coord_itr_2( coord_itr_1 + 1);
          coord_itr_2 != coord_itr_end;
          ++coord_itr_2
        )
        {
          const linal::Vector3D delta_xyz( linal::AbsoluteDifference( **coord_itr_1, **coord_itr_2));

          // if either x, y or z coordinate difference is larger that the maximal distance, go to next coordinate
          if( delta_xyz.X() > maxdist || delta_xyz.Y() > maxdist || delta_xyz.Z() > maxdist)
          {
            continue;
          }

          // check the actual (square) distance - if too big or too small, skip
          const double dist2( delta_xyz.SquareNorm());
          if( dist2 >= upper_tresh2 || dist2 <= lower_tresh2)
          {
            continue;
          }
          // find the interval it belongs to
          else if( dist2 < low2)
          {
            // short_point_pairs
            current_short_list.InsertElement( *coord_itr_2);
            pointpairs.First()[ *coord_itr_2].InsertElement( *coord_itr_1);
            count( 0)++;
          }
          else if( dist2 >= middle2)
          {
            // long_point_pairs
            current_long_list.InsertElement( *coord_itr_2);
            pointpairs.Third()[ *coord_itr_2].InsertElement( *coord_itr_1);
            count( 2)++;
          }
          else
          {
            // middle_point_pairs
            current_middle_list.InsertElement( *coord_itr_2);
            pointpairs.Second()[ *coord_itr_2].InsertElement( *coord_itr_1);
            count( 1)++;
          }
        }
      }

      // tells user that calculating pairwise distances has finished
      BCL_MessageStd
      (
        std::string( "CalculatePointPairs( COORDINATES) calculate pairwise distance has finshed\n") +
        "number of pairs of points within lower distance: "  + util::Format()( count( 0)) + "\n"
        "number of pairs of points within middle distance: " + util::Format()( count( 1)) + "\n"
        "number of pairs of points within higher distance: " + util::Format()( count( 2))
      );

      // end
      return pointpairs;
    }

    //! @brief determine triangular bases given a map of all point pairs
    //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
    //! @return a list of possible triangles
    storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > >
    GeometricHashing::DeterminePossibleBases
    (
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS
    )
    {
      // store all possible bases for a random search
      storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > > bases;

      // loop over all pairs to find bases within threshold
      for
      (
        storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
          first_itr( POINT_PAIRS.Third().Begin()),
          first_itr_end( POINT_PAIRS.Third().End());
        first_itr != first_itr_end;
        ++first_itr
      )
      {
        for
        (
          storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
            itr0( first_itr->second.Begin()), itr0_end( first_itr->second.End());
          itr0 != itr0_end; ++itr0
        )
        {
          storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
            current_map_itr1( POINT_PAIRS.Second().Find( *itr0));

          for
          (
            storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
              itr1( current_map_itr1->second.Begin()), itr1_end( current_map_itr1->second.End());
            itr1 != itr1_end;
            ++itr1)
          {
            storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
              current_map_itr2( POINT_PAIRS.First().Find( *itr1));

            // serach a matching point to close the triangular base
            const storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
              itr2( std::find( current_map_itr2->second.Begin(), current_map_itr2->second.End(), first_itr->first));

            // check that triangular base is closed
            if( itr2 != current_map_itr2->second.End())
            {
              bases.PushBack( storage::VectorND< 3, util::SiPtr< const linal::Vector3D> >( *itr0, *itr1, *itr2));
            }
          }
        }
      }

      // end
      return bases;
    }

    //! @brief pick a subset of random bases from a list of triangular bases which are equally distributed over different distances from the given center
    //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
    //! @param NUMBER size of the resulting set
    //! @param CENTER the center to which the bases should be distributed equally distant
    //! @param NUMBER_BINS number of bin to distribute bases equally over
    //! @return List of transformations and their centers
    storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
    GeometricHashing::GatherEquallyDistributedBases
    (
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS,
      const size_t NUMBER,
      const linal::Vector3D &CENTER,
      const size_t NUMBER_BINS
    )
    {
      BCL_MessageStd( "determine all transformations and centers for all possible bases!");

      // calculate all bases and their centers
      const storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
        bases_center( ConvertBasesToTransformationsAndCenter( DeterminePossibleBases( POINT_PAIRS)));

      BCL_MessageStd( "number of possible bases: " + util::Format()( bases_center.GetSize()));

      // store the distance from the center with the base and its center
      BCL_MessageStd( "determine the distance of all bases from given center");
      storage::List< storage::Pair< double, util::SiPtr< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > > > distance_base_center;

      for
      (
        storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >::const_iterator
          itr( bases_center.Begin()), itr_end( bases_center.End());
        itr != itr_end;
        ++itr
      )
      {
        distance_base_center.PushBack
        (
          storage::Pair< double, util::SiPtr< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > >
          (
            ( CENTER - itr->Second()).Norm(),
            util::SiPtr< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >( &( *itr))
          )
        );
      }

      BCL_MessageStd( "sort all bases by their distance form the center");
      // sort by distance
      distance_base_center.Sort
      (
        storage::PairBinaryPredicateFirst< double, util::SiPtr< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > >
        (
          **math::Comparisons< double>::GetEnums().e_Less
        )
      );

      // acquire the shortest and largest distance
      const storage::VectorND< 2, double> short_long
      (
        distance_base_center.FirstElement().First(), distance_base_center.LastElement().First()
      );

      // calculate the binsize
      const double bin_size( ( short_long.Second() - short_long.First()) / NUMBER_BINS);

      BCL_MessageStd
      (
        "bin all bases by their distance in " + util::Format()( NUMBER) + " bins in a range: " +
        util::Format()( short_long)
      );
      // bin all bases into that map
      storage::Map< size_t, util::SiPtrList< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > >
        binned_bases;

      for
      (
        storage::List< storage::Pair< double, util::SiPtr< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > > >::const_iterator
          itr( distance_base_center.Begin()), itr_end( distance_base_center.End());
        itr != itr_end;
        ++itr
      )
      {
        binned_bases[ size_t( itr->First() / bin_size)].PushBack( itr->Second());
      }

      BCL_MessageStd
      (
        "acquire random bases from all bins"
      );
      // instantiate the bases
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > final_bases_center;

      for
      (
        storage::Map< size_t, util::SiPtrList< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > >::iterator
          bin_itr( binned_bases.Begin()), bin_itr_end( binned_bases.End());
        bin_itr != bin_itr_end && final_bases_center.GetSize() < NUMBER;
        ++bin_itr
      )
      {
        util::SiPtrList< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > &current_list( bin_itr->second);

        // iterate to the maximal number, or the number of available bases
        for( size_t i( 0), i_end( std::min( NUMBER, current_list.GetSize())); i < i_end; ++i)
        {
          util::SiPtrList< const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >::iterator current_base_itr
          (
            random::GetGlobalRandom().Iterator( current_list.Begin(), current_list.End(), current_list.GetSize())
          );

          // insert the transformation derived from the triangular base
          final_bases_center.PushBack( **current_base_itr);

          // remove the iterator
          current_list.Remove( current_base_itr);
        }
      }

      BCL_MessageStd
      (
        "list containing equally distributed bases contains " + util::Format()( final_bases_center.GetSize()) +
        " bases"
      );

      // end
      return final_bases_center;
    }

    //! @brief pick a subset of random bases from a list of triangular bases which are equally distributed over different distances from the given center
    //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
    //! @param NUMBER size of the resulting set
    //! @param CENTER the center to which the bases should be distributed equally distant
    //! @param NUMBER_BINS number of bin to distribute bases equally over
    //! @return List of transformations and their centers
    storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
    GeometricHashing::GatherEquallyDistributedBases3D
    (
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS,
      const size_t NUMBER,
      const linal::Vector3D &CENTER,
      const size_t NUMBER_BINS
    ) const
    {
      BCL_MessageStd( "determine all transformations and centers for all possible bases!");

      // calculate all bases and their centers
      const storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
        bases_center( ConvertBasesToTransformationsAndCenter( DeterminePossibleBases( POINT_PAIRS)));

      BCL_MessageStd( "number of possible bases: " + util::Format()( bases_center.GetSize()));

      // collect all centers
      util::SiPtrVector< const linal::Vector3D> centers;
      for
      (
        storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >::const_iterator
          itr( bases_center.Begin()), itr_end( bases_center.End());
        itr != itr_end;
        ++itr
      )
      {
        centers.PushBack( util::ToSiPtr( itr->Second()));
      }

      const storage::VectorND< 2, linal::Vector3D> grid_corners( density::Mask3d::DetermineGridCorners( centers));
      const double feature_distance( m_HashMap->GetMinimalDistance());
      double grid_resolution( feature_distance);

      // store the bases in a grid
      BCL_MessageStd( "associate all transformations with points on a grid");
      storage::Map< size_t, storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > > bases_on_grid;

      while( bases_on_grid.GetSize() < NUMBER && grid_resolution > feature_distance * 0.2)
      {
        bases_on_grid.Reset();
        for
        (
          storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >::const_iterator
            itr( bases_center.Begin()), itr_end( bases_center.End());
          itr != itr_end;
          ++itr
        )
        {
          linal::Vector3D grid_position( itr->Second() - grid_corners.First());
          grid_position /= feature_distance;

          const size_t key( size_t( grid_position( 0)) * 1000000 + size_t( grid_position( 1)) * 1000 + size_t( grid_position( 2)));

          bases_on_grid[ key].PushBack( *itr);
        }
        grid_resolution *= 0.9;
      }

      BCL_MessageStd( "number of occupied grid elements: " + util::Format()( bases_on_grid.GetSize()));

      // instantiate the bases
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > final_bases;

      while( final_bases.GetSize() < NUMBER)
      {
        for
        (
          storage::Map< size_t, storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > >::iterator
            grid_itr( bases_on_grid.Begin()), grid_itr_end( bases_on_grid.End());
          grid_itr != grid_itr_end;
          ++grid_itr
        )
        {
          storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >::iterator current_itr
          (
            random::GetGlobalRandom().Iterator( grid_itr->second.Begin(), grid_itr->second.End(), grid_itr->second.GetSize())
          );

          // insert the transformation derived from the triangular base
          final_bases.PushBack( *current_itr);
        }
      }

      BCL_MessageStd
      (
        "list containing equally distributed bases contains " + util::Format()( final_bases.GetSize()) +
        " bases"
      );

      while( final_bases.GetSize() > NUMBER)
      {
        final_bases.Remove( random::GetGlobalRandom().Iterator( final_bases.Begin(), final_bases.End(), final_bases.GetSize()));
      }

      // end
      return final_bases;
    }

    //! @brief convert a list of triangular bases into transformations and their triangular centers
    storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
    GeometricHashing::ConvertBasesToTransformationsAndCenter
    (
      const storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > > &BASES
    )
    {
      // a randomly selected list of bases
      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > trans_center;

      for
      (
        storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > >::const_iterator
          current_base_points_itr( BASES.Begin()), current_base_points_itr_end( BASES.End());
        current_base_points_itr != current_base_points_itr_end;
        ++current_base_points_itr
      )
      {
        // insert the transformation derived from the triangular base
        trans_center.PushBack
        (
          CalculateBasisTransformationAndCenter
          (
            *current_base_points_itr->First(),
            *current_base_points_itr->Second(),
            *current_base_points_itr->Third()
          )
        );
      }

      // end
      return trans_center;
    }

    //! @brief for a set of coordinates, return the maximal number of bases
    //! @param COORDINATES the coordinates the number of bases is determined for
    //! @param THRESHOLD the threshold used for the triangular base
    size_t
    GeometricHashing::NumberPossibleBases
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const storage::VectorND< 4, double> &THRESHOLD
    )
    {
      BCL_Assert( COORDINATES.GetSize() > 2, "needs at least 3 points");

      //calculate  pairwise distances
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > >
        pointpairs( CalculatePointPairs( COORDINATES, THRESHOLD));

      // determine the number of possible bases
      return NumberPossibleBases( pointpairs);
    }

    //! @brief for pointpairs, return the maximal number of bases
    //! @param POINT_PAIRS point pairs for lower, middle an upper threshold
    //! @return number of bases which fullfill the given THRESHOLD
    size_t
    GeometricHashing::NumberPossibleBases
    (
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > > &POINT_PAIRS
    )
    {
      size_t number_bases( 0);
      // loop over all pairs to find bases within threshold
      // point1
      for
      (
        storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
          first_itr( POINT_PAIRS.Third().Begin()),
          first_itr_end( POINT_PAIRS.Third().End());
        first_itr != first_itr_end;
        ++first_itr
      )
      {
        // line point1->point2
        for
        (
          storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
            itr0( first_itr->second.Begin()), itr0_end( first_itr->second.End());
          itr0 != itr0_end; ++itr0
        )
        {
          // set all point3s connecting to point2
          const storage::Set< util::SiPtr< const linal::Vector3D> >
            &current_map1( POINT_PAIRS.Second().Find( *itr0)->second);
          // point2->point3
          for
          (
            storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
              itr1( current_map1.Begin()), itr1_end( current_map1.End());
            itr1 != itr1_end;
            ++itr1
          )
          {
            const storage::Set< util::SiPtr< const linal::Vector3D> >
              &current_map2( POINT_PAIRS.First().Find( *itr1)->second);

            // search a matching point to close the triangular base
            const storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
              itr2( std::find( current_map2.Begin(), current_map2.End(), first_itr->first));

            // check that triangular base is closed
            if( itr2 != current_map2.End())
            {
              ++number_bases;
              break;
            }
          }
        }
      }

      return number_bases;
    }

    //! build Hashmap from PointCloud and lower and upper Threshold for the basis - by default m_Threshold
    void GeometricHashing::BuildHash( const PointCloud &COORDINATES)
    {
      // if hash map was already calculated, return
      if( m_HashMap->IsReadOnly())
      {
        BCL_MessageStd( "skip building since hash map already exists");
        return;
      }

      BCL_Assert( COORDINATES.GetSize() > 2, "needs at least 3 points");
      const storage::VectorND< 4, double> threshold( m_HashMap->GetThreshold());
      const double radius( m_HashMap->GetRadius());
      const util::SiPtrVector< const linal::Vector3D> siptr_vec_coord( util::ConvertToConstSiPtrVector( COORDINATES.GetData()));

      // calculate  pairwise distances
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > >
        pointpairs( CalculatePointPairs( siptr_vec_coord, threshold));

      // calculate the number of possible bases
      const size_t number_long_distances( pointpairs.Third().GetSize());
      const size_t msg_fraction( 10);
      const size_t count_long_msg( number_long_distances / msg_fraction);
      size_t long_distance_count( 0);

      BCL_MessageStd( "building hash with nr long distances: " + util::Format()( number_long_distances));

      // loop over all pairs to find bases within threshold - iterate over list of long triangle sides
      for
      (
        storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
          first_itr( pointpairs.Third().Begin()),
          first_itr_end( pointpairs.Third().End());
        first_itr != first_itr_end;
        ++first_itr
      )
      {
        // iterate over all points that are connected by a long side
        for
        (
          storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
            itr0( first_itr->second.Begin()), itr0_end( first_itr->second.End());
          itr0 != itr0_end;
          ++itr0
        )
        {
          // get the list of points that are connected by a middle length distance
          storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
            current_map_itr1( pointpairs.Second().Find( *itr0));

          // iterate over all points connected to itr0 by a middle length side
          for
          (
            storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
              itr1( current_map_itr1->second.Begin()), itr1_end( current_map_itr1->second.End());
            itr1 != itr1_end;
            ++itr1)
          {
            // get the list of points that are connected by a short length distance
            storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > >::const_iterator
              current_map_itr2( pointpairs.First().Find( *itr1));

            // search a matching point to close the triangular base
            const storage::Set< util::SiPtr< const linal::Vector3D> >::const_iterator
              itr2( std::find( current_map_itr2->second.Begin(), current_map_itr2->second.End(), first_itr->first));

            // check that triangular base is closed
            if( itr2 == current_map_itr2->second.End())
            {
              // triangle is not closed - try next point
              continue;
            }

            // build transformation matrix and center
            const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> trans_center
            (
              // itr0 is connected to itr1 by long side to itr2 by short side
              // itr1 is connected to itr2 by middle side to itr0 by long side
              // itr2 is connected to itr0 by short side to itr1 by middle side
              CalculateBasisTransformationAndCenter( **itr0, **itr1, **itr2)
            );

            // generate a list of neighbors within radius for the center of the base
            const util::SiPtrVector< const linal::Vector3D> neighbors
            (
              DetermineNeighbors( trans_center.Second(), siptr_vec_coord, radius)
            );

            // store hash keys with Transformation matrix in the HashMap
            m_HashMap->Store( neighbors, trans_center.First());

            // triangle was closed
            break;
          }
        }

        // message every x percent
        ++long_distance_count;
        if( long_distance_count % count_long_msg == 0)
        {
          BCL_MessageStd
          (
            "building hash " + util::Format()( long_distance_count / count_long_msg * msg_fraction) + "% done"
          );
        }
      }

      // finalize hash map - depending on the implementation of the hash map, if program breaks before that point
      // the hash map knows, that it was not finished, and will clean itself
      m_HashMap->Finalize();

      BCL_MessageStd
      (
        "building hash has finished\nNumber of Hash keys: " + util::Format()( m_HashMap->NumberOfHashkeys()) +
        "\tNumber of TransformationMatrices: " + util::Format()( m_HashMap->NumberOfTransformationMatrices())
      );
    }

    //! @brief search hash map and return storage list of pairs of Transformation matrices for the highest counts
    //! @param COORDINATES the coordinates that will be searched in the hash map
    //! @param SAVE_BEST the maximal number of possible transformations returned
    //! @param SEARCH_TRIALS is the number of random searched bases, SAVEBEST is the number of transformationmatrices which are stored
    //! @param DIFFERENCE_ROT_TRANS has two values for the max deviation of a new found transformation in rotation and translation in comparison to already found transformations
    storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >
    GeometricHashing::SearchTarget
    (
      const PointCloud &COORDINATES,
      const size_t &SAVE_BEST,
      const size_t &SEARCH_TRIALS,
      const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
    ) const
    {
      BCL_Assert( m_HashMap->IsReadOnly(), "cannot fit with a map, that is not completely built");

      // create a siptr vector on the point cloud
      const util::SiPtrVector< const linal::Vector3D> siptr_vec_coord( util::ConvertToConstSiPtrVector( COORDINATES.GetData()));

      // copy of the threshold
      const storage::VectorND< 4, double> threshold( m_HashMap->GetThreshold());

      // tells user that the building of the geometric hash map starts
      BCL_MessageStd( "start of fitting with threshold of " + util::Format()( threshold));

      // gather all possible bases with their distance to the position of the center of those bases
      const storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > random_bases_center
      (
        GatherEquallyDistributedBases3D
        (
          CalculatePointPairs( siptr_vec_coord, threshold),
          SEARCH_TRIALS,
          COORDINATES.GetCenter(),
          200
        )
      );

      // check that there are any bases in the target
      if( random_bases_center.IsEmpty())
      {
         BCL_MessageCrt( "There are no possible bases in this target!");

         // return an empty list
         return storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >();
      }

      // report the number of possible bases in target
      BCL_MessageStd
      (
        "There are " + util::Format()( random_bases_center.GetSize()) + " possible bases to be searched!"
      );

      // number of random bases
      const size_t number_random_bases( random_bases_center.GetSize());

      // issue a warning when there are less possible bases than SEARCHTRIALS
      if( random_bases_center.GetSize() < SEARCH_TRIALS)
      {
        BCL_MessageCrt
        (
          "There are less possible bases in this target than user defined search trials! " +
          util::Format()( number_random_bases) + " < " + util::Format()( SEARCH_TRIALS)
        );
      }

      // map to store jobs an their result
      storage::List
      <
        storage::Pair
        <
          util::ShPtr< sched::JobInterface>,
          storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
        >
      >
      jobs_result;

      // iterate over all randomly selected bases
      for
      (
        storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >::const_iterator
          current_base_itr( random_bases_center.Begin()), current_base_itr_end( random_bases_center.End());
        current_base_itr != current_base_itr_end;
        ++current_base_itr
      )
      {
        // create new job trans_center result pair
        jobs_result.PushBack
        (
          storage::Pair
          <
            util::ShPtr< sched::JobInterface>,
            storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
          >()
        );

        // create reference to the just create triplet
        storage::Pair
        <
          util::ShPtr< sched::JobInterface>,
          storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
        > &current_job_result( jobs_result.LastElement());

        // initialize the actual job object
        current_job_result.First() =
        util::ShPtr< sched::JobInterface>
        (
          new sched::TertiaryFunctionJobWithData
          <
            const storage::Pair< math::TransformationMatrix3D, linal::Vector3D>,
            const util::SiPtrVector< const linal::Vector3D>,
            const size_t,
            storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >,
            GeometricHashing
          >
          (
            0,
            *this,
            &GeometricHashing::ReturnBestCounts,
            *current_base_itr,
            siptr_vec_coord,
            SAVE_BEST,
            sched::JobInterface::e_READY,
            &current_job_result.Second()
          )
        );

        // submit job to the scheduler
        sched::GetScheduler().RunJob( current_job_result.First());

        if( jobs_result.GetSize() % ( number_random_bases / 10) == 1)
        {
          BCL_MessageStd
          (
            "starting to process base " + util::Format()( jobs_result.GetSize()) + " of " + util::Format()( number_random_bases) + " bases for fitting!"
          );
        }
      }

      // instantiate storage::Vector storing all best transformations
      storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> > all_best_transformations;

      // join all threads
      for
      (
        storage::List
        <
          storage::Pair
          <
            util::ShPtr< sched::JobInterface>,
            storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
          >
        >::iterator
          job_result_itr( jobs_result.Begin()), job_result_itr_end( jobs_result.End());
        job_result_itr != job_result_itr_end;
        ++job_result_itr
      )
      {
        sched::GetScheduler().Join( job_result_itr->First());
        all_best_transformations.Append( job_result_itr->Second());
      }

      // filter similar transformations and save these with the highest counts
      return FilterSimilarTransformations( all_best_transformations, SAVE_BEST, DIFFERENCE_ROT_TRANS);
    }

    //! @brief returns all possible bases within the given pointcloud
    //! @param COORDINATES the coordinates in the hash map
    storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> >
    GeometricHashing::ReturnPossibleBases
    (
      const PointCloud &COORDINATES
    ) const
    {
      // create a siptr vector on the point cloud
      const util::SiPtrVector< const linal::Vector3D> siptr_vec_coord( util::ConvertToConstSiPtrVector( COORDINATES.GetData()));

      // copy of the threshold
      const storage::VectorND< 4, double> threshold( m_HashMap->GetThreshold());

      // calculate  pairwise distances
      const storage::VectorND< 3, storage::Map< util::SiPtr< const linal::Vector3D>, storage::Set< util::SiPtr< const linal::Vector3D> > > >
        pointpairs( CalculatePointPairs( siptr_vec_coord, threshold));

      // store all possible bases
      storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > >
        bases( DeterminePossibleBases( pointpairs));

      storage::List< storage::Pair< math::TransformationMatrix3D, linal::Vector3D> > possible_bases;

      // check that there are any bases in the target
      if( bases.IsEmpty())
      {
         BCL_MessageCrt( "There are no possible bases in this target!");

         return possible_bases;
      }

      // iterate over all triangles
      for
      (
        storage::List< storage::VectorND< 3, util::SiPtr< const linal::Vector3D> > >::iterator
          current_base_points_itr( bases.Begin()), current_base_points_itr_end( bases.End());
        current_base_points_itr != current_base_points_itr_end;
        ++current_base_points_itr
      )
      {
        // insert the transformation derived from the triangular base
        possible_bases.PushBack
        (
          CalculateBasisTransformationAndCenter
          (
            *current_base_points_itr->First(),
            *current_base_points_itr->Second(),
            *current_base_points_itr->Third()
          )
        );
      }

      // end
      return possible_bases;
    }

    //! @brief return the best transformations and their hashscore for one base within a set of coordinates
    //! @param TRANSFORMATION_CENTER the base and center to be considered within the COORDINATES
    //! @param COORDINATES the set of coordinates
    //! @param SAVE_BEST the maximal number of results to be returned
    //! @return a list of Transformations and their hashscore
    storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
    GeometricHashing::ReturnBestCounts
    (
      const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> &TRANSFORMATION_CENTER,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const size_t &SAVE_BEST
    ) const
    {
      const double radius( m_HashMap->GetRadius());
      const util::SiPtrVector< const linal::Vector3D> neighbors( DetermineNeighbors( TRANSFORMATION_CENTER.Second(), COORDINATES, radius));

      // search the m_HashMap for the best transformations
      storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
        best_transformations( m_HashMap->ReturnBestCounts( neighbors, TRANSFORMATION_CENTER.First(), SAVE_BEST));

      // apply the current base transformation to the best transformations
      for
      (
        storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >::iterator
          best_tm_itr( best_transformations.Begin()), best_tm_itr_end( best_transformations.End());
        best_tm_itr != best_tm_itr_end;
        ++best_tm_itr
      )
      {
        ( *best_tm_itr->First()) = math::TransformationMatrix3D( TRANSFORMATION_CENTER.First())( math::Inverse( *best_tm_itr->First()));
      }

      // return the best transformations for the current base transformation
      return best_transformations;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write HashMap to std::ostream using the given util::Format
    std::ostream &GeometricHashing::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_HashMap, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

    //! read HashMap from io::IFStream
    std::istream &GeometricHashing::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_HashMap, ISTREAM);

      // return
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate Transformation, so that 3 Vectors are the Basis for a coordinate System and their center
    //! @param VECTOR_1 is connected to VECTOR_2 by long side to VECTOR_3 by short side
    //! @param VECTOR_2 is connected to VECTOR_3 by middle side to VECTOR_1 by long side
    //! @param VECTOR_3 is connected to VECTOR_1 by short side to VECTOR_2 by middle side
    //! @return the Transformation matrix, that puts center of triangle in origin, Vector2 on x-axis and vector 3 in xy-plane and the center of the untransformed vectors
    storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
    GeometricHashing::CalculateBasisTransformationAndCenter
    (
      const linal::Vector3D &VECTOR_1,
      const linal::Vector3D &VECTOR_2,
      const linal::Vector3D &VECTOR_3
    )
    {
      // initialize
      linal::Vector3D data[ 3] = { VECTOR_1, VECTOR_2, VECTOR_3};

      util::SiPtrVector< linal::Vector3D> vect( 3, data);

      storage::Pair< math::TransformationMatrix3D, linal::Vector3D> trans_center
      (
        math::TransformationMatrix3D(),                            // default transformation
        ( ( *vect( 0)) + ( *vect( 1)) + ( *vect( 2))) / double( 3) // calculate arithmetic middle
      );

      // move origin to arithmetic middle point of triangular base
      trans_center.First()( -trans_center.Second());
      TransformCoordinates( vect, trans_center.First());

      // rotate vector 2 in XZ-plane
      math::TransformationMatrix3D temp;
      const double alpha_x( atan2( vect( 1)->Y(), vect( 1)->Z()));
      trans_center.First()( GetAxes().e_X, alpha_x);
      temp( GetAxes().e_X, alpha_x);
      TransformCoordinates( vect, temp);

      // rotate vector 2 onto X-axis
      temp.SetUnit();
      const double beta_y( atan2( vect( 1)->Z(), vect( 1)->X()));
      trans_center.First()( GetAxes().e_Y, beta_y);
      temp( GetAxes().e_Y, beta_y);
      TransformCoordinates( vect, temp);

      // rotate vector 3 in XY-plane ( will also rotate vector 1 into XY-plane)
      temp.SetUnit();
      const double gamma_x( -atan2( vect( 2)->Z(), vect( 2)->Y()));
      trans_center.First()( GetAxes().e_X, gamma_x);
      temp( GetAxes().e_X, gamma_x);
      TransformCoordinates( vect, temp);

      // end
      return trans_center;
    }

    //! filter similar TransforamtionMatrices according to DIFFERENCE_ROT_TRANS
    storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >
    GeometricHashing::FilterSimilarTransformations
    (
      const storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> > &BESTTRANSFORMATIONS,
      const size_t &KEEPBEST,
      const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
    )
    {
      // convert input list, into list of siptr to transformation matrix and counts
      std::multiset< std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, GeometricHashStorageHashMap::Sp_CountCompare> sorted_basis_count;
      for
      (
        storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >::const_iterator
          sp_tm_itr( BESTTRANSFORMATIONS.Begin()),
          sp_tm_itr_end( BESTTRANSFORMATIONS.End());
        sp_tm_itr != sp_tm_itr_end; ++sp_tm_itr
      )
      {
        sorted_basis_count.insert
        (
          std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>
          (
            util::SiPtr< const math::TransformationMatrix3D>( sp_tm_itr->First()),
            sp_tm_itr->Second()
          )
        );
      }

      // storage for filtered transformations
      storage::List< storage::Pair< math::TransformationMatrix3D, size_t> > filtered_transformations;

      // iterate from the transformation with the highest count (reverse iterator) over all transformations, till there
      // is no transformation left, or as many transformations are stored as required by the user through the argument KEEPBEST
      for
      (
        std::multiset
        <
          std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>,
          GeometricHashStorageHashMap::Sp_CountCompare
        >::const_reverse_iterator
          set_revitr( sorted_basis_count.rbegin()), set_revitr_end( sorted_basis_count.rend());
        set_revitr != set_revitr_end && filtered_transformations.GetSize() != KEEPBEST;
        ++set_revitr
      )
      {
        bool similar( false);
        for
        (
          storage::List< storage::Pair< math::TransformationMatrix3D, size_t> >::const_iterator
            tm_itr( filtered_transformations.Begin()), tm_itr_end( filtered_transformations.End());
          tm_itr != tm_itr_end;
          ++tm_itr
        )
        {
          math::TransformationMatrix3D diff( tm_itr->First());
          diff( math::Inverse( *set_revitr->first));

          // check for similar Translation and similar rotation
          if
          (
               diff.GetTranslation().Norm() < DIFFERENCE_ROT_TRANS.Second()                // translation
            && diff.GetRotation().EffectiveRotationAngle() < DIFFERENCE_ROT_TRANS.First()  // rotation
          )
          {
            // the matrix from the multiset is similar to an already push backed transformation - will not be stored
            similar = true;
            break;
          }
        }

        // if transformation matrix was similar to already stored transformation get next transformation from the multiset
        if( !similar)
        {
          filtered_transformations.PushBack
          (
            storage::Pair< math::TransformationMatrix3D, size_t>( *set_revitr->first, set_revitr->second)
          );
        }
      }

      // end
      return filtered_transformations;
    }

    //! @brief determine all neighbors within given a radius
    //! @param POINT_OF_INTEREST the point for which all neighbors within the FEATURE_RADIUS is returned
    //! @param COORDINATES the coordinates that are considered
    //! @param FEATURE_RADIUS
    //! @return vector of coordinates within the FEATURE_RADIUS of the POINT_OF_INTEREST
    util::SiPtrVector< const linal::Vector3D>
    GeometricHashing::DetermineNeighbors
    (
      const linal::Vector3D &POINT_OF_INTEREST,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const double FEATURE_RADIUS
    )
    {
      // initialize the vector that stores the neighbors
      util::SiPtrVector< const linal::Vector3D> neighbors;

      // save a square radius, saves time for length of connecting vector comparisons
      const double radius_square( math::Sqr( FEATURE_RADIUS));

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        // skip the point of interest
        if( coord_itr->GetPointer() == &POINT_OF_INTEREST)
        {
          continue;
        }

        // vector that connects the points
        const linal::Vector3D connecting_vector( linal::AbsoluteDifference( POINT_OF_INTEREST, **coord_itr));

        // check if it is within radius
        if
        (
             connecting_vector.X() < FEATURE_RADIUS
          && connecting_vector.Y() < FEATURE_RADIUS
          && connecting_vector.Z() < FEATURE_RADIUS
          && connecting_vector.SquareNorm() < radius_square
        )
        {
          neighbors.PushBack( *coord_itr);
        }
      }

      // end
      return neighbors;
    }

    //! @brief calculate the distance distribution
    //! @param COORDINATES the coordinates that are considered
    //! @param THRESHOLD that determines the triangle
    //! @param NUMBER_BINS number of bins in the histogram generated
    //! @return histogram of for all the distances
    math::Histogram
    GeometricHashing::CalculateDistanceDistribution
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const storage::VectorND< 2, double> &THRESHOLD,
      const size_t NUMBER_BINS
    )
    {
      const double bin_size( ( THRESHOLD.Second() - THRESHOLD.First()) / double( NUMBER_BINS));
      math::Histogram distance_distribution( THRESHOLD.First(), bin_size, NUMBER_BINS);

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr_1( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr_1 != coord_itr_end;
        ++coord_itr_1
      )
      {
        const linal::Vector3D &coord1( **coord_itr_1);
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator
            coord_itr_2( coord_itr_1 + 1);
          coord_itr_2 != coord_itr_end;
          ++coord_itr_2
        )
        {
          distance_distribution.PushBack( linal::Distance( coord1, **coord_itr_2));
        }
      }

      // end
      return distance_distribution;
    }

    //! @brief calculate the distance distribution
    //! @param COORDINATES the coordinates that are considered
    //! @param BIN_SIZE the size of the bins
    //! @param MAX_DISTANCE the maximal distance to be considered
    //! @return histogram of for all the distances
    math::Histogram
    GeometricHashing::CalculateDistanceDistribution
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const double BIN_SIZE,
      const math::Range< double> &DISTANCE_RANGE
    )
    {
      math::Histogram distance_distribution( DISTANCE_RANGE.GetMin(), BIN_SIZE, DISTANCE_RANGE.GetWidth() / BIN_SIZE);

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr_1( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr_1 != coord_itr_end;
        ++coord_itr_1
      )
      {
        const linal::Vector3D &coord1( **coord_itr_1);
        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator
            coord_itr_2( coord_itr_1 + 1);
          coord_itr_2 != coord_itr_end;
          ++coord_itr_2
        )
        {
          distance_distribution.PushBack( linal::Distance( coord1, **coord_itr_2));
        }
      }

      // end
      return distance_distribution;
    }

    //! @brief create three equidistant ranges from two thresholds
    //! @param THRESHOLD lower and upper boundaries for ranges
    //! @return four threshold bounding three equivistant intervals
    storage::VectorND< 4, double>
    GeometricHashing::CreateEquidistantIntervals( const storage::VectorND< 2, double> &THRESHOLD)
    {
      const double interval_length( ( THRESHOLD.Second() - THRESHOLD.First()) * 1.0 / 3.0);
      return storage::VectorND< 4, double>
             (
               THRESHOLD.First(),
               THRESHOLD.First() + 1 * interval_length,
               THRESHOLD.First() + 2 * interval_length,
               THRESHOLD.Second()
             );
    }

    //! @brief create three ranges containing approximately equal number of distances
    //! @param DISTANCE_DISTRIBUTION histogram of distance distributions
    //! @return four threshold bounding three intervals with equal number of distances
    storage::VectorND< 4, double>
    GeometricHashing::CreateEqualOccupiedIntervals( const math::Histogram &DISTANCE_DISTRIBUTION)
    {
      // histogram without boundaries
      const linal::Vector< double> histogram( DISTANCE_DISTRIBUTION.GetHistogram());
      const size_t total_counts( size_t( histogram.Sum()));

      // desired number of counts in each interval
      const size_t desired_count( total_counts / 3);

      double threshold1( 0);
      double threshold2( 0);

      // iterate over all bins to find the correct cutoff
      for( size_t index( 0), max_index( histogram.GetSize()), sum( 0); index < max_index; ++index)
      {
        sum += histogram( index);
        if( sum > desired_count)
        {
          threshold1 = DISTANCE_DISTRIBUTION.GetBoundaries().First() + ( index + 1) * DISTANCE_DISTRIBUTION.GetBinSize();
          break;
        }
      }

      for( size_t index( 0), max_index( histogram.GetSize()), sum( 0); index < max_index; ++index)
      {
        sum += histogram( index);
        if( sum > 2 * desired_count)
        {
          threshold2 = DISTANCE_DISTRIBUTION.GetBoundaries().First() + ( index + 1) * DISTANCE_DISTRIBUTION.GetBinSize();
          break;
        }
      }

      return storage::VectorND< 4, double>
             (
               DISTANCE_DISTRIBUTION.GetBoundaries().First(),
               threshold1,
               threshold2,
               DISTANCE_DISTRIBUTION.GetBoundaries().Second()
             );
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_geometric_hash_storage_classes.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "coord/bcl_coord_geometric_hash_storage_hash_map.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GeometricHashStorageClasses::GeometricHashStorageClasses() :
      e_HashMap( AddEnum( "HashMap", util::ShPtr< GeometricHashStorageInterface>( new GeometricHashStorageHashMap())))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GeometricHashStorageClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief flag for choosing the geometris hash storage over the commandline
    const util::ShPtr< command::FlagInterface> &
    GeometricHashStorageClasses::GetGeometricHashStorageClassFlag() const
    {
      static const util::ShPtr< command::FlagInterface> s_geometric_hash_storage_class_flag
      (
        new command::FlagStatic
        (
          "hash_storage",
          "choice of storage class",
          command::Parameter
          (
            "hash_storage_class",
            "name for one of the available hash storage classes - depending on size and speed requirements",
            command::ParameterCheckEnumerate< GeometricHashStorageClasses>(),
            e_HashMap
          )
        )
      );

      // end
      return s_geometric_hash_storage_class_flag;
    }

    //! @brief construct a geometric hash storage class form the commandline argument
    util::ShPtr< GeometricHashStorageInterface>
    GeometricHashStorageClasses::ConstructFromCommandline
    (
      const std::string &MRC_NAME,
      const double MRC_RESOLUTION,
      const double MRC_VOXELSIZE,
      const storage::VectorND< 3, size_t> &MRC_EXTENSION,
      const size_t NUMBER_POINTS,
      const double FEATURE_DISTANCE,
      const double RATIO_INTENSITY_GRADIENT,
      const storage::VectorND< 4, double> &THRESHOLD,
      const double FEATURE_RADIUS,
      const PointToKeyInterface &POINT_TO_KEY,
      const size_t MIN_NUMBER_NEIGHBORS,
      const double MIN_NEIGHBOR_DISTANCE
    ) const
    {
      return util::ShPtr< GeometricHashStorageInterface>
        (
          (
            *GeometricHashStorageClass( GetGeometricHashStorageClassFlag()->GetFirstParameter()->GetValue())
          )->Construct
          (
            MRC_NAME,
            MRC_RESOLUTION,
            MRC_VOXELSIZE,
            MRC_EXTENSION,
            NUMBER_POINTS,
            FEATURE_DISTANCE,
            RATIO_INTENSITY_GRADIENT,
            THRESHOLD,
            FEATURE_RADIUS,
            POINT_TO_KEY,
            MIN_NUMBER_NEIGHBORS,
            MIN_NEIGHBOR_DISTANCE
          )
        );
    }

    //! @brief construct on access function for all GeometricHashStorageClasses
    //! @return reference to only instances of GeometricHashStorageClasses
    GeometricHashStorageClasses &GetGeometricHashStorageClasses()
    {
      return GeometricHashStorageClasses::GetEnums();
    }

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< coord::GeometricHashStorageInterface>, coord::GeometricHashStorageClasses>;

  } // namespace util
} // namespace bcl
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
#include "coord/bcl_coord_geometric_hash_storage_hash_map.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> GeometricHashStorageHashMap::s_Instance
    (
      GetObjectInstances().AddInstance( new GeometricHashStorageHashMap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default private constructor
    GeometricHashStorageHashMap::GeometricHashStorageHashMap() :
      GeometricHashStorageInterface(),
      m_MRCName(),
      m_MRCResolution(),
      m_MRCVoxelSize(),
      m_MRCExtension(),
      m_NumberPoints(),
      m_MinimalDistance(),
      m_Threshold(),
      m_Radius(),
      m_PointToKey(),
      m_MinNumberNeighbors(),
      m_MinNeighborDistance(),
      m_Final( false)
    {
    }

    //! @brief construct from all informations
    GeometricHashStorageHashMap::GeometricHashStorageHashMap
    (
      const std::string &MRC_NAME,
      const double MRC_RESOLUTION,
      const double MRC_VOXELSIZE,
      const storage::VectorND< 3, size_t> &MRC_EXTENSION,
      const size_t NUMBER_POINTS,
      const double FEATURE_DISTANCE,
      const double RATIO_INTENSITY_GRADIENT,
      const storage::VectorND< 4, double> &THRESHOLD,
      const double FEATURE_RADIUS,
      const PointToKeyInterface &POINT_TO_KEY,
      const size_t MIN_NUMBER_NEIGHBORS,
      const double MIN_NEIGHBOR_DISTANCE
    ) :
      GeometricHashStorageInterface(),
      m_MRCName( MRC_NAME),
      m_MRCResolution( MRC_RESOLUTION),
      m_MRCVoxelSize( MRC_VOXELSIZE),
      m_MRCExtension( MRC_EXTENSION),
      m_NumberPoints( NUMBER_POINTS),
      m_MinimalDistance( FEATURE_DISTANCE),
      m_RatioIntensityGradient( RATIO_INTENSITY_GRADIENT),
      m_Threshold( THRESHOLD),
      m_Radius( FEATURE_RADIUS),
      m_PointToKey( POINT_TO_KEY.Clone()),
      m_MinNumberNeighbors( MIN_NUMBER_NEIGHBORS),
      m_MinNeighborDistance( MIN_NEIGHBOR_DISTANCE),
      m_Final( false)
    {
    }

    //! virtual constructor from all information
    GeometricHashStorageInterface *GeometricHashStorageHashMap::Construct
    (
      const std::string &MRC_NAME,
      const double MRC_RESOLUTION,
      const double MRC_VOXELSIZE,
      const storage::VectorND< 3, size_t> &MRC_EXTENSION,
      const size_t NUMBER_POINTS,
      const double FEATURE_DISTANCE,
      const double RATIO_INTENSITY_GRADIENT,
      const storage::VectorND< 4, double> &THRESHOLD,
      const double FEATURE_RADIUS,
      const PointToKeyInterface &POINT_TO_KEY,
      const size_t MIN_NUMBER_NEIGHBORS,
      const double MIN_NEIGHBOR_DISTANCE
    ) const
    {
      return new GeometricHashStorageHashMap
                 (
                   MRC_NAME,
                   MRC_RESOLUTION,
                   MRC_VOXELSIZE,
                   MRC_EXTENSION,
                   NUMBER_POINTS,
                   FEATURE_DISTANCE,
                   RATIO_INTENSITY_GRADIENT,
                   THRESHOLD,
                   FEATURE_RADIUS,
                   POINT_TO_KEY,
                   MIN_NUMBER_NEIGHBORS,
                   MIN_NEIGHBOR_DISTANCE
                 );
    }

    //! @brief copy constructor
    //! @return pointer to a new GeometricHashStorageHashMap copied from this one
    GeometricHashStorageHashMap *GeometricHashStorageHashMap::Clone() const
    {
      return new GeometricHashStorageHashMap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GeometricHashStorageHashMap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief store a set of coordinates with the according Transformatiomatrix
    //! @param COORDINATES the set of coordinates that is associated with the transformation
    //! @param TRANSFORMATIONMATRIX3D the transformation for that set of coordinates
    void
    GeometricHashStorageHashMap::Store
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D
    )
    {
      BCL_Assert( !m_Final, "unable to store into a map, that is already final");
      // store the transformation matrix
      m_TransformationMatrices.PushBack( util::ShPtr< math::TransformationMatrix3D>( TRANSFORMATIONMATRIX3D.Clone()));

      // util::SiPtr to the Last Transformationmatrix3D in this Vector
      const util::SiPtr< const math::TransformationMatrix3D> sp_transformationmatrix3D( m_TransformationMatrices.LastElement());

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        // quantize the coordinate according to the point to key function
        const storage::Triplet< int, int, int> triplet_hash_key
        (
          m_PointToKey->operator()( linal::Vector3D( **coord_itr).Transform( TRANSFORMATIONMATRIX3D))
        );

        // generate hash key
        const size_t hash_key( ConvertTripletToKey( triplet_hash_key));

        // increase the count for that key and that base by 1
        ++m_Hash[ hash_key][ sp_transformationmatrix3D];
      }
    }

    //! @brief return the best counting bases
    //! @param COORDINATES the coordinates to be considered
    //! @param NUMBERBESTCOUNTS the max number of transformations to return
    //! @return the list of transformations with their hash score
    storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> >
    GeometricHashStorageHashMap::ReturnBestCounts
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
      const size_t NUMBERBESTCOUNTS
    ) const
    {
      BCL_Assert( m_Final, "unable to count within an unfinished map");

      // store the number of occurrences for each of the keys that come from the set of coordinates
      storage::Map< size_t, size_t> keys_count;

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( COORDINATES.Begin()), coord_itr_end( COORDINATES.End());
        coord_itr != coord_itr_end;
        ++coord_itr
      )
      {
        // quantize the coordinate according to the point to key function
        const storage::Triplet< int, int, int> triplet_hash_key
        (
          m_PointToKey->operator()( linal::Vector3D( **coord_itr).Transform( TRANSFORMATIONMATRIX3D))
        );

        // generate hash key
        const size_t hash_key( ConvertTripletToKey( triplet_hash_key));

        ++keys_count[ hash_key];
      }

      std::map< util::SiPtr< const math::TransformationMatrix3D>, size_t> basis_count;

      // iterate over all the keys in the target
      for
      (
        storage::Map< size_t, size_t>::const_iterator
          keys_count_itr( keys_count.Begin()), keys_count_itr_end( keys_count.End());
        keys_count_itr != keys_count_itr_end;
        ++keys_count_itr
      )
      {
        // find hash key and count the matches in the histogram
        const storage::HashMap< size_t, storage::Map< util::SiPtr< const math::TransformationMatrix3D>, size_t> >::const_iterator
          itr( m_Hash.Find( keys_count_itr->first));

        // add for each key the number of bases associated with it
        if( itr != m_Hash.End())
        {
          for
          (
            storage::Map< util::SiPtr< const math::TransformationMatrix3D>, size_t>::const_iterator
              sp_tm_itr( itr->second.Begin()), sp_tm_itr_end( itr->second.End());
            sp_tm_itr != sp_tm_itr_end;
            ++sp_tm_itr
          )
          {
            basis_count[ sp_tm_itr->first] += std::min( sp_tm_itr->second, keys_count_itr->second);
          }
        }
      }

      // these will be the best transformationmatrices sorted from highest to lowest count
      storage::List< storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t> > transformationmatrices;

      const std::multiset< std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, Sp_CountCompare>
        sorted_basis_count( basis_count.begin(), basis_count.end());
      size_t count( 0);

      // insert the best NUMBERBESTCOUNTS matrices to the best matrices beginning from the highest count
      for
      (
        std::multiset
        <
          std::pair< util::SiPtr< const math::TransformationMatrix3D>, size_t>, Sp_CountCompare
        >::const_reverse_iterator
          sorted_basis_revitr( sorted_basis_count.rbegin()),
          sorted_basis_revitr_end( sorted_basis_count.rend());
        sorted_basis_revitr != sorted_basis_revitr_end && count < NUMBERBESTCOUNTS; ++sorted_basis_revitr, ++count)
      {
        transformationmatrices.PushBack
        (
          storage::Pair< util::ShPtr< math::TransformationMatrix3D>, size_t>
          (
            util::ShPtr< math::TransformationMatrix3D>( sorted_basis_revitr->first->Clone()),
            sorted_basis_revitr->second
          )
        );
      }

      //return best transformationmatrices
      return transformationmatrices;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read Hash from io::IFStream
    std::istream &GeometricHashStorageHashMap::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TransformationMatrices, ISTREAM);
      io::Serialize::Read( m_Hash, ISTREAM);
      io::Serialize::Read( m_MRCName, ISTREAM);
      io::Serialize::Read( m_MRCResolution, ISTREAM);
      io::Serialize::Read( m_MRCVoxelSize, ISTREAM);
      io::Serialize::Read( m_MRCExtension, ISTREAM);
      io::Serialize::Read( m_NumberPoints, ISTREAM);
      io::Serialize::Read( m_MinimalDistance, ISTREAM);
      io::Serialize::Read( m_RatioIntensityGradient, ISTREAM);
      io::Serialize::Read( m_Threshold, ISTREAM);
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_PointToKey, ISTREAM);
      io::Serialize::Read( m_MinNumberNeighbors, ISTREAM);
      io::Serialize::Read( m_MinNeighborDistance, ISTREAM);

      // set final to true
      m_Final = true;

      // return
      return ISTREAM;
    }

    //! write Hash to std::ostream using the given util::Format
    std::ostream &GeometricHashStorageHashMap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // check that this is a final geometric hash storage
      BCL_Assert( m_Final, "writing of unfinished map not supported");

      // write member
      io::Serialize::Write( m_TransformationMatrices, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Hash, OSTREAM, INDENT)                   << '\n';
      io::Serialize::Write( m_MRCName, OSTREAM, INDENT)                << '\n';
      io::Serialize::Write( m_MRCResolution, OSTREAM, INDENT)          << '\n';
      io::Serialize::Write( m_MRCVoxelSize, OSTREAM, INDENT)           << '\n';
      io::Serialize::Write( m_MRCExtension, OSTREAM, INDENT)           << '\n';
      io::Serialize::Write( m_NumberPoints, OSTREAM, INDENT)           << '\n';
      io::Serialize::Write( m_MinimalDistance, OSTREAM, INDENT)        << '\n';
      io::Serialize::Write( m_RatioIntensityGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Threshold, OSTREAM, INDENT)              << '\n';
      io::Serialize::Write( m_Radius, OSTREAM, INDENT)                 << '\n';
      io::Serialize::Write( m_PointToKey, OSTREAM, INDENT)             << '\n';
      io::Serialize::Write( m_MinNumberNeighbors, OSTREAM, INDENT)     << '\n';
      io::Serialize::Write( m_MinNeighborDistance, OSTREAM, INDENT)    << '\n';

      //return
      return OSTREAM;
    }

  } // namespace coord
} // namespace bcl

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
#include "coord/bcl_coord_geometric_hash_storage_interface.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_transformation_matrix_3d.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    //! This function converts triplet of ints to a hashkey by shifting each integer by certain amount, so that in the binary representation these three integers do not overlap, given that the ints are smaller than pow(2,s_SingleKeyShift)-1
    size_t GeometricHashStorageInterface::ConvertTripletToKey( const storage::Triplet< int, int, int> &HASHKEYS)
    {
      return ( size_t( HASHKEYS.First()) << ( 2 * s_SingleKeyShift)) ^
             ( size_t( HASHKEYS.Second()) <<       s_SingleKeyShift)  ^
               size_t( HASHKEYS.Third());
    }

    //! checks if the TRANSFORMATIONMATRIX3D is similar to any stored transformation matrix according to DIFFERENCE
    bool GeometricHashStorageInterface::IsSimilarTransformation
    (
      const util::ShPtrVector< math::TransformationMatrix3D> &TRANSFORAMTIONS,
      const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D,
      const storage::VectorND< 2, double> &DIFFERENCE_ROT_TRANS
    )
    {
      // iterate over all stored transformation matrices
      for
      (
        util::ShPtrVector< math::TransformationMatrix3D>::const_iterator
          matrix_itr( TRANSFORAMTIONS.Begin()), matrix_itr_end( TRANSFORAMTIONS.End());
        matrix_itr != matrix_itr_end;
        ++matrix_itr
      )
      {
        // check if current matrix is similar within given tolerance to given TRANSFORMATIONMATRIX3D
        if
        (
          math::SimilarWithinTolerance
          (
            **matrix_itr,
            TRANSFORMATIONMATRIX3D,
            DIFFERENCE_ROT_TRANS.Second(),
            DIFFERENCE_ROT_TRANS.First()
          )
        )
        {
          return true;
        }
      }

      // no similar transformation matrix found
      return false;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_geometry_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    //! @brief boolean operator GEOMETRY_LHS == GEOMETRY_RHS
    //! @param GEOMETRY_LHS first GEOMETRY
    //! @param GEOMETRY_RHS second GEOMETRY
    //! @return whether GEOMETRY_LHS is equal to GEOMETRY_RHS
    bool operator ==( const GeometryInterface &GEOMETRY_LHS, const GeometryInterface &GEOMETRY_RHS)
    {
      return
      (
        GEOMETRY_LHS.GetOrientation() == GEOMETRY_RHS.GetOrientation() &&
        math::EqualWithinTolerance( GEOMETRY_LHS.GetMainAxis().GetLength(), GEOMETRY_RHS.GetMainAxis().GetLength()) &&
        math::EqualWithinTolerance( GEOMETRY_LHS.GetExtent( GetAxes().e_X), GEOMETRY_RHS.GetExtent( GetAxes().e_X)) &&
        math::EqualWithinTolerance( GEOMETRY_LHS.GetExtent( GetAxes().e_Y), GEOMETRY_RHS.GetExtent( GetAxes().e_Y)) &&
        GEOMETRY_LHS.GetGeometries().GetSize() == GEOMETRY_RHS.GetGeometries().GetSize()
      );
    }

    //! @brief function to check whether two GEOMETRYs are equal within tolerance
    //! @param GEOMETRY_LHS first GEOMETRY
    //! @param GEOMETRY_RHS second GEOMETRY
    //! @return whether GEOMETRY_LHS is equal to GEOMETRY_RHS within tolerance
    bool EqualWithinTolerance
    (
      const GeometryInterface &GEOMETRY_LHS,
      const GeometryInterface &GEOMETRY_RHS,
      const double &RELATIVE_TOLERANCE
    )
    {
      return
      (
        math::SimilarWithinTolerance
        (
          GEOMETRY_LHS.GetOrientation(), GEOMETRY_RHS.GetOrientation(), 0.001, 0.001
        ) &&
        math::EqualWithinTolerance
        (
          GEOMETRY_LHS.GetMainAxis().GetLength(), GEOMETRY_RHS.GetMainAxis().GetLength(), RELATIVE_TOLERANCE
        ) &&
        math::EqualWithinTolerance
        (
          GEOMETRY_LHS.GetExtent( GetAxes().e_X),
          GEOMETRY_RHS.GetExtent( GetAxes().e_X), RELATIVE_TOLERANCE
        ) &&
        math::EqualWithinTolerance
        (
          GEOMETRY_LHS.GetExtent( GetAxes().e_Y),
          GEOMETRY_RHS.GetExtent( GetAxes().e_Y), RELATIVE_TOLERANCE
        ) &&
        GEOMETRY_LHS.GetGeometries().GetSize() == GEOMETRY_RHS.GetGeometries().GetSize()
      );
    }

    //! @brief calculate the shortest connecting linesegment between two geometries
    //! @param GEOMETRY_A first geometry
    //! @param GEOMETRY_B second geometry
    //! @param MINIMAL_INTERFACE_LENGTH minimal interface length
    //! @return the shortest connecting linesegment between two geometries
    storage::Pair< LineSegment3D, bool> ShortestConnectionBetweenGeometries
    (
      const GeometryInterface &GEOMETRY_A, const GeometryInterface &GEOMETRY_B, const double MINIMAL_INTERFACE_LENGTH
    )
    {
      // make copies of the main axes
      LineSegment3D geometry_a_line( GEOMETRY_A.GetMainAxis());
      LineSegment3D geometry_b_line( GEOMETRY_B.GetMainAxis());

      // shorten them
      geometry_a_line.Shorten( MINIMAL_INTERFACE_LENGTH);
      geometry_b_line.Shorten( MINIMAL_INTERFACE_LENGTH);

      // end
      return ShortestConnectionBetweenLineSegments3D( geometry_a_line, geometry_b_line);
    }

    //! @brief returns the x, y and z extents
    //! @return the requested extents as linal::Vector3D
    linal::Vector3D GeometryInterface::GetExtents() const
    {
      return
        linal::Vector3D
        (
          GetExtent( GetAxes().e_X),
          GetExtent( GetAxes().e_Y),
          GetExtent( GetAxes().e_Z)
        );
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_moment_of_inertia.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MomentOfInertia::s_Instance
    (
      GetObjectInstances().AddInstance( new MomentOfInertia())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MomentOfInertia::MomentOfInertia()
    {
    }

    //! @brief Clone function
    //! @return pointer to new MomentOfInertia
    MomentOfInertia *MomentOfInertia::Clone() const
    {
      return new MomentOfInertia( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MomentOfInertia::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate transformation that translates into the center of weights and that sorts principal axes of inertia according to principal moments of inertia x - smallest, z - largest
    //! @param COORDINATES_WEIGHT_MATRIX matrix with coordinate and weight (4 cols) in number points rows
    //! @return transformation matrix and principal moments of inertias
    storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
    MomentOfInertia::TransformationAndMoments
    (
      const linal::MatrixConstInterface< double> &COORDINATES_WEIGHT_MATRIX
    ) const
    {
      BCL_Assert( COORDINATES_WEIGHT_MATRIX.GetNumberCols() == 4, "need 4 columns");

      const storage::Pair< math::TransformationMatrix3D, linal::Vector3D> error_return
      (
        math::TransformationMatrix3D( util::UndefinedObject()),
        linal::Vector3D( util::GetUndefined< double>())
      );

      // center of mass
      const linal::Vector3D center_of_mass( CenterOfMass( COORDINATES_WEIGHT_MATRIX));

      // error
      if( !center_of_mass.IsDefined())
      {
        return error_return;
      }

      double ixx( 0);
      double iyy( 0);
      double izz( 0);
      double ixy( 0);
      double ixz( 0);
      double iyz( 0);

      // iterate over all amino acids
      for
      (
        const double *row( COORDINATES_WEIGHT_MATRIX.Begin()), *row_end( COORDINATES_WEIGHT_MATRIX.End());
        row != row_end;
        row += 4
      )
      {
        linal::Vector3D current_coord( row);
        current_coord -= center_of_mass;

        const double current_weight( row[ 3]);

        // diagonal
        ixx += current_weight * ( math::Sqr( current_coord.Y()) + math::Sqr( current_coord.Z()));
        iyy += current_weight * ( math::Sqr( current_coord.X()) + math::Sqr( current_coord.Z()));
        izz += current_weight * ( math::Sqr( current_coord.X()) + math::Sqr( current_coord.Y()));

        // off diagonal; products of inertia
        ixy -= current_weight * current_coord.X() * current_coord.Y();
        ixz -= current_weight * current_coord.X() * current_coord.Z();
        iyz -= current_weight * current_coord.Y() * current_coord.Z();
      }

      // return the tensor
      const double elements[ 3 * 3] = { ixx, ixy, ixz, ixy, iyy, iyz, ixz, iyz, izz};
      linal::Matrix3x3< double> moment( elements);

      linal::Matrix3x3< double> eigen_vectors( 0.0);
      linal::Vector< double> eigenvalues( 3);
      moment.EigenVectorsSymmetric( eigen_vectors, eigenvalues);
      // transpose so that the eigenvectors are in the rows
      eigen_vectors.Transpose();
      // sort eigenvectors by eigenvalues
      eigen_vectors.SortRowsAndVector( eigenvalues);

      eigen_vectors.SwapRows( 0, 2);
      std::swap( eigenvalues( 0), eigenvalues( 2));

      // orthogonalize
      eigen_vectors.ReplaceRow
      (
        2,
        linal::CrossProduct( linal::Vector3D( eigen_vectors[ 0]), linal::Vector3D( eigen_vectors[ 1]))
      );

      // shift and rotate molecule
      math::TransformationMatrix3D transformation;
      transformation( math::RotationMatrix3D( eigen_vectors));
      transformation( center_of_mass);

      // end
      return
        storage::Pair< math::TransformationMatrix3D, linal::Vector3D>
        (
          math::Inverse( transformation),
          linal::Vector3D( eigenvalues)
        );
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns the ordered principal moments of inertia
    //! @param COORDINATE_WEIGHT_MATRIX 4*n matrix with three coordinates and one weight as ciolumns, and n rows
    //! @return ordered principal moments of inertia
    linal::Vector3D MomentOfInertia::operator()( const linal::MatrixConstInterface< double> &COORDINATE_WEIGHT_MATRIX) const
    {
      return TransformationAndMoments( COORDINATE_WEIGHT_MATRIX).Second();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MomentOfInertia::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MomentOfInertia::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the center of mass
    //! @param COORDINATES_WEIGHT_MATRIX matrix with coordinate and weight (4 cols) in number points rows
    //! @return the center of mass
    linal::Vector3D
    MomentOfInertia::CenterOfMass( const linal::MatrixConstInterface< double> &COORDINATES_WEIGHT_MATRIX)
    {
      BCL_Assert( COORDINATES_WEIGHT_MATRIX.GetNumberCols() == 4, "need 4 columns: 3 coordinates and 1 weight");

      // sum of weight
      double weight_sum( 0);

      // sum of weighted positions
      linal::Vector3D weighted_positions( 0.0);

      // iterate over all rows
      for
      (
        const double *row( COORDINATES_WEIGHT_MATRIX.Begin()), *row_end( COORDINATES_WEIGHT_MATRIX.End());
        row != row_end;
        row += 4
      )
      {
        // current coordinate
        linal::Vector3D current_coord( row);

        // current weight
        const double current_weight( row[ 3]);

        // weigh position
        current_coord *= current_weight;

        // add position
        weighted_positions += current_coord;

        // sum weights
        weight_sum += current_weight;
      }

      // average
      return weighted_positions / weight_sum;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_movable_eccentric.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MovableEccentric::s_Instance
    (
      GetObjectInstances().AddInstance( new MovableEccentric())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MovableEccentric::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the geometric center of the hinge
    //! @return the geometric center of the hinge
    linal::Vector3D MovableEccentric::GetCenter() const
    {
      return m_Hinge->GetCenter();
    }

    //! @brief returns the orientation of the hinge
    //! @param AXIS The requested Axis for which orientation will be returned
    //! @return the orientation of the hinge
    linal::Vector3D MovableEccentric::GetAxis( const Axis &AXIS) const
    {
      return m_Hinge->GetAxis( AXIS);
    }

    //! @brief returns the orientation and Position as TransformationMatrix3D
    //! @return the orientation and Position as TransformationMatrix3D
    const math::TransformationMatrix3D MovableEccentric::GetOrientation() const
    {
      return m_Hinge->GetOrientation();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION Vector along which the translation will occur
    void MovableEccentric::Translate( const linal::Vector3D &TRANSLATION)
    {
      m_Data->Translate( TRANSLATION);
    }

    //! @brief transform the object by a given TRANSFORMATIONMATRIX3D using m_Hinge as origin
    //! @param TRANSFORMATIONMATRIX3D TransformationMatrix3D that will be applied using m_Hinge as origin
    void MovableEccentric::Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
    {
      // apply the transformation to MovableData
      m_Data->Transform( TRANSFORMATIONMATRIX3D);
    }

    //! @brief rotate the object by a given ROTATIONMATRIX3D
    //! @brief ROTATIONMATRIX3D RotationMatrix3D that defines the rotation to be applied
    void MovableEccentric::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      // apply the rotation to MovableData
      m_Data->Rotate( ROTATIONMATRIX3D);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator =
    //! @param MOVABLE_ECCENTRIC MovableEccentric to be copied
    //! @return MovableEccentric object with members updated to MOVABLE_ECCENTRIC
    MovableEccentric &MovableEccentric::operator =( const MovableEccentric &MOVABLE_ECCENTRIC)
    {
      // update data members
      m_Data = MOVABLE_ECCENTRIC.m_Data;
      m_Hinge = MOVABLE_ECCENTRIC.m_Hinge;

      // end
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MovableEccentric::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &MovableEccentric::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }
  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_movable_interface.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief transform randomly by inner coordinates
    //! @param MAX_TRANSLATION maximal move distance
    //! @param MAX_ROTATION    maximal rotation angle
    void MovableInterface::RandomTransformation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      BCL_MessageStd( "center: " + util::Format()( GetCenter()));
      // transform this object according to the matrix
      Transform( GenerateRandomTransformationAroundCenter( MAX_TRANSLATION, MAX_ROTATION, GetCenter()));
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_combine.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoveCombine::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveCombine())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveCombine::MoveCombine() :
      m_Moves()
    {
    }

    //! @brief constructor from a list of moves
    //! @param MOVES_LIST List of moves
    MoveCombine::MoveCombine( const util::ShPtrList< MoveInterface> &MOVES_LIST) :
      m_Moves()
    {
      for( auto move_it( MOVES_LIST.Begin()); move_it != MOVES_LIST.End(); ++move_it)
      {
        m_Moves.PushBack( **move_it);
      }
    }

    //! @brief Clone function
    //! @return pointer to new MoveCombine
    MoveCombine *MoveCombine::Clone() const
    {
      return new MoveCombine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoveCombine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveCombine::GetAlias() const
    {
      static const std::string s_name( "MoveCombine");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveCombine::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Combines a list of moves and applies them consecutively.");
      serializer.AddInitializer
      (
        "moves",
        "list of moves to be applied",
        io::Serialization::GetAgent( &m_Moves)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveCombine::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // iterate over the moves
      for( auto move_itr( m_Moves.Begin()), move_itr_end( m_Moves.End()); move_itr != move_itr_end; ++move_itr)
      {
        // apply the move
        ( *move_itr)->Move( MOVEABLE_OBJECT);
      }

      // end
      return MOVEABLE_OBJECT;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_rotate_defined.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveRotateDefined::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveRotateDefined())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveRotateDefined::MoveRotateDefined() :
      m_RotationAngle(),
      m_RotationAxis(),
      m_RotateInternal()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param ROTATION_AXIS axis to rotate around
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateDefined::MoveRotateDefined
    (
      const double ROTATION_ANGLE_RAD,
      const Axis &ROTATION_AXIS,
      const bool INTERNAL
    ) :
      m_RotationAngle( ROTATION_ANGLE_RAD),
      m_RotationAxis( ROTATION_AXIS),
      m_RotateInternal( INTERNAL)
    {
    }

    //! @brief copy constructor
    //! @param MOVE_ROTATE_DEFINED MoveRotateDefined to be copied
    MoveRotateDefined::MoveRotateDefined( const MoveRotateDefined &MOVE_ROTATE_DEFINED) :
      m_RotationAngle( MOVE_ROTATE_DEFINED.m_RotationAngle),
      m_RotationAxis( MOVE_ROTATE_DEFINED.m_RotationAxis),
      m_RotateInternal( MOVE_ROTATE_DEFINED.m_RotateInternal)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveRotateDefined *MoveRotateDefined::Clone() const
    {
      return new MoveRotateDefined( *this);
    }

    //! @brief destructor
    MoveRotateDefined::~MoveRotateDefined()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveRotateDefined::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveRotateDefined::GetAlias() const
    {
      static const std::string s_name( "MoveRotateDefined");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveRotateDefined::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Rotate object around given axis and angle.");
      serializer.AddInitializer
      (
        "angle",
        "rotation angle",
        io::Serialization::GetAgent( &m_RotationAngle)
      );
      serializer.AddInitializer
      (
        "axis",
        "axis of rotation",
        io::Serialization::GetAgent( &m_RotationAxis)
      );
      serializer.AddInitializer
      (
        "internal",
        "move objects to origin before rotating them",
        io::Serialization::GetAgent( &m_RotateInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveRotateDefined::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_RotateInternal)
      {
        // transform with generated transformationmatrix
        MOVEABLE_OBJECT.Transform( GenerateTransformationInternal( MOVEABLE_OBJECT));
      }
      else
      {
        // rotate with generate rotationmatrix
        MOVEABLE_OBJECT.Rotate( GenerateRotation());
      }

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix3D for internal coordinates (center of MOVEABLE_OBJECT as origin)
    //! @param MOVEABLE_OBJECT
    //! @return a translation vector
    const math::TransformationMatrix3D MoveRotateDefined::GenerateTransformationInternal( const MovableInterface &MOVEABLE_OBJECT) const
    {
      // apply transformation to origin
      math::TransformationMatrix3D new_matrix( math::Inverse( MOVEABLE_OBJECT.GetOrientation()));

      // generate random rotation
      new_matrix( GenerateRotation());

      // move back to original position in space
      new_matrix( MOVEABLE_OBJECT.GetOrientation());

      // end
      return new_matrix;
    }

    //! @brief generate a RotationMatrix3D
    //! @return a translation vector
    const math::RotationMatrix3D MoveRotateDefined::GenerateRotation() const
    {
      return math::RotationMatrix3D( m_RotationAxis, m_RotationAngle);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    //! @param MOVE_ROTATE_DEFINED MoveRotateDefined to be copied
    MoveRotateDefined &MoveRotateDefined::operator =( const MoveRotateDefined &MOVE_ROTATE_DEFINED)
    {
      // update members
      m_RotationAngle  = MOVE_ROTATE_DEFINED.m_RotationAngle;
      m_RotationAxis   = MOVE_ROTATE_DEFINED.m_RotationAxis;
      m_RotateInternal = MOVE_ROTATE_DEFINED.m_RotateInternal;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns an instance of MoveRotateDefined set to 180 degrees which forms a flip move around given AXIS
    //! @param AXIS Axis around which the flip is going to be applied
    //! @return an instance of MoveRotateDefined set to 180 degrees which forms a flip move around given AXIS
    MoveRotateDefined MoveRotateDefined::GetFlipMove( const Axis &AXIS)
    {
      // return MoveRotateDefined that rotates 180 degrees (flips) around provided AXIS
      return MoveRotateDefined( math::g_Pi, AXIS);
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_rotate_random.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveRotateRandom::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveRotateRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveRotateRandom::MoveRotateRandom()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_RotateInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_RotateInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const double MIN_ROTATION_ANGLE_RAD,
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_RotateInternal( INTERNAL)
    {

    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLES_RAD minimal angle in radian for rotation for each axis
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param INTERNAL if true, object will be centered in origin before applying rotation and move it back
    MoveRotateRandom::MoveRotateRandom
    (
      const linal::Vector3D &MIN_ROTATION_ANGLES_RAD,
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const bool INTERNAL
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLES_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_RotateInternal( INTERNAL)
    {

    }

    //! @brief copy constructor
    //! @param MOVE_ROTATE_RANDOM MoveRotateRandom to be copied
    MoveRotateRandom::MoveRotateRandom( const MoveRotateRandom &MOVE_ROTATE_RANDOM) :
      m_MinRotationAngles( MOVE_ROTATE_RANDOM.m_MinRotationAngles),
      m_MaxRotationAngles( MOVE_ROTATE_RANDOM.m_MaxRotationAngles),
      m_RotateInternal( MOVE_ROTATE_RANDOM.m_RotateInternal)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveRotateRandom *MoveRotateRandom::Clone() const
    {
      return new MoveRotateRandom( *this);
    }

    //! @brief destructor
    MoveRotateRandom::~MoveRotateRandom()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveRotateRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveRotateRandom::GetAlias() const
    {
      static const std::string s_name( "MoveRotateRandom");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveRotateRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Random rotation of an object.");
      serializer.AddInitializer
      (
        "min angles",
        "minimum angles for rotation",
        io::Serialization::GetAgent( &m_MinRotationAngles),
        "(0,0,0)"
      );
      serializer.AddInitializer
      (
        "max angles",
        "maximum angles for rotation",
        io::Serialization::GetAgent( &m_MaxRotationAngles)
      );
      serializer.AddInitializer
      (
        "internal",
        "move object to origin before rotation",
        io::Serialization::GetAgent( &m_RotateInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveRotateRandom::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_RotateInternal)
      {
        // transform with generated transformationmatrix
        MOVEABLE_OBJECT.Transform( GenerateRandomTransformationInternal( MOVEABLE_OBJECT));
      }
      else
      {
        // rotate with generate rotationmatrix
        MOVEABLE_OBJECT.Rotate( GenerateRandomRotation());
      }

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix3D for internal coordinates (center of MOVEABLE_OBJECT as origin)
    //! @param MOVEABLE_OBJECT
    //! @return a translation vector
    const math::TransformationMatrix3D
    MoveRotateRandom::GenerateRandomTransformationInternal
    (
      const MovableInterface &MOVEABLE_OBJECT
    ) const
    {
      // apply transformation to origin
      math::TransformationMatrix3D new_matrix( math::Inverse( MOVEABLE_OBJECT.GetOrientation()));

      // generate random rotation
      new_matrix( GenerateRandomRotation());

      // move back to original position in space
      new_matrix( MOVEABLE_OBJECT.GetOrientation());

      // end
      return new_matrix;
    }

    //! @brief generate a RotationMatrix3D
    //! @return a translation vector
    const math::RotationMatrix3D MoveRotateRandom::GenerateRandomRotation() const
    {
      return GenerateRandomRotation( m_MaxRotationAngles, m_MinRotationAngles);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate a RotationMatrix3D given the max rotation angles for each axis
    //! @return a translation vector
    const math::RotationMatrix3D MoveRotateRandom::GenerateRandomRotation
    (
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const linal::Vector3D &MIN_ROTATION_ANGLES_RAD
    )
    {
      // calculate the difference
      const linal::Vector3D diff_angles( MAX_ROTATION_ANGLES_RAD - MIN_ROTATION_ANGLES_RAD);

      double rotation_around_x
      (
        random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0)) * diff_angles.X()
      );
      double rotation_around_y
      (
        random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0)) * diff_angles.Y()
      );
      double rotation_around_z
      (
        random::GetGlobalRandom().Random< double>( double( -1.0), double( 1.0)) * diff_angles.Z()
      );

      // add the minimum rotations
      rotation_around_x += rotation_around_x < 0.0 ? -MIN_ROTATION_ANGLES_RAD.X() : MIN_ROTATION_ANGLES_RAD.X();
      rotation_around_y += rotation_around_y < 0.0 ? -MIN_ROTATION_ANGLES_RAD.Y() : MIN_ROTATION_ANGLES_RAD.Y();
      rotation_around_z += rotation_around_z < 0.0 ? -MIN_ROTATION_ANGLES_RAD.Z() : MIN_ROTATION_ANGLES_RAD.Z();

      // construct the rotatio matrix from these three rotations
      math::RotationMatrix3D rotation( GetAxes().e_X, rotation_around_x);
      rotation *= math::RotationMatrix3D( GetAxes().e_Y, rotation_around_y);
      rotation *= math::RotationMatrix3D( GetAxes().e_Z, rotation_around_z);

      // end
      return rotation;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_rotate_random_external_reference.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

    //! @brief MethodType as string
    //! @param METHOD_TYPE the MethodType
    //! @return the string for the MethodType
    const std::string &MoveRotateRandomExternalReference::GetMethodDescriptor( const MethodType &METHOD_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "INTERNAL", "INTERNAL_ROTATE", "INTERNAL_TRANSLATE", "EXTERNAL", GetStaticClassName< MethodType>()
      };

      return s_descriptors[ METHOD_TYPE];
    }

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveRotateRandomExternalReference::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveRotateRandomExternalReference())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference() :
      m_MinRotationAngles(),
      m_MaxRotationAngles(),
      m_ReferenceOrientation(),
      m_Method( e_Internal)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const double MAX_ROTATION_ANGLE_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const double MIN_ROTATION_ANGLE_RAD,
      const double MAX_ROTATION_ANGLE_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_ROTATION_ANGLES_RAD minimal angle in radian for rotation for each axis
    //! @param MAX_ROTATION_ANGLES_RAD maximal angle in radian for rotation for each axis
    //! @param REFERENCE_ORIENTATION external reference orientation
    //! @param ROTATION_METHOD rotation method to be used
    MoveRotateRandomExternalReference::MoveRotateRandomExternalReference
    (
      const linal::Vector3D &MIN_ROTATION_ANGLES_RAD,
      const linal::Vector3D &MAX_ROTATION_ANGLES_RAD,
      const math::TransformationMatrix3D &REFERENCE_ORIENTATION,
      const MethodType ROTATION_METHOD
    ) :
      m_MinRotationAngles( MIN_ROTATION_ANGLES_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLES_RAD),
      m_ReferenceOrientation( REFERENCE_ORIENTATION),
      m_Method( ROTATION_METHOD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoveRotateRandomExternalReference
    MoveRotateRandomExternalReference *MoveRotateRandomExternalReference::Clone() const
    {
      return new MoveRotateRandomExternalReference( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoveRotateRandomExternalReference::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveRotateRandomExternalReference::GetAlias() const
    {
      static const std::string s_name( "MoveRotateRandomExternalReference");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveRotateRandomExternalReference::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies a random rotation relative to an external coordinate system.");
      serializer.AddInitializer
      (
        "min rotation",
        "minimum rotation angle",
        io::Serialization::GetAgent( &m_MinRotationAngles)
      );
      serializer.AddInitializer
      (
        "max rotation",
        "maximum rotation angles",
        io::Serialization::GetAgent( &m_MaxRotationAngles)
      );
      // serializer.AddInitializer
      // (
      //   "reference orientation",
      //   "external reference orientation",
      //   io::Serialization::GetAgent( &m_ReferenceOrientation)
      // );
      serializer.AddInitializer
      (
        "method",
        "rotation method",
        io::Serialization::GetAgent( &m_Method)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveRotateRandomExternalReference::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // construct a transformation matrix to be applied to the object
      math::TransformationMatrix3D transform;

      // apply the move using the method defined by "m_Method"
      switch( m_Method)
      {
        case e_Internal:
          // move to the external reference
          transform( math::Inverse( MOVEABLE_OBJECT.GetOrientation()));
          transform( m_ReferenceOrientation);

          // generate random rotation
          transform( GenerateRandomRotation());

          // move back to the original position
          transform( MOVEABLE_OBJECT.GetOrientation());

          break;

        case e_InternalRotate:
          // rotate to match the external reference
          transform( math::Inverse( MOVEABLE_OBJECT.GetOrientation().GetRotation()));
          transform( m_ReferenceOrientation.GetRotation());

          // generate random rotation
          transform( GenerateRandomRotation());

          // rotate back to the original position
          transform( MOVEABLE_OBJECT.GetOrientation().GetRotation());

          break;

        case e_InternalTranslate:
          // translate to the external reference
          transform( -MOVEABLE_OBJECT.GetCenter());
          transform( m_ReferenceOrientation.GetOrigin());

          // generate random rotation
          transform( GenerateRandomRotation());

          // move back to the original position
          transform( MOVEABLE_OBJECT.GetCenter());

          break;

        case e_External:
          // generate random rotation
          transform( GenerateRandomRotation());

          break;
        default: break;
      }

      // apply the transformation
      MOVEABLE_OBJECT.Transform( transform);

      // end
      return MOVEABLE_OBJECT;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate a random rotation
    //! @return transformation matrix for the rotation
    math::TransformationMatrix3D MoveRotateRandomExternalReference::GenerateRandomRotation() const
    {
      // initialize transformation matrix by moving reference to origin
      math::TransformationMatrix3D transform( math::Inverse( m_ReferenceOrientation));

      // apply a random rotation
      transform( MoveRotateRandom::GenerateRandomRotation( m_MaxRotationAngles, m_MinRotationAngles));

      // move back to the reference position
      transform( m_ReferenceOrientation);

      // end
      return transform;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_transform_random.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "coord/bcl_coord_move_rotate_random.h"
#include "coord/bcl_coord_move_translate_random.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveTransformRandom::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTransformRandom())
    );

    //! @brief default constructor
    MoveTransformRandom::MoveTransformRandom()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATION        maximal distance for translation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const double MAX_TRANSLATION,
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATION        maximal distance for translation for each direction
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radians for rotation for each direction
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const linal::Vector3D &MAX_TRANSLATION,
      const linal::Vector3D &MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles(),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATION        minimal distance for translation
    //! @param MAX_TRANSLATION        maximal distance for translation
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const double MIN_TRANSLATION,
      const double MAX_TRANSLATION,
      const double MIN_ROTATION_ANGLE_RAD,
      const double MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATION        minimal distance for translation for each direction
    //! @param MAX_TRANSLATION        maximal distance for translation for each direction
    //! @param MIN_ROTATION_ANGLE_RAD minimal angle in radian for rotation for each direction
    //! @param MAX_ROTATION_ANGLE_RAD maximal angle in radian for rotation for each direction
    //! @param INTERNAL is set to true, rotation and translation is performed relative to the center
    MoveTransformRandom::MoveTransformRandom
    (
      const linal::Vector3D &MIN_TRANSLATION,
      const linal::Vector3D &MAX_TRANSLATION,
      const linal::Vector3D &MIN_ROTATION_ANGLE_RAD,
      const linal::Vector3D &MAX_ROTATION_ANGLE_RAD,
      const bool INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_MinRotationAngles( MIN_ROTATION_ANGLE_RAD),
      m_MaxRotationAngles( MAX_ROTATION_ANGLE_RAD),
      m_TransformInternal( INTERNAL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveTransformRandom *MoveTransformRandom::Clone() const
    {
      return new MoveTransformRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveTransformRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTransformRandom::GetAlias() const
    {
      static const std::string s_name( "MoveTransformRandom");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTransformRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies a random transformation.");
      serializer.AddInitializer
      (
        "min translation",
        "minimum translation",
        io::Serialization::GetAgent( &m_MinTranslation)
      );
      serializer.AddInitializer
      (
        "max translation",
        "maximum translation",
        io::Serialization::GetAgent( &m_MaxTranslation)
      );
      serializer.AddInitializer
      (
        "min rotation",
        "minimum rotation",
        io::Serialization::GetAgent( &m_MinRotationAngles)
      );
      serializer.AddInitializer
      (
        "max rotation",
        "maximum rotation",
        io::Serialization::GetAgent( &m_MaxRotationAngles)
      );
      serializer.AddInitializer
      (
        "internal transformation",
        "move objects to origin before applying transformation",
        io::Serialization::GetAgent( &m_TransformInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTransformRandom::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // transform with generated transformationmatrix
      MOVEABLE_OBJECT.Transform( GenerateRandomTransformation( MOVEABLE_OBJECT));

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix
    //! @param MOVEABLE_OBJECT
    //! @return Transformationmatrix
    const math::TransformationMatrix3D MoveTransformRandom::GenerateRandomTransformation( const MovableInterface &MOVEABLE_OBJECT) const
    {
      // new transformation matrix
      math::TransformationMatrix3D transformation;
      math::TransformationMatrix3D orientation;

      // random rotation and translation
      const math::RotationMatrix3D rand_rot( MoveRotateRandom::GenerateRandomRotation( m_MaxRotationAngles, m_MinRotationAngles));
      linal::Vector3D rand_trans( MoveTranslateRandom::GenerateRandomTranslation( m_MaxTranslation, m_MinTranslation));

      // move Movable to origin
      if( m_TransformInternal)
      {
        // get the start orientation
        orientation = MOVEABLE_OBJECT.GetOrientation();
        // move to origin
        transformation( math::Inverse( orientation));
      }

      // apply random rotation
      transformation( rand_rot);

      // move Movable back to original position and orientation
      if( m_TransformInternal)
      {
        // move back to original position
        transformation( orientation);
        // correct the translation
        rand_trans.Rotate( rand_rot);
      }

      // apply random transformation
      transformation( rand_trans);

      // end
      return transformation;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_translate_defined.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveTranslateDefined::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTranslateDefined())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveTranslateDefined::MoveTranslateDefined() :
      m_Translation(),
      m_TranslateInternal( true)
    {
    }

    //! @brief construct from  translation and internal coordinates or not
    //! @param TRANSLATION a translation vector
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateDefined::MoveTranslateDefined( const linal::Vector3D &TRANSLATION, const bool TRANSLATE_INTERNAL) :
      m_Translation( TRANSLATION),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveTranslateDefined *MoveTranslateDefined::Clone() const
    {
      return new MoveTranslateDefined( *this);
    }

    //! @brief virtual destructor
    MoveTranslateDefined::~MoveTranslateDefined()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveTranslateDefined::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTranslateDefined::GetAlias() const
    {
      static const std::string s_name( "MoveTranslateDefined");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTranslateDefined::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "translate object");
      serializer.AddInitializer
      (
        "translation",
        "translation vector",
        io::Serialization::GetAgent( &m_Translation)
      );
      serializer.AddInitializer
      (
        "internal",
        "move to the origin before applying the translation",
        io::Serialization::GetAgent( &m_TranslateInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT reference on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTranslateDefined::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_TranslateInternal)
      {
        // transform with generated translation relative to internal orientation
        MOVEABLE_OBJECT.Translate( linal::Vector3D( m_Translation).Rotate( MOVEABLE_OBJECT.GetOrientation().GetRotation()));
      }
      else
      {
        // transform with generated translation relative to normal coordinatesystem
        MOVEABLE_OBJECT.Translate( m_Translation);
      }

      // end
      return MOVEABLE_OBJECT;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_move_translate_external_axis.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_range.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MoveTranslateExternalAxis::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTranslateExternalAxis())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveTranslateExternalAxis::MoveTranslateExternalAxis() :
      m_MinTranslation( util::GetUndefinedDouble()),
      m_MaxTranslation( util::GetUndefinedDouble()),
      m_ExternalAxis( GetAxes().e_Undefined)
    {
    }

    //! @brief construct from translation range and axis
    //! @param MIN_TRANSLATION minimal distance for translation
    //! @param MAX_TRANSLATION maximal distance for translation
    //! @param EXTERNAL_AXIS axis to translate towards or away from
    MoveTranslateExternalAxis::MoveTranslateExternalAxis
    (
      const double MIN_TRANSLATION,
      const double MAX_TRANSLATION,
      const Axis &EXTERNAL_AXIS
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_ExternalAxis( EXTERNAL_AXIS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoveTranslateExternalAxis
    MoveTranslateExternalAxis *MoveTranslateExternalAxis::Clone() const
    {
      return new MoveTranslateExternalAxis( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &MoveTranslateExternalAxis::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTranslateExternalAxis::GetAlias() const
    {
      static const std::string s_name( "MoveTranslateExternalAxis");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTranslateExternalAxis::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Translate along an external axis.");
      serializer.AddInitializer
      (
        "min translation",
        "minimum translation distance",
        io::Serialization::GetAgent( &m_MinTranslation)
      );
      serializer.AddInitializer
      (
        "max translation",
        "maximum translation distance",
        io::Serialization::GetAgent( &m_MaxTranslation)
      );
      serializer.AddInitializer
      (
        "external axis",
        "translation axis",
        io::Serialization::GetAgent( &m_ExternalAxis)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTranslateExternalAxis::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      // get a random translation distance
      double distance
      (
        random::GetGlobalRandom().Double( math::Range< double>( m_MinTranslation, m_MaxTranslation))
      );

      // randomly decide to flip the sign
      if( random::GetGlobalRandom().Boolean())
      {
        distance *= -1.0;
      }

      // initialize translation vector
      linal::Vector3D translate;

      // if the axis is the x-axis
      if( m_ExternalAxis == GetAxes().e_X)
      {
        // get the vector from the x-axis to the object
        translate = linal::Vector3D( 0.0, MOVEABLE_OBJECT.GetCenter().Y(), MOVEABLE_OBJECT.GetCenter().Z());
      }
      // if the axis is the y-axis
      else if( m_ExternalAxis == GetAxes().e_Y)
      {
        // get the vector from the y-axis to the object
        translate = linal::Vector3D( MOVEABLE_OBJECT.GetCenter().X(), 0.0, MOVEABLE_OBJECT.GetCenter().Z());
      }
      // if the axis is the z-axis
      else if( m_ExternalAxis == GetAxes().e_Z)
      {
        // get the vector from the z-axis to the object
        translate = linal::Vector3D( MOVEABLE_OBJECT.GetCenter().X(), MOVEABLE_OBJECT.GetCenter().Y(), 0.0);
      }

      // normalize the vector
      translate.Normalize();

      // translate the object by the chosen distance
      MOVEABLE_OBJECT.Translate( distance * translate);

      // end
      return MOVEABLE_OBJECT;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace coord

} // namespace bcl
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
#include "coord/bcl_coord_move_translate_random.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_movable_interface.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MoveTranslateRandom::s_Instance
    (
      util::Enumerated< MoveInterface>::AddInstance( new MoveTranslateRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoveTranslateRandom::MoveTranslateRandom()
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATION        maximal distance for translation
    MoveTranslateRandom::MoveTranslateRandom
    (
      const double MAX_TRANSLATION,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATION),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MAX_TRANSLATIONS maximal distance for translation in each direction
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateRandom::MoveTranslateRandom
    (
      const linal::Vector3D &MAX_TRANSLATIONS,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation(),
      m_MaxTranslation( MAX_TRANSLATIONS),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATION minimal distance for translation in each direction (x, y, z)
    //! @param MAX_TRANSLATION maximal distance for translation in each direction (x, y, z)
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateRandom::MoveTranslateRandom
    (
      const double MIN_TRANSLATION,
      const double MAX_TRANSLATION,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATION),
      m_MaxTranslation( MAX_TRANSLATION),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {
    }

    //! @brief construct from undirected maximal translation and rotation
    //! @param MIN_TRANSLATIONS minimal distances for translation in each direction (x, y, z)
    //! @param MAX_TRANSLATIONS maximal distances for translation in each direction (x, y, z)
    //! @param TRANSLATE_INTERNAL maximal distance for translation in each direction
    MoveTranslateRandom::MoveTranslateRandom
    (
      const linal::Vector3D &MIN_TRANSLATIONS,
      const linal::Vector3D &MAX_TRANSLATIONS,
      const bool TRANSLATE_INTERNAL
    ) :
      m_MinTranslation( MIN_TRANSLATIONS),
      m_MaxTranslation( MAX_TRANSLATIONS),
      m_TranslateInternal( TRANSLATE_INTERNAL)
    {

    }

    //! @brief Clone function
    //! @return pointer to new Move
    MoveTranslateRandom *MoveTranslateRandom::Clone() const
    {
      return new MoveTranslateRandom( *this);
    }

    //! @brief destructor
    MoveTranslateRandom::~MoveTranslateRandom()
    {}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoveTranslateRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &MoveTranslateRandom::GetAlias() const
    {
      static const std::string s_name( "MoveTranslateRandom");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoveTranslateRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Applies a random translation.");
      serializer.AddInitializer
      (
        "min translation",
        "minimum translation",
        io::Serialization::GetAgent( &m_MinTranslation)
      );
      serializer.AddInitializer
      (
        "max translation",
        "maximum translation",
        io::Serialization::GetAgent( &m_MaxTranslation)
      );
      serializer.AddInitializer
      (
        "internal",
        "move to origin before applying translation",
        io::Serialization::GetAgent( &m_TranslateInternal),
        "true"
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the MovableInterface derived object
    //! @param MOVEABLE_OBJECT refernce on MovableInterface derived object
    //! @return reference to movable object
    MovableInterface &MoveTranslateRandom::Move( MovableInterface &MOVEABLE_OBJECT) const
    {
      if( m_TranslateInternal)
      {
        // transform with generated translation relative to internal orientation
        MOVEABLE_OBJECT.Translate( GenerateRandomTranslation().Rotate( MOVEABLE_OBJECT.GetOrientation().GetRotation()));
      }
      else
      {
        // transform with generated translation relative to normal coordinatesystem
        MOVEABLE_OBJECT.Translate( GenerateRandomTranslation());
      }

      // end
      return MOVEABLE_OBJECT;
    }

    //! @brief generate a TransformationMatrix
    //! @return a translation vector
    linal::Vector3D MoveTranslateRandom::GenerateRandomTranslation() const
    {
      return GenerateRandomTranslation( m_MaxTranslation, m_MinTranslation);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate a TransformationMatrix
    //! @param MIN_TRANSLATIONS
    //! @param MAX_TRANSLATIONS
    //! @return a translation vector
    linal::Vector3D MoveTranslateRandom::GenerateRandomTranslation
    (
      const linal::Vector3D &MAX_TRANSLATIONS,
      const linal::Vector3D &MIN_TRANSLATIONS
    )
    {
      // initialize a translation vector and set to random translation
      linal::Vector3D translation;
      translation.SetRandomTranslation( linal::Vector3D( MAX_TRANSLATIONS - MIN_TRANSLATIONS));

      // now add the min translations, if random translation is negative then subtract the translation otherwise add
      translation.X() += translation.X() < 0.0 ? -MIN_TRANSLATIONS.X() : MIN_TRANSLATIONS.X();
      translation.Y() += translation.Y() < 0.0 ? -MIN_TRANSLATIONS.Y() : MIN_TRANSLATIONS.Y();
      translation.Z() += translation.Z() < 0.0 ? -MIN_TRANSLATIONS.Z() : MIN_TRANSLATIONS.Z();

      // end
      return translation;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_orientation_interface.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate a random transformation around a given center
    //! @param MAX_TRANSLATION maximal translation
    //! @param MAX_ROTATION    maximal rotation angle (in rad)
    //! @param CENTER          origin around which the transformation will occur
    //! @return generated TransformationMatrix
    math::TransformationMatrix3D
    OrientationInterface::GenerateRandomTransformationAroundCenter
    (
      const double MAX_TRANSLATION,
      const double MAX_ROTATION,
      const linal::Vector3D &CENTER
    )
    {
      //move center of mass in origin
      math::TransformationMatrix3D new_matrix( -CENTER);

      //apply random rotation
      new_matrix( math::RotationMatrix3D().SetRand( MAX_ROTATION));

      //apply random translation
      new_matrix( linal::Vector3D().SetRandomTranslation( MAX_TRANSLATION));

      //move center of mass back to original position
      new_matrix( CENTER);

      // end
      return new_matrix;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_point_cloud.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointCloud::s_Instance
    (
      GetObjectInstances().AddInstance( new PointCloud())
    );

    //! Translate the entire pointcloud to the center of mass
    void PointCloud::CenterPointCloud()
    {

      // determine the center
      const linal::Vector3D center( GetCenter());

      // construct SiPtVector on all corrdinate
      util::SiPtrVector< linal::Vector3D> coords( util::ConvertToSiPtrVector( m_Data));

      // transform all coordinates - move them to the center
      TransformCoordinates( coords, math::TransformationMatrix3D( -center));
    }

    storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
    PointCloud::DetermineNeighbors( const double DISTANCE) const
    {
      const double square_distance( math::Sqr( DISTANCE));
      storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> > neighbors;

      for( const_iterator itr_a( Begin()), itr_end( End()); itr_a != itr_end; ++itr_a)
      {
        const util::SiPtr< const linal::Vector3D> point_a( *itr_a);

        for( const_iterator itr_b( itr_a + 1); itr_b != itr_end; ++itr_b)
        {
          const linal::Vector3D delta_xyz( linal::AbsoluteDifference( *itr_a, *itr_b));
          if( delta_xyz.X() < DISTANCE && delta_xyz.Y() < DISTANCE && delta_xyz.Z() < DISTANCE)
          {
            double current_square_distance = delta_xyz.SquareNorm();
            if( current_square_distance < square_distance)
            {
              const util::SiPtr< const linal::Vector3D> point_b( *itr_b);

              // point b is neighbor of point a
              neighbors[ point_a].PushBack( point_b);
              // point a is neighbor of point b
              neighbors[ point_b].PushBack( point_a);
            }
          }
        }
      }

      return neighbors;
    }

    size_t PointCloud::RemoveSingles( const size_t MIN_NEIGHBORS, const double DISTANCE)
    {
      const storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        neighbors( DetermineNeighbors( DISTANCE));

      PointCloud cleaned_cloud;
      for( const_iterator itr( Begin()), itr_end( End()); itr != itr_end; ++itr)
      {
        storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >::const_iterator
          find_itr( neighbors.Find( util::ToSiPtr( *itr)));
        if( find_itr == neighbors.End() || find_itr->second.GetSize() >= MIN_NEIGHBORS)
        {
          cleaned_cloud.PushBack( *itr);
        }
      }

      const size_t number_removed_points( GetSize() - cleaned_cloud.GetSize());
      *this = cleaned_cloud;

      return number_removed_points;
    }

    //! write this PointCloud to a given FILE
    //! @param FILENAME
    //! pdb file where pointcloud should be written to
    void PointCloud::WriteToPDB( const std::string &FILENAME) const
    {
      //instatiate write stream and check that it was possible to open
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      size_t atom_serial( 1);
      size_t res_serial( 1);
      pdb::Handler pdb;

      //iterate over all points in this PointCloud
      for
      (
        storage::Vector< linal::Vector3D>::const_iterator point_itr( Begin()), point_itr_end( End());
          point_itr != point_itr_end;
        ++point_itr, ++res_serial, ++atom_serial
      )
      {
        //instantiate new atom line
        util::ShPtr< pdb::Line> atom_line( new pdb::Line( pdb::GetLineTypes().ATOM));

        //write pseuda Glycine with CA containing the current point to pdb line
        atom_line->Put( pdb::GetEntryTypes().ATOMSerial, atom_serial);
        atom_line->Put( pdb::GetEntryTypes().ATOMName, "CA");
        atom_line->Put( pdb::GetEntryTypes().ATOMResidueName, "GLY");
        atom_line->Put( pdb::GetEntryTypes().ATOMChainID, 'A');
        atom_line->Put( pdb::GetEntryTypes().ATOMResidueSequenceID, res_serial);
        atom_line->PutCoordinates( *point_itr);

        // write lines in pdb
        pdb.PushBack( atom_line);
      }

      pdb.WriteLines( write);
      io::File::CloseClearFStream( write);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PointCloud::Read( std::istream &ISTREAM)
    {
      // read base classes
      io::Serialize::Read( m_Data, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &PointCloud::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write base classes
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace coord
} // namespace bcl
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

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "coord/bcl_coord_axes.h"
#include "coord/bcl_coord_point_to_key_classes.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeyCartesian
    //! @brief operator will quantize point to three ints according to a Cartesian frame
    //!
    //! @see @link example_coord_point_to_key_cartesian.cpp @endlink
    //! @author woetzen
    //! @date Aug 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API PointToKeyCartesian :
      public PointToKeyInterface
    {
    private:
    //////////
    // data //
    //////////

      static const std::string &GetCoordinateSystemDescriptor();

      //! @brief distance resolution
      double m_DistanceResolution;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! instance of the PointToKey enum
      static PointToKey e_Cartesian;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from Resolution
      PointToKeyCartesian( const double DISTANCE_RESOLUTION);

      //! virtual copy constructor
      PointToKeyCartesian *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the CoordinateSystem as string
      //! @return string describing the coordinate system
      const std::string &GetCoordinateSystem() const;

      //! @brief get angular resolution
      //! @return angular resolution
      double GetAngularResolution() const;

      //! @brief set angular resolution
      //! @param ANGULAR_RESOLUTION for this convert function
      void SetAngularResolution( const double ANGULAR_RESOLUTION);

      //! @brief get distance resolution
      //! @return distance resolution
      double GetDistanceResolution() const;

      //! @brief set distance resolution
      //! @param DISTANCE_RESOLUTION for this convert function
      void SetDistanceResolution( const double DISTANCE_RESOLUTION);

      //! operator will quantize point to three ints within a coordinate frame
      storage::Triplet< int, int, int> operator()( const linal::Vector3D &POINT) const;

      //! operator converting Triplet to a point
      linal::Vector3D operator()( const storage::Triplet< int, int, int> &TRIPLET) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! write to OSTREAM
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read from ISTREAM
      std::istream &Read( std::istream &ISTREAM);

    }; // class PointToKeyCartesian

  //////////
  // data //
  //////////

    const std::string &PointToKeyCartesian::GetCoordinateSystemDescriptor()
    {
      static const std::string s_coordinate_system( "Cartesian");
      return s_coordinate_system;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointToKeyCartesian::s_Instance
    (
      GetObjectInstances().AddInstance( new PointToKeyCartesian( PointToKeyInterface::GetDefaultDistanceResolution()))
    );

    //! instance of the PointToKey enum
    PointToKey PointToKeyCartesian::e_Cartesian
    (
      GetPointToKeyClasses().AddEnum( PointToKeyCartesian::GetCoordinateSystemDescriptor(), util::ShPtr< PointToKeyInterface>( new PointToKeyCartesian( PointToKeyInterface::GetDefaultDistanceResolution())))
    );

    //! construct from Resolution
    PointToKeyCartesian::PointToKeyCartesian( const double DISTANCE_RESOLUTION) :
      m_DistanceResolution( DISTANCE_RESOLUTION)
    {
    }

    //! virtual copy constructor
    PointToKeyCartesian *PointToKeyCartesian::Clone() const
    {
      return new PointToKeyCartesian( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeyCartesian::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the CoordinateSystem as string
    //! @return string describing the coordinate system
    const std::string &PointToKeyCartesian::GetCoordinateSystem() const
    {
      // return
      return GetCoordinateSystemDescriptor();
    }

    //! @brief get angular resolution
    //! @return angular resolution
    double PointToKeyCartesian::GetAngularResolution() const
    {
      return util::GetUndefinedSize_t();
    }

    //! @brief set angular resolution
    //! @param ANGULAR_RESOLUTION for this convert function
    void PointToKeyCartesian::SetAngularResolution( const double ANGULAR_RESOLUTION)
    {
      return;
    }

    //! @brief get distance resolution
    //! @return distance resolution
    double PointToKeyCartesian::GetDistanceResolution() const
    {
      return m_DistanceResolution;
    }

    //! @brief set distance resolution
    //! @param DISTANCE_RESOLUTION for this convert function
    void PointToKeyCartesian::SetDistanceResolution( const double DISTANCE_RESOLUTION)
    {
      m_DistanceResolution = DISTANCE_RESOLUTION;
    }

    //! operator will quantize point to three ints within a coordinate frame
    storage::Triplet< int, int, int> PointToKeyCartesian::operator()( const linal::Vector3D &POINT) const
    {
      // convert x y z to size_t and convert to hash key
      return storage::Triplet< int, int, int>
             (
               int( POINT.X() / m_DistanceResolution),
               int( POINT.Y() / m_DistanceResolution),
               int( POINT.Z() / m_DistanceResolution)
             );
    }

    //! operator converting Triplet to a point
    linal::Vector3D PointToKeyCartesian::operator()( const storage::Triplet< int, int, int> &TRIPLET) const
    {
      return linal::Vector3D
             (
               double( TRIPLET.First())  * m_DistanceResolution,
               double( TRIPLET.Second()) * m_DistanceResolution,
               double( TRIPLET.Third())  * m_DistanceResolution
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write to OSTREAM
    std::ostream &PointToKeyCartesian::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_DistanceResolution, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read from ISTREAM
    std::istream &PointToKeyCartesian::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_DistanceResolution, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_point_to_key_classes.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_point_to_key_spherical_radius.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PointToKeyClasses::PointToKeyClasses() :
      e_SphericalRadius( AddEnum( PointToKeySphericalRadius::GetCoordinateSystemDescriptor(), util::ShPtr< PointToKeyInterface>( new PointToKeySphericalRadius( PointToKeyInterface::GetDefaultAngularResolution()))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeyClasses::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get point to key function for given name
    //! @param NAME
    //! @param ANGULAR_RESOLUTION
    //! @param DISTANCE_RESOLUTION
    //! @return ShPtr to PointToKey function
    util::ShPtr< PointToKeyInterface>
    PointToKeyClasses::GetPointToKeyFunction
    (
      const std::string &NAME,
      const double ANGULAR_RESOLUTION,
      const double DISTANCE_RESOLUTION
    ) const
    {
      // get enum with that name and make clone of point to key function
      util::ShPtr< PointToKeyInterface> function( ( *EnumType( NAME))->Clone());

      // set resolution
      function->SetAngularResolution( ANGULAR_RESOLUTION);
      function->SetDistanceResolution( DISTANCE_RESOLUTION);

      // end
      return function;
    }

    PointToKeyClasses &GetPointToKeyClasses()
    {
      return PointToKeyClasses::GetEnums();
    }

  } // namespace coord

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< coord::PointToKeyInterface>, coord::PointToKeyClasses>;

  } // namespace util
} // namespace bcl
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
#include "coord/bcl_coord_point_to_key_interface.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {

  //////////
  // data //
  //////////

    //! @brief default angular resolution
    double PointToKeyInterface::GetDefaultAngularResolution()
    {
      return double( 10.0);
    }

    //! @brief default distance resolution
    double PointToKeyInterface::GetDefaultDistanceResolution()
    {
      return double( 2.0);
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! virtual destructor
    PointToKeyInterface::~PointToKeyInterface()
    {
    }

  } // namespace coord
} // namespace bcl
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

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord.h"
#include "coord/bcl_coord_axes.h"
#include "coord/bcl_coord_point_to_key_classes.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PointToKeySpherical
    //! @brief perator will quantize point to three ints according to a Spherical coordinate frame
    //!
    //! @see @link example_coord_point_to_key_spherical.cpp @endlink
    //! @author woetzen
    //! @date Aug 17, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API PointToKeySpherical :
      public PointToKeyInterface
    {
    private:

    //////////
    // data //
    //////////

      static const std::string &GetCoordinateSystemDescriptor();

      //! @brief angular resolution
      double m_AngularResolution;

      //! @brief distance resolution
      double m_DistanceResolution;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! instance of the PointToKey enum
      static PointToKey e_Spherical;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! construct from Resolution
      PointToKeySpherical( const double ANGULAR_RESOLUTION);

      //! copy constructor
      PointToKeySpherical *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the CoordinateSystem as string
      //! @return string describing the coordinate system
      const std::string &GetCoordinateSystem() const;

      //! @brief get angular resolution
      //! @return angular resolution
      double GetAngularResolution() const;

      //! @brief set angular resolution
      //! @param ANGULAR_RESOLUTION for this convert function
      void SetAngularResolution( const double ANGULAR_RESOLUTION);

      //! @brief get distance resolution
      //! @return distance resolution
      double GetDistanceResolution() const;

      //! @brief set distance resolution
      //! @param DISTANCE_RESOLUTION for this convert function
      void SetDistanceResolution( const double DISTANCE_RESOLUTION);

      //! operator will quantize point to three ints within a coordinate frame
      storage::Triplet< int, int, int> operator()( const linal::Vector3D &POINT) const;

      //! operator converting Triplet to a point
      linal::Vector3D operator()( const storage::Triplet< int, int, int> &TRIPLET) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! write to OSTREAM
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! read from ISTREAM
      std::istream &Read( std::istream &ISTREAM);

    }; // class PointToKeySpherical

  //////////
  // data //
  //////////

    const std::string &PointToKeySpherical::GetCoordinateSystemDescriptor()
    {
      static const std::string s_coordinate_system( "Spherical");
      return s_coordinate_system;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointToKeySpherical::s_Instance
    (
      GetObjectInstances().AddInstance( new PointToKeySpherical( PointToKeyInterface::GetDefaultAngularResolution()))
    );

    //! instance of the PointToKey enum
    PointToKey PointToKeySpherical::e_Spherical
    (
      GetPointToKeyClasses().AddEnum( PointToKeySpherical::GetCoordinateSystemDescriptor(), util::ShPtr< PointToKeyInterface>( new PointToKeySpherical( PointToKeyInterface::GetDefaultAngularResolution())))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct from Resolution
    PointToKeySpherical::PointToKeySpherical( const double ANGULAR_RESOLUTION) :
      m_AngularResolution( ANGULAR_RESOLUTION),
      m_DistanceResolution( PointToKeyInterface::GetDefaultDistanceResolution())
    {
    }

    //! virtual copy constructor
    PointToKeySpherical *PointToKeySpherical::Clone() const
    {
      return new PointToKeySpherical( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeySpherical::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the CoordinateSystem as string
    //! @return string describing the coordinate system
    const std::string &PointToKeySpherical::GetCoordinateSystem() const
    {
      // return
      return GetCoordinateSystemDescriptor();
    }

    //! @brief get angular resolution
    //! @return angular resolution
    double PointToKeySpherical::GetAngularResolution() const
    {
      return m_AngularResolution;
    }

    //! @brief set angular resolution
    //! @param ANGULAR_RESOLUTION for this convert function
    void PointToKeySpherical::SetAngularResolution( const double ANGULAR_RESOLUTION)
    {
      m_AngularResolution = ANGULAR_RESOLUTION;
    }

    //! @brief get distance resolution
    //! @return distance resolution
    double PointToKeySpherical::GetDistanceResolution() const
    {
      return m_DistanceResolution;
    }

    //! @brief set distance resolution
    //! @param DISTANCE_RESOLUTION for this convert function
    void PointToKeySpherical::SetDistanceResolution( const double DISTANCE_RESOLUTION)
    {
      m_DistanceResolution = DISTANCE_RESOLUTION;
    }

    //! operator will quantize point to three ints within a coordinate frame
    storage::Triplet< int, int, int> PointToKeySpherical::operator()( const linal::Vector3D &POINT) const
    {
      // derive double values
      const double len( POINT.Norm());                   // length
      const double costheta( POINT.Z() / len);           // costheta
      const double phi( atan2( POINT.Y(), POINT.X()));   // phi
      const double theta( acos( costheta) / math::g_Pi); // theta

      // convert to integers
      const size_t lenint = size_t( len / m_DistanceResolution);
      const size_t thetaint = size_t( theta * m_AngularResolution);
      const size_t points = size_t( sin( math::g_Pi * thetaint / m_AngularResolution) * m_AngularResolution); // number phi points
      const size_t phiint = size_t( ( phi / math::g_Pi) * ( points + 1));

      // convert to hash key
      return storage::Triplet< int, int, int>( int( lenint), int( thetaint), int( phiint));
    }

    //! operator converting Triplet to a point
    linal::Vector3D PointToKeySpherical::operator()( const storage::Triplet< int, int, int> &TRIPLET) const
    {
      // get integer values
      const double lenint( TRIPLET.First());
      const double thetaint( TRIPLET.Second());
      const double phiint( TRIPLET.Third());

      // Calculate angles and distances
      const double len( lenint * m_DistanceResolution);
      const double theta( thetaint / m_AngularResolution);
      const size_t points( size_t( sin( math::g_Pi * thetaint / m_AngularResolution) * m_AngularResolution));
      const double phi( phiint * math::g_Pi / double( points + 1));

      // create the point and set the coordinates
      linal::Vector3D point;
      const double costheta( cos( theta * math::g_Pi));
      const double z( costheta * len);
      const double hypo( math::Sqrt( math::Sqr( len) - math::Sqr( z))); // = r*sin(theta) just faster to calculate that way when z and len is given
      const double x( hypo * cos( phi));
      const double y( hypo * sin( phi));

      //end
      return linal::Vector3D( x, y, z);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write to OSTREAM
    std::ostream &PointToKeySpherical::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AngularResolution, OSTREAM, INDENT);
      io::Serialize::Write( m_DistanceResolution, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read from ISTREAM
    std::istream &PointToKeySpherical::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AngularResolution, ISTREAM);
      io::Serialize::Read( m_DistanceResolution, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace coord
} // namespace bcl
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
#include "coord/bcl_coord_point_to_key_spherical_radius.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace coord
  {
  //////////
  // data //
  //////////

    const std::string &PointToKeySphericalRadius::GetCoordinateSystemDescriptor()
    {
      static const std::string s_coordinate_system( "SphericalRadius");
      return s_coordinate_system;
    }

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PointToKeySphericalRadius::s_Instance
    (
      GetObjectInstances().AddInstance( new PointToKeySphericalRadius( PointToKeyInterface::GetDefaultAngularResolution()))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! construct from Resolution
    PointToKeySphericalRadius::PointToKeySphericalRadius( const double ANGULAR_RESOLUTION) :
      m_AngularResolution( ANGULAR_RESOLUTION),
      m_DistanceResolution( PointToKeyInterface::GetDefaultDistanceResolution())
    {
    }

    //! virtual copy constructor
    PointToKeySphericalRadius *PointToKeySphericalRadius::Clone() const
    {
      return new PointToKeySphericalRadius( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PointToKeySphericalRadius::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the CoordinateSystem as string
    //! @return string describing the coordinate system
    const std::string &PointToKeySphericalRadius::GetCoordinateSystem() const
    {
      // return
      return GetCoordinateSystemDescriptor();
    }

    //! @brief get angular resolution
    //! @return angular resolution
    double PointToKeySphericalRadius::GetAngularResolution() const
    {
      return m_AngularResolution;
    }

    //! @brief set angular resolution
    //! @param ANGULAR_RESOLUTION for this convert function
    void PointToKeySphericalRadius::SetAngularResolution( const double ANGULAR_RESOLUTION)
    {
      m_AngularResolution = ANGULAR_RESOLUTION;
    }

    //! @brief get distance resolution
    //! @return distance resolution
    double PointToKeySphericalRadius::GetDistanceResolution() const
    {
      return m_DistanceResolution;
    }

    //! @brief set distance resolution
    //! @param DISTANCE_RESOLUTION for this convert function
    void PointToKeySphericalRadius::SetDistanceResolution( const double DISTANCE_RESOLUTION)
    {
      m_DistanceResolution = DISTANCE_RESOLUTION;
    }

    //! operator will quantize point to three ints within a coordinate frame
    storage::Triplet< int, int, int> PointToKeySphericalRadius::operator()( const linal::Vector3D &POINT) const
    {
      // derive double values
      const double len( POINT.Norm());                   // length
      const double costheta( POINT.Z() / len);           // costheta
      const double phi( atan2( POINT.Y(), POINT.X()));   // phi
      const double theta( acos( costheta) / math::g_Pi); // theta

      // convert to integers
      const size_t lenint = size_t( math::Sqrt( len / m_DistanceResolution));
      const size_t thetaint = size_t( theta * m_AngularResolution);
      const size_t points = size_t( sin( math::g_Pi * thetaint / m_AngularResolution) * m_AngularResolution); // number phi points
      const size_t phiint = size_t( ( phi / math::g_Pi) * ( points + 1));

      // convert to hash key
      return storage::Triplet< int, int, int>( int( lenint), int( thetaint), int( phiint));
    }

    //! operator converting Triplet to a point
    linal::Vector3D PointToKeySphericalRadius::operator()( const storage::Triplet< int, int, int> &TRIPLET) const
    {
      // get integer values
      const double lenint( TRIPLET.First());
      const double thetaint( TRIPLET.Second());
      const double phiint( TRIPLET.Third());

      // Calculate angles and distances
      const double len( math::Sqr( lenint) * m_DistanceResolution);
      const double theta( thetaint / m_AngularResolution);
      const size_t points( size_t( sin( math::g_Pi * thetaint / m_AngularResolution) * m_AngularResolution));
      const double phi( phiint * math::g_Pi / double( points + 1));

      // create the point and set the coordinates
      linal::Vector3D point;
      const double costheta( cos( theta * math::g_Pi));
      const double z( costheta * len);
      const double hypo( math::Sqrt( math::Sqr( len) - math::Sqr( z))); // = r*sin(theta) just faster to calculate that way when z and len is given
      const double x( hypo * cos( phi));
      const double y( hypo * sin( phi));

      //end
      return linal::Vector3D( x, y, z);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! write to OSTREAM
    std::ostream &PointToKeySphericalRadius::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_AngularResolution, OSTREAM, INDENT);
      io::Serialize::Write( m_DistanceResolution, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! read from ISTREAM
    std::istream &PointToKeySphericalRadius::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_AngularResolution, ISTREAM);
      io::Serialize::Read( m_DistanceResolution, ISTREAM);

      // end
      return ISTREAM;
    }

  } // namespace coord
} // namespace bcl
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
