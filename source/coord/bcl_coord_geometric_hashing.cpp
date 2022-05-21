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
