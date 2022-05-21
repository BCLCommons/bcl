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
#include "quality/bcl_quality_maxsub.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "quality/bcl_quality_rmsd.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace quality
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MaxSub::s_Instance
    (
      GetObjectInstances().AddInstance( new MaxSub())
    );

    //! @brief returns default RMSD cutoff
    //! @return default RMSD cutoff
    double MaxSub::GetDefaultRMSDCutoff()
    {
      return 3.5;
    }

    //! @brief returns default length of the seed subset of coordinates
    //! @return default seed length of the seed subset of coordinates
    size_t MaxSub::GetDefaultSeedLength()
    {
      return 3;
    }

    //! @brief returns default number of iterations to extend an individual seed subset
    //! @return default number of iterations to extend an individual seed subset
    size_t MaxSub::GetDefaultNumberIterations()
    {
      return 4;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a RMSD cutoff, seed length and number of iterations
    //! @param RMSD_CUTOFF RMSD cutoff for the longest subset of superimposed coordinates
    //! @param SEED_LENGTH length of the seed subset of coordinates
    //! @param NUMBER_ITERATIONS number of iterations to extend an individual seed subset
    MaxSub::MaxSub
    (
      const double RMSD_CUTOFF,
      const size_t SEED_LENGTH,
      const size_t NUMBER_ITERATIONS
    ) :
      m_RMSDCutoff( RMSD_CUTOFF),
      m_SeedLength( SEED_LENGTH),
      m_NumberIterations( NUMBER_ITERATIONS)
    {
      // check the given seed length
      BCL_Assert
      (
        SEED_LENGTH >= 3,
        "The seed length has to be at least 3 for transformations to work; not " + util::Format()( SEED_LENGTH)
      );
    }

    //! @brief Clone function
    //! @return pointer to new MaxSub
    MaxSub *MaxSub::Clone() const
    {
      return new MaxSub( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MaxSub::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &MaxSub::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief find a larger subset by extending the given one that has a RMSD below the cutoff
    //! @param SUBSET subset of coordinates that are used to be as seed and to be extended
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return transformation for the largest subset
    math::TransformationMatrix3D MaxSub::ExtendSubset
    (
      storage::List< size_t> &SUBSET,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());

      // make a copy of the seed subset
      const storage::List< size_t> seed_subset( SUBSET);

      // iterate over coordinates to extend the subset
      for( size_t nr_itr( 1); nr_itr <= m_NumberIterations; ++nr_itr)
      {
        // double check that there are at least 3 members in the subset
        if( SUBSET.GetSize() < 3)
        {
          BCL_MessageVrb
          (
            "the size of subset is smaller than 3: " + util::Format()( SUBSET.GetSize()) +
            " skipping more iterations"
          )
          return math::TransformationMatrix3D();
        }

        // store allowed distance cutoff and its squared value
        const double this_distance_cutoff( double( nr_itr) * m_RMSDCutoff / double( m_NumberIterations));
        const double distance_cutoff_squared( this_distance_cutoff * this_distance_cutoff);

        // superimpose the coordinates for the subset of coordinates
        const math::TransformationMatrix3D transformation
        (
          RMSD::SuperimposeCoordinates
          (
            CollectCoordinatesSubset( SUBSET, REFERENCE_COORDINATES),
            CollectCoordinatesSubset( SUBSET, COORDINATES)
          )
        );

        // if the transformation was not defined
        if( !transformation.IsDefined())
        {
          BCL_MessageVrb
          (
            "undefined transformation in Maxsub calculation at iteration " + util::Format()( nr_itr) +
            " skipping more iterations"
          )
          return transformation;
        }

        // initialize the extended subset to be returned
        storage::List< size_t> extended_subset;

        // now iterate over all residues
        for( size_t index( 0); index < number_points; ++index)
        {
          // if it's already
          if( index >= seed_subset.FirstElement() && index <= seed_subset.LastElement())
          {
            extended_subset.PushBack( index);
            continue;
          }

          // make a copy of the coordinates from COORDINATES_B
          linal::Vector3D coordinate_copy( *COORDINATES( index));

          // calculate the distance squared between the coordinate A and transformed coordinate B
          const double this_distance_squared
          (
            ( ( *REFERENCE_COORDINATES( index)) - coordinate_copy.Transform( transformation)).SquareNorm()
          );

          // if the distance between the transformed coordinates is smaller than the RMSD cutoff
          if( this_distance_squared < distance_cutoff_squared)
          {
            // insert it into the new subset
            extended_subset.PushBack( index);
          }
        } // end index

        // update the subset the extended subset
        SUBSET = extended_subset;

      } // end iterations

      // now recompute the transformation with the latest SUBSET
      // superimpose the coordinates for the subset of coordinates
      const math::TransformationMatrix3D final_transformation
      (
        RMSD::SuperimposeCoordinates
        (
          CollectCoordinatesSubset( SUBSET, REFERENCE_COORDINATES),
          CollectCoordinatesSubset( SUBSET, COORDINATES)
        )
      );

      // create a new subset
      storage::List< size_t> filtered_subset;

      // calculate the cutoff square
      const double rmsd_cutoff_squared( m_RMSDCutoff * m_RMSDCutoff);

      // iterate over the indices in SUBSET
      for
      (
        storage::List< size_t>::iterator index_itr( SUBSET.Begin()), index_itr_end( SUBSET.End());
        index_itr != index_itr_end;
        ++index_itr
      )
      {
        // make a copy of the coordinates from COORDINATES_B
        linal::Vector3D coordinate_copy( *COORDINATES( *index_itr));

        // calculate the distance squared between coordinate A and transformed coordinate B
        const double this_distance_squared
        (
          ( ( *REFERENCE_COORDINATES( *index_itr)) - coordinate_copy.Transform( final_transformation)).SquareNorm()
        );

        // if the distance is not greater than the cutoff
        if( this_distance_squared <= ( rmsd_cutoff_squared))
        {
          // add it to filtered subset
          filtered_subset.PushBack( *index_itr);
        }
      }

      // assign the filtered subset
      SUBSET = filtered_subset;

      // return final transformation
      return final_transformation;
    }

    //! @brief collects the subset of coordinates specified by the given list of indices
    //! @param SUBSET indices of coordinates that define the subset of coordinates to be collected
    //! @param COORDINATES util::SiPtrVector< Vector3D> which will be measured for agreement with COORDINATES_B
    //! @return vector of coordinates that correspond to the requested subset
    util::SiPtrVector< const linal::Vector3D> MaxSub::CollectCoordinatesSubset
    (
      const storage::List< size_t> &SUBSET,
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES
    ) const
    {
      // initialize SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates_subset;
      coordinates_subset.AllocateMemory( SUBSET.GetSize());

      // iterate over the subset
      for
      (
        storage::List< size_t>::const_iterator index_itr( SUBSET.Begin()), index_itr_end( SUBSET.End());
        index_itr != index_itr_end;
        ++index_itr
      )
      {
        // insert the coordinates at the given index to the subset to be returned
        coordinates_subset.PushBack( COORDINATES( *index_itr));
      }

      // end
      return coordinates_subset;
    }

    //! @brief calculates MaxSub between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return MaxSub between given coordinates
    double MaxSub::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateMaxSubAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).First();
    }

    //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D MaxSub::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateMaxSubAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).Second();
    }

    //! @brief calculates the MaxSub and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return pair of the MaxSub and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    storage::Pair< double, math::TransformationMatrix3D> MaxSub::CalculateMaxSubAndSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // make sure the sizes of the coordinates are the same
      if( COORDINATES.GetSize() != REFERENCE_COORDINATES.GetSize())
      {
        BCL_MessageStd( "The number of coordinates given to MaxSub calculation differ!");
        return
          storage::Pair< double, math::TransformationMatrix3D>
          (
            util::GetUndefinedDouble(), math::TransformationMatrix3D()
          );
      }

      // initialize the variable to hold the indices for the largest subset found so far
      storage::List< size_t> largest_subset;

      // initialize best transformation
      math::TransformationMatrix3D best_transformation;

      // store number of points and number of seed fragments
      const size_t number_points( COORDINATES.GetSize());
      const size_t number_seed_subsets( number_points - m_SeedLength + 1);

      // iterate over coordinates to create seed segments
      for( size_t index( 0); index < number_seed_subsets; ++index)
      {
        // initialize seed subset
        storage::List< size_t> seed_subset;

        // form the seed subset
        for( size_t seed_index( index); seed_index < index + m_SeedLength; ++seed_index)
        {
          seed_subset.PushBack( seed_index);
        }

        // extend the seed and store the transformation
        const math::TransformationMatrix3D transformation
        (
          ExtendSubset( seed_subset, COORDINATES, REFERENCE_COORDINATES)
        );

        // if the extended subset leads to size smaller than 3
        if( seed_subset.GetSize() < 3)
        {
          continue;
        }

        // if the length of the extended seed is larger then the longest seed so far
        if( seed_subset.GetSize() > largest_subset.GetSize())
        {
          // update this
          largest_subset = seed_subset;

          // also update the best transformation
          best_transformation = transformation;
        }
      }
      // calculate MaxSub value
      const double maxsub( OptimalValue() * double( largest_subset.GetSize()) / double( number_points));

      // return the MaxSub value and the corresponding transformation matrix
      return storage::Pair< double, math::TransformationMatrix3D>( maxsub, best_transformation);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MaxSub::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RMSDCutoff      , ISTREAM);
      io::Serialize::Read( m_SeedLength      , ISTREAM);
      io::Serialize::Read( m_NumberIterations, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MaxSub::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RMSDCutoff      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeedLength      , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberIterations, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace quality
} // namespace bcl
