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
#include "model/bcl_model_dtree_data_partition_function_interface.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_comparisons.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get all split ratings for this partition
    //! @param INPUTDATASET a dataset
    //! @param SHOW_STATUS whether to show current status
    //! @return a vector containing ratings for all columns
    linal::Vector< float> DtreeDataPartitionFunctionInterface::GetAllPartitionRatings
    (
      const storage::Vector< FeatureResultAndState> &INPUTDATASET,
      const bool &SHOW_STATUS
    ) const
    {
      // dataset information
      const size_t descriptor_count( INPUTDATASET.FirstElement().GetFeature().GetSize());

      linal::Vector< float> ratings( descriptor_count, -std::numeric_limits< float>::max());

      // cache # of items in each state on the right side of the partition (initially everything)
      linal::Vector< size_t> total_items_in_state( DetermineTotalStateCounts( INPUTDATASET));
      if( total_items_in_state.GetSize() == 0)
      {
        // pure node, return an undefined split
        return ratings;
      }

      const size_t status_bar_length( 20);

      for
      (
        size_t current_descriptor_index( 0);
        current_descriptor_index < descriptor_count;
        ++current_descriptor_index
      )
      {
        if( SHOW_STATUS)
        {
          // determine progress percent
          const size_t percent( float( current_descriptor_index) * 100.0 / float( descriptor_count));

          // determine number of stars in the status bar
          const size_t number_stars( percent * status_bar_length / 100);

          const std::string status
          (
            "["
            + std::string( number_stars, '*')
            + std::string( status_bar_length - number_stars, ' ')
            + "] "
            + util::Format()( percent) + "% "
            + util::Format()( current_descriptor_index) + "/" + util::Format()( descriptor_count)
            + " columns partitioned"
          );
          util::GetLogger().LogStatus( status);
        }
        DtreeBinaryPartition split_for_this_descriptor
        (
          DetermineSplitValueAndRating
          (
            ExtractAndSortDescriptorAtIndex( current_descriptor_index, INPUTDATASET),
            current_descriptor_index,
            total_items_in_state
          )
        );
        ratings( current_descriptor_index) = split_for_this_descriptor.GetSplitRating();
      }

      return ratings;
    }

    //! @brief extracts a column of descriptors at index INDEX, sorted by index
    //! @param INDEX descriptor index
    //! @param INPUTDATASET ShPtr to a dataset interface
    //! @return list of pairs of descriptor values and a FeatureReference to the results
    //! @note sorts descriptor values will be in ascending order
    storage::Vector< storage::Pair< float, FeatureReference< size_t> > >
      DtreeDataPartitionFunctionInterface::ExtractAndSortDescriptorAtIndex
    (
      const size_t INDEX,
      const storage::Vector< FeatureResultAndState> &INPUTDATASET
    )
    {
      // list of activities and expected results
      storage::Vector< storage::Pair< float, FeatureReference< size_t> > > descriptor_activity_pairs;

      // calculates the data set size once
      const size_t dataset_size( INPUTDATASET.GetSize());
      descriptor_activity_pairs.AllocateMemory( dataset_size);

      // Goes through the dataset and extracts the values for one descriptor and all associated activities into a
      // vector of pairs so that they can be sorted
      for( size_t dataset_counter( 0); dataset_counter < dataset_size; ++dataset_counter)
      {
        descriptor_activity_pairs.PushBack
        (
          storage::Pair< float, FeatureReference< size_t> >
          (
            INPUTDATASET( dataset_counter).GetFeature()( INDEX),
            INPUTDATASET( dataset_counter).GetResultState()
          )
        );
      }

      //! sort pairs by first value (Descriptor) in ascending order
      descriptor_activity_pairs.Sort
      (
        std::less< storage::Pair< float, FeatureReference< size_t> > >()
      );

      // return list of activities and expected results
      return descriptor_activity_pairs;
    } // ExtractAndSortDescriptorAtIndex

    //! @brief determines the best split value and rating
    //! @param SORTED_DESCRIPTORS sorted descriptor values for that feature
    //! @param FEATURE_INDEX index of the feature to split based off of
    //! @param TOTAL_STATE_COUNTS counts of values in each state
    //! @return the best split value and rating
    DtreeBinaryPartition
    DtreeDataPartitionFunctionInterface::DetermineSplitValueAndRating
    (
      const storage::Vector< storage::Pair< float, FeatureReference< size_t> > > &SORTED_DESCRIPTORS,
      const size_t &FEATURE_INDEX,
      const linal::Vector< size_t> &TOTAL_STATE_COUNTS
    ) const
    {
      // cache # of items in each state on the right side of the partition (initially everything)
      if( TOTAL_STATE_COUNTS.GetSize() == 0)
      {
        // only pure states were found, no need for a split
        return DtreeBinaryPartition();
      }

      // object to hold the best split
      DtreeBinaryPartition best_split;

      // get the rating of the complete partition
      const float initial_rating( RatePartition( TOTAL_STATE_COUNTS, SORTED_DESCRIPTORS.GetSize()));

      // iterator for considering candidate descriptor values
      storage::Vector< storage::Pair< float, FeatureReference< size_t> > >::const_iterator
        itr_left( SORTED_DESCRIPTORS.Begin()), itr_right( SORTED_DESCRIPTORS.Begin());

      // move itr_right to the right, so it points to the element on the right side of the partition
      ++itr_right;

      // initialize count of elements on the left and right
      // Initially, all elements are on the right of the partition
      // so initialize the right counts with the total counts
      linal::Vector< size_t> state_counts_left( itr_left->Second());
      linal::Vector< size_t> state_counts_right( TOTAL_STATE_COUNTS - itr_left->Second());

      // keep track of whether the loop just skipped over a segment of descriptors with identical values for which there
      // were multiple states, in which case we need to re-evaluate the partition at the end of the string
      // otherwise, there is no need to re-evaluate the partition until the state changes
      bool found_identical_descriptors_with_different_states( false);

      //iterate over candidate set, which includes boundaries between possible split values,
      // keeping track of the size of left and right hand sides
      for
      (
        size_t left_size( 1), right_size( SORTED_DESCRIPTORS.GetSize() - 1);
        right_size > 0;
        ++itr_left, ++itr_right, ++left_size, --right_size
      )
      {
        // check that the partition values differ
        if( itr_left->First() == itr_right->First())
        {
          // partition fences are the same, check whether the states differed so that we can update
          // found_identical_descriptors_with_different_states
          if( found_identical_descriptors_with_different_states || !( itr_left->Second() == itr_right->Second()))
          {
            found_identical_descriptors_with_different_states = true;
          }
        }
        // partition values differed
        // If the states were identical, or found_identical_descriptors_with_different_states was true, then we need to
        // evaluate this as a possible partition
        else if( found_identical_descriptors_with_different_states || !( itr_left->Second() == itr_right->Second()))
        {
          // reset found_identical_descriptors_with_different_states, since these feature values differed
          found_identical_descriptors_with_different_states = false;

          // calculate partition_rating
          const float partition_rating
          (
            RatePartition( state_counts_left, left_size)
            + RatePartition( state_counts_right, right_size)
            - initial_rating
          );

          // keep best information gain
          if( partition_rating > best_split.GetSplitRating())
          {
            BCL_MessageDbg
            (
              "partition_rating: " + util::Format()( partition_rating)
              + " states: L/R/T/LT/RT: " + util::Format()( left_size) + "/"
              + util::Format()( right_size) + "/"
              + util::Format()( SORTED_DESCRIPTORS.GetSize()) + "/"
              + util::Format()( state_counts_left( 0)) + "/"
              + util::Format()( state_counts_right( 0))
            );
            const float best_split_value( ( itr_left->First() + itr_right->First()) / 2.0);
            best_split = DtreeBinaryPartition( FEATURE_INDEX, best_split_value, partition_rating);
          }
        }

        // update state counts
        state_counts_left += itr_right->Second();
        state_counts_right -= itr_right->Second();
      }

      return best_split;
    }

    //! @brief determines the total state counts
    //! @param DATA dataset
    linal::Vector< size_t> DtreeDataPartitionFunctionInterface::DetermineTotalStateCounts
    (
      const storage::Vector< FeatureResultAndState> &DATA
    )
    {
      // handle the special case off an empty dataset
      if( DATA.IsEmpty())
      {
        return linal::Vector< size_t>();
      }

      // dataset information
      const size_t number_states( DATA.FirstElement().GetResultState().GetSize());
      const size_t total_dataset_elements( DATA.GetSize());

      // cache # of items in each state on the right side of the partition (initially everything)
      linal::Vector< size_t> total_items_in_state( number_states, size_t( 0));
      for( size_t data_counter( 0); data_counter < total_dataset_elements; ++data_counter)
      {
        // add the state vector to the total_items_in_state vector
        total_items_in_state += DATA( data_counter).GetResultState();
      }

      BCL_MessageVrb( "Total items in state: " + util::Format()( total_items_in_state));

      // determine whether this node is already pure, e.g. all items are in a single state
      const size_t number_used_states
      (
        total_items_in_state.GetSize()
        - std::count( total_items_in_state.Begin(), total_items_in_state.End(), size_t( 0))
        - std::count( total_items_in_state.Begin(), total_items_in_state.End(), total_dataset_elements)
      );

      // return if the dataset was already pure
      if( number_used_states == 0)
      {
        return linal::Vector< size_t>();
      }

      return total_items_in_state;
    }

  } // namespace model
} // namespace bcl
