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
#include "model/bcl_model_dtree_information_gain_data_partition_function.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DtreeInformationGainDataPartitionFunction::s_Instance
    (
      util::Enumerated< DtreeDataPartitionFunctionInterface>::AddInstance( new DtreeInformationGainDataPartitionFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! copy constructor
    DtreeInformationGainDataPartitionFunction *DtreeInformationGainDataPartitionFunction::Clone() const
    {
      return new DtreeInformationGainDataPartitionFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DtreeInformationGainDataPartitionFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &DtreeInformationGainDataPartitionFunction::GetAlias() const
    {
      static const std::string s_Name( "InformationGain");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator determines the best partition of the data
    //! @param DATA vector of feature result references to consider
    //! @return function returns how to partition to the data
    DtreeBinaryPartition DtreeInformationGainDataPartitionFunction::operator()
    (
      const storage::Vector< FeatureResultAndState> &DATA
    ) const
    {
      // dataset information
      const size_t descriptor_count( DATA.FirstElement().GetFeature().GetSize());

      // cache # of items in each state on the right side of the partition (initially everything)
      linal::Vector< size_t> total_items_in_state( DetermineTotalStateCounts( DATA));
      if( total_items_in_state.GetSize() == 0)
      {
        // pure node, return an undefined split
        return DtreeBinaryPartition();
      }

      // object to hold the best split
      DtreeBinaryPartition best_split;

      for
      (
        size_t current_descriptor_index( 0);
        current_descriptor_index < descriptor_count;
        ++current_descriptor_index
      )
      {
        DtreeBinaryPartition split_for_this_descriptor
        (
          DetermineSplitValueAndRating
          (
            ExtractAndSortDescriptorAtIndex( current_descriptor_index, DATA),
            current_descriptor_index,
            total_items_in_state
          )
        );

        BCL_MessageVrb
        (
          "split_for_this_descriptor index,value,info gain: "
          + util::Format()( split_for_this_descriptor.GetFeatureIndex()) + ","
          + util::Format()( split_for_this_descriptor.GetSplitValue()) + ","
          + util::Format()( split_for_this_descriptor.GetSplitRating())
        );
        if
        (
          util::IsDefined( split_for_this_descriptor.GetFeatureIndex()) &&
          util::IsDefined( split_for_this_descriptor.GetSplitValue()) &&
          util::IsDefined( split_for_this_descriptor.GetSplitRating()) &&
          split_for_this_descriptor.GetSplitRating() > best_split.GetSplitRating()
        )
        {
          best_split = split_for_this_descriptor;
        }
      }

      // determine accuracy, which may be used to select which node to expand next
      best_split.DetermineAccuracy( DATA);
      return best_split;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DtreeInformationGainDataPartitionFunction::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses raw information gain to decide which feature index to make the next decision based on see http://en.wikipedia.org/wiki/Information_gain"
      );
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates information gain from counts of states in each segment
    //! @param STATE_COUNTS counts of each state (e.g. true, false, or 0,1,2, etc.) on one segment of the partition
    //! @param SIZE total number of items in the partition
    //! @return the information gain metric (unnormalized)
    float DtreeInformationGainDataPartitionFunction::RatePartition
    (
      const linal::Vector< size_t> &STATE_COUNTS,
      const size_t &SIZE
    ) const
    {
      // initialize a sum for the information gain of each state
      float information_gain( 0.0);

      // get the total number of elements in that state
      const float size_as_float( SIZE);

      for
      (
        auto itr( STATE_COUNTS.Begin()), itr_end( STATE_COUNTS.End());
        itr != itr_end;
        ++itr
      )
      {
        // only calculate information gain for states that were present in this partition
        if( *itr != 0 && *itr != SIZE) // information gain for presence of state
        {
          // information gain for each state (S) is
          // count( items in S) * log( probability of item in S in this partition)
          information_gain += *itr * std::log( *itr / size_as_float);

          // information gain for absence of state
          const float count_not_in_state( SIZE - *itr);
          information_gain += count_not_in_state * std::log( count_not_in_state / size_as_float);
        }
      }

      // return the information gain for this partition
      return information_gain;
    }

  } // namespace model
} // namespace bcl
