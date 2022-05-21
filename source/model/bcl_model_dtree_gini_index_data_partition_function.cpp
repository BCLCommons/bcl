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
#include "model/bcl_model_dtree_gini_index_data_partition_function.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DtreeGiniIndexDataPartitionFunction::s_Instance
    (
      util::Enumerated< DtreeDataPartitionFunctionInterface>::AddInstance( new DtreeGiniIndexDataPartitionFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! copy constructor
    DtreeGiniIndexDataPartitionFunction *DtreeGiniIndexDataPartitionFunction::Clone() const
    {
      return new DtreeGiniIndexDataPartitionFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DtreeGiniIndexDataPartitionFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &DtreeGiniIndexDataPartitionFunction::GetAlias() const
    {
      static const std::string s_Name( "Gini");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator determines the best partition of the data
    //! @param DATA vector of feature result references to consider
    //! @return function returns how to partition to the data
    DtreeBinaryPartition DtreeGiniIndexDataPartitionFunction::operator()
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

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DtreeGiniIndexDataPartitionFunction::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses gini index to decide which feature index to make the next decision based on see http://en.wikipedia.org/wiki/Decision_tree_learning#Gini_impurity"
      );
      return parameters;
    }

    //! @brief calculates information gain from counts of states in each segment
    //! @param STATE_COUNTS counts of each state (e.g. true, false, or 0,1,2, etc.) on one segment of the partition
    //! @param SIZE total number of items in the partition
    //! @return the information gain metric (unnormalized)
    float DtreeGiniIndexDataPartitionFunction::RatePartition
    (
      const linal::Vector< size_t> &STATE_COUNTS,
      const size_t &SIZE
    ) const
    {
      // compute the gini index
      float sum_gini( 0.0);
      for( const size_t *itr( STATE_COUNTS.Begin()), *itr_end( STATE_COUNTS.End()); itr != itr_end; ++itr)
      {
        sum_gini += math::Sqr( float( *itr) / float( SIZE));
      }
      return sum_gini * SIZE;
    }

  } // namespace model
} // namespace bcl
