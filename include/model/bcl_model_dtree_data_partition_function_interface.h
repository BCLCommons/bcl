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

#ifndef BCL_MODEL_DTREE_DATA_PARTITION_FUNCTION_INTERFACE_H_
#define BCL_MODEL_DTREE_DATA_PARTITION_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_dtree_binary_partition.h"
#include "bcl_model_feature_reference.h"
#include "bcl_model_feature_result_and_state.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DtreeDataPartitionFunctionInterface
    //! @brief is a Function Interface class for defining
    //!        the individual tree splitting technique. Specific Implementations are
    //!        ROC-tree splitting, and splitting according to information gain.
    //!
    //! @remarks example unnecessary
    //!
    //! @see @link example_model_dtree_data_partition_function_interface.h @endlink
    //! @author lemmonwa
    //! @date 07/07/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DtreeDataPartitionFunctionInterface :
      public util::SerializableInterface
    {

    public:

    ///////////////////////////////////
    // construction and desctruction //
    ///////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new DtreeBinaryPartition
      virtual DtreeDataPartitionFunctionInterface *Clone() const = 0;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator determines the best partition of the data
      //! @param DATA vector of feature result references to consider
      //! @return function returns how to partition to the data
      virtual DtreeBinaryPartition operator()
      (
        const storage::Vector< FeatureResultAndState> &DATA
      ) const = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get all split ratings for this partition
      //! @param INPUTDATASET a dataset
      //! @param SHOW_STATUS whether to show current status
      //! @return a vector containing ratings for all columns
      virtual linal::Vector< float> GetAllPartitionRatings
      (
        const storage::Vector< FeatureResultAndState> &INPUTDATASET,
        const bool &SHOW_STATUS = false
      ) const;

      //! @brief extracts a column of descriptors at index INDEX, sorted by index
      //! @param INDEX descriptor index
      //! @param INPUTDATASET a dataset
      //! @return list of pairs of descriptor values and a FeatureReference to the results
      //! @note sorts descriptor values will be in ascending order
      static storage::Vector< storage::Pair< float, FeatureReference< size_t> > >
      ExtractAndSortDescriptorAtIndex
      (
        const size_t INDEX,
        const storage::Vector< FeatureResultAndState> &INPUTDATASET
      );

    protected:

      //! @brief determines the best split value and rating
      //! @param SORTED_DESCRIPTORS sorted descriptor values for that feature
      //! @param FEATURE_INDEX index of the feature to split based off of
      //! @param TOTAL_STATE_COUNTS counts of values in each state
      //! @return the best split value and rating
      DtreeBinaryPartition DetermineSplitValueAndRating
      (
        const storage::Vector< storage::Pair< float, FeatureReference< size_t> > > &SORTED_DESCRIPTORS,
        const size_t &FEATURE_INDEX,
        const linal::Vector< size_t> &TOTAL_STATE_COUNTS
      ) const;

      //! @brief determines the total state counts
      //! @param DATA dataset
      static linal::Vector< size_t> DetermineTotalStateCounts
      (
        const storage::Vector< FeatureResultAndState> &DATA
      );

      //! @brief abstract function to rate a partition (higher ratings should be better)
      //! @param STATE_COUNTS counts of each state (e.g. true, false, or 0,1,2, etc.) on one segment of the partition
      //! @param SIZE total number of items in the partition
      //! @return the rating
      virtual float RatePartition
      (
        const linal::Vector< size_t> &STATE_COUNTS,
        const size_t &SIZE
      ) const = 0;

    }; // class DtreeDataPartitionFunctionInterface

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_DTREE_DATA_PARTITION_FUNCTION_INTERFACE_H_
