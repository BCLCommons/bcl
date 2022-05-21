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
#include "model/bcl_model_dtree_roc_data_partition_function.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DtreeRocDataPartitionFunction::s_Instance
    (
      util::Enumerated< DtreeDataPartitionFunctionInterface>::AddInstance( new DtreeRocDataPartitionFunction())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! m_MinimumAuROC = 0.5, m_MaximumAuROC = 0.95 according to (Hossain 2006)
    DtreeRocDataPartitionFunction::DtreeRocDataPartitionFunction() :
      m_MinimumAuROC( 0.5),
      m_MaximumAuROC( 0.95)
    {
    }

    //! @brief constructor taking all parameters
    //! @param MAX maximum of allowed ROC AUC
    //! @param MIN minimum of allowed ROC AUC
    DtreeRocDataPartitionFunction::DtreeRocDataPartitionFunction( float MIN, float MAX)
    {
      BCL_Assert( MIN < MAX, "MIN must be less than MAX");
      m_MinimumAuROC = MIN;
      m_MaximumAuROC = MAX;
    }

    //! copy constructor
    DtreeRocDataPartitionFunction *DtreeRocDataPartitionFunction::Clone() const
    {
      return new DtreeRocDataPartitionFunction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DtreeRocDataPartitionFunction::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &DtreeRocDataPartitionFunction::GetAlias() const
    {
      static const std::string s_Name( "ROC");
      return s_Name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator determines the best partition of the data
    //! @param DATA vector of feature result references to consider
    //! @return function returns how to partition to the data
    DtreeBinaryPartition DtreeRocDataPartitionFunction::operator()
    (
      const storage::Vector< FeatureResultAndState> &DATA
    ) const
    {
      // get the total counts of each state
      const linal::Vector< size_t> total_state_counts( DetermineTotalStateCounts( DATA));

      // make sure this node wasn't already pure, if it was, just return the undefined split
      if( total_state_counts.GetSize() == 0)
      {
        return DtreeBinaryPartition();
      }

      // holds the best descriptor and AuROC as a pair and initializes to zero
      DtreeBinaryPartition best_split;
      float best_auroc( 0.0);

      // number of descriptors
      const size_t total_descriptors( DATA.FirstElement().GetFeature().GetSize());

      // number of results
      const size_t number_results( DATA.FirstElement().GetResultState().GetSize());

      // list of descriptor value and result bool for roc curve
      storage::Vector< storage::Pair< float, FeatureReference< size_t> > > descriptor_activity_pairs;

      // loop over all descriptors
      for( size_t current_index( 0); current_index < total_descriptors; ++current_index)
      {
        // list of descriptor value and result bool for roc curve
        descriptor_activity_pairs = ExtractAndSortDescriptorAtIndex( current_index, DATA);

        // compute minimum, maximum, and average roc curve values
        double min_roc( 1.0), max_roc( 0.0);
        math::RunningAverage< double> ave_roc;

        for( size_t result_number( 0); result_number < number_results; ++result_number)
        {
          storage::List< storage::Pair< double, bool> > descriptor_activity_pairs_double;
          // casting to double from float because that's what roc curve wants... will be changed in the future XXX
          for( size_t feature_id( 0), number_features( DATA.GetSize()); feature_id < number_features; ++feature_id)
          {
            descriptor_activity_pairs_double.PushBack();
            // store feature value
            descriptor_activity_pairs_double.LastElement().First() = descriptor_activity_pairs( feature_id).First();
            // store result as bool
            descriptor_activity_pairs_double.LastElement().Second()
              = descriptor_activity_pairs( feature_id).Second()( result_number);
          }

          // calculate AuROC for the current descriptor list and compare it to the saved max AuROC
          math::ROCCurve temp_roc_curve( descriptor_activity_pairs_double);
          // calculate AUC
          float temp_Au_roc_curve( temp_roc_curve.Integral());

          // fpr undefined roc curve values, continue
          if( !util::IsDefined( temp_Au_roc_curve))
          {
            continue;
          }

          // take max of distance from diagonal in ROC curve
          temp_Au_roc_curve = std::max( temp_Au_roc_curve, float( 1) - temp_Au_roc_curve);

          if( temp_Au_roc_curve < min_roc)
          {
            min_roc = temp_Au_roc_curve;
          }
          if( temp_Au_roc_curve > max_roc)
          {
            max_roc = temp_Au_roc_curve;
          }
          ave_roc += temp_Au_roc_curve;
        }
        // save a new descriptor with the new max AuROC
        if( ave_roc.GetAverage() > best_auroc && max_roc <= m_MaximumAuROC && min_roc >= m_MinimumAuROC)
        {
          best_split =
            DetermineSplitValueAndRating
            (
              descriptor_activity_pairs,
              current_index,
              total_state_counts
            );
          best_auroc = ave_roc.GetAverage();
        }
      }

      // return best descriptor auc
      BCL_MessageDbg( "best split is :" + util::Format()( best_split));

      // determine accuracy, which may be used to select which node to expand next
      best_split.DetermineAccuracy( DATA);

      // return the best split
      return best_split;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief rates by gini-index (denormalized
    //! @param STATE_COUNTS counts of each state (e.g. true, false, or 0,1,2, etc.) on one segment of the partition
    //! @param SIZE total number of items in the partition
    //! @return the # of correct classifications on this side
    float DtreeRocDataPartitionFunction::RatePartition
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

      // return the maximum value for this partition (the max # of correct classifications)
      // this was the method chosen in the hussan paper, unfortunately it results in all the actives being split into
      // separate nodes, resulting in a tree roughly 2x the size of the # of actives
      // It is ONLY appropriate for categorical data, it utterly fails when numerical data is used
      //return math::Statistics::MaximumValue( STATE_COUNTS.Begin(), STATE_COUNTS.End());
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DtreeRocDataPartitionFunction::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses information gain based on the receiver operator curve to decide which feature index to make the next decision with"
      );

      parameters.AddInitializer
      (
        "min area under ROC",
        "minimum area under the ROC curve before termination of iterations",
        io::Serialization::GetAgent( &m_MinimumAuROC),
        "0.5"
      );

      parameters.AddInitializer
      (
        "max area under ROC",
        "maximum area under the ROC curve (AuROC) before termination of iterations, based on Hossain paper to avoid overfitting",
        io::Serialization::GetAgent( &m_MaximumAuROC),
        "0.95"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
