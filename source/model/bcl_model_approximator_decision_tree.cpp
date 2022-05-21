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
#include "model/bcl_model_approximator_decision_tree.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ApproximatorDecisionTree::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorDecisionTree())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ApproximatorDecisionTree::ApproximatorDecisionTree() :
      m_DataPartitioner(),
      m_ActivityCutoff(),
      m_UnexpandedSubTrees(),
      m_DecisionTree( new DecisionTree()),
      m_MinIncorrect( 0)
    {
    }

    //! @brief constructor from all necessary parameters
    //! @param DATA_PARTITIONER ShPtr to DtreeDataPartitionFunctionInterface
    //! @param OBJECTIVE_FUNCTION ShPtr to objective function
    //! @param ACTIVITY_CUTOFF the cutoff for binary classification
    //! @param TRAINING_DATA the dataset used to build the tree
    ApproximatorDecisionTree::ApproximatorDecisionTree
    (
      const util::Implementation< DtreeDataPartitionFunctionInterface> &DATA_PARTITIONER,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
      const float ACTIVITY_CUTOFF,
      const util::ShPtr< descriptor::Dataset> &TRAINING_DATA
    ) :
      ApproximatorBase( OBJECTIVE_FUNCTION),
      m_DataPartitioner( DATA_PARTITIONER),
      m_ActivityCutoff( ACTIVITY_CUTOFF),
      m_PartitionScoreType( DtreeBinaryPartition::e_SplitRating),
      m_UnexpandedSubTrees(),
      m_DecisionTree( new DecisionTree()),
      m_MinIncorrect( 0)
    {
      if( TRAINING_DATA.IsDefined())
      {
        util::ShPtr< descriptor::Dataset> training_data( TRAINING_DATA);
        SetTrainingData( training_data);
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorDecisionTree::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorDecisionTree::GetAlias() const
    {
      static const std::string s_Name( "DecisionTree");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorDecisionTree::SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA)
    {
      m_TrainingData = DATA;

      //! set up processed data set for future use, contains bools for activities
      util::ShPtr< storage::Vector< FeatureResultAndState> > processed_dataset
      (
        CreateDataSetReferences
        (
          m_TrainingData->GetFeaturesPtr()->GetMatrix(),
          m_TrainingData->GetResultsPtr()->GetMatrix()
        )
      );

      m_UnexpandedSubTrees.Reset();

      // determine the root node's activity and add it as the next node to expand.
      AddExpandableNode( m_DecisionTree, processed_dataset);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorDecisionTree::GetCurrentModel() const
    {
      if( m_UnexpandedSubTrees.IsEmpty())
      {
        // no more iteration to do, so just return the final decision tree
        return m_DecisionTree;
      }

      // more iteration to do, need to deeply clone the tree
      util::ShPtr< DecisionTree> deep_cloned_tree( m_DecisionTree->DeepClone());
      Prune( *deep_cloned_tree);
      return deep_cloned_tree;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorDecisionTree::GetCurrentApproximation() const
    {
      util::ShPtr< Interface> model( GetCurrentModel());

      // return the model & objective function value
      return
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>( model, m_ObjectiveFunction->operator()( model))
        );
    }

    //! @brief Prune a decision tree model based on size & cutoff (if applicable)
    //! @param MODEL the decision tree to prune
    void ApproximatorDecisionTree::Prune( DecisionTree &MODEL) const
    {
      // prune the tree if the objective function is a strict cutoff-based objective function
      if( m_ObjectiveFunction->GetGoalType() == ObjectiveFunctionInterface::e_Classification)
      {
        MODEL.Prune( m_ActivityCutoff, m_MinIncorrect);
      }
      else
      {
        // otherwise, just prune by size
        MODEL.Prune( util::GetUndefined< double>(), m_MinIncorrect);
      }
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorDecisionTree::Next()
    {
      // if there are still unexpanded nodes left
      if( !m_UnexpandedSubTrees.IsEmpty())
      {
        // expand the node with the best partition rating
        storage::Triplet
        <
          util::ShPtr< DecisionTree>,                             // the tree to be expanded
          util::ShPtr< storage::Vector< FeatureResultAndState> >, // references of feature results to expand
          DtreeBinaryPartition                                    // best binary partition for this dataset
        > best_partition_info( m_UnexpandedSubTrees.LastElement());

        // remove that element from the unexpanded subtrees list
        m_UnexpandedSubTrees.PopBack();

        // expand the given nodee_Debug
        Expand( best_partition_info.First(), best_partition_info.Second(), best_partition_info.Third());

        if( m_UnexpandedSubTrees.IsEmpty())
        {
          // last node was just expanded, prune the tree, if necessary
          Prune( *m_DecisionTree);
        }
      }

      this->GetTracker().Track( GetCurrentApproximation());
    } // Next

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorDecisionTree::Read( std::istream &ISTREAM)
    {
      // read members
      ApproximatorBase::Read( ISTREAM);
      io::Serialize::Read( m_DataPartitioner, ISTREAM);
      io::Serialize::Read( m_ActivityCutoff, ISTREAM);
      io::Serialize::Read( m_PartitionScoreType, ISTREAM);
      io::Serialize::Read( m_UnexpandedSubTrees, ISTREAM);
      io::Serialize::Read( m_DecisionTree, ISTREAM);
      io::Serialize::Read( m_TrainingData, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorDecisionTree::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      ApproximatorBase::Write( OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DataPartitioner, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ActivityCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PartitionScoreType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UnexpandedSubTrees, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DecisionTree, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TrainingData, OSTREAM, INDENT);
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorDecisionTree::CanContinue() const
    {
      return !m_UnexpandedSubTrees.IsEmpty();
    }

    //! @brief constructs children nodes, adds them to current node, and queues children for visiting
    //! @param CURRENT_NODE a ShPtr to the node to expand
    //! @param DATA data from the parent node
    //! @param BINARY_PARTITION partition to use
    void ApproximatorDecisionTree::Expand
    (
      util::ShPtr< DecisionTree> &CURRENT_NODE,
      const util::ShPtr< storage::Vector< FeatureResultAndState> > &DATA,
      const DtreeBinaryPartition &BINARY_PARTITION
    )
    {
      //! makes the new node for the returned decision function/dataset pair and sets the activity for the newly
      //! created node based on the accompanying dataset
      util::ShPtr< DecisionTree> left_tree_sp( new DecisionTree()), right_tree_sp( new DecisionTree());

      // partition the dataset into the data that would be attached to the left node and data that would be attached
      // to the right node
      util::ShPtr< storage::Vector< FeatureResultAndState> > left_data_set( new storage::Vector< FeatureResultAndState>);
      util::ShPtr< storage::Vector< FeatureResultAndState> > right_data_set( new storage::Vector< FeatureResultAndState>);

      const size_t split_index( BINARY_PARTITION.GetFeatureIndex());
      const float split_value( BINARY_PARTITION.GetSplitValue());

      // get a reference to the data
      const storage::Vector< FeatureResultAndState> &data( *DATA);

      //! calculates the data set size once
      const size_t total_data_set_elements( data.GetSize());

      BCL_MessageVrb
      (
        "Splitting index: " + util::Format()( split_index)
        + " total elements: " + util::Format()( total_data_set_elements)
        + " split ranking: " + util::Format()( BINARY_PARTITION.GetSplitRating())
      );

      //! Goes through the dataset and partitions it into the right and left datasets as appropriate
      for( size_t dataset_counter( 0); dataset_counter < total_data_set_elements; ++dataset_counter)
      {
        if( data( dataset_counter).GetFeature()( split_index) <= split_value)
        {
          left_data_set->PushBack( data( dataset_counter));
        }
        else
        {
          right_data_set->PushBack( data( dataset_counter));
        }
      }

      // if either partitions are empty, then there is no point in expanding the tree, so return
      if( left_data_set->IsEmpty() || right_data_set->IsEmpty())
      {
        return;
      }

      CURRENT_NODE->ConfigureBranch( split_index, split_value, left_tree_sp, right_tree_sp);

      if( BINARY_PARTITION.GetFinalNumIncorrect() >= m_MinIncorrect)
      {
        AddExpandableNode( left_tree_sp, left_data_set);
        AddExpandableNode( right_tree_sp, right_data_set);
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorDecisionTree::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Trains a decision tree using an arbitrary method of choosing which feature index to partition next"
      );

      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each batch step",
        io::Serialization::GetAgent( &GetObjectiveFunction()->GetImplementation()),
        "Accuracy"
      );

      parameters.AddInitializer
      (
        "partitioner",
        "method of deciding which component of the feature vector to use to add a new node to the decision tree",
        io::Serialization::GetAgent( &m_DataPartitioner),
        "InformationGain"
      );

      parameters.AddInitializer
      (
        "activity cutoff",
        "result value that separates correct from incorrect results",
        io::Serialization::GetAgent( &m_ActivityCutoff),
        "0.5"
      );
      parameters.AddInitializer
      (
        "node score",
        "choice of methods by which to pick the next node to expand",
        io::Serialization::GetAgent( &m_PartitionScoreType),
        "SplitRating"
      );
      parameters.AddInitializer
      (
        "min split",
        "Minimum number of incorrect classifications in a node that will be split",
        io::Serialization::GetAgent( &m_MinIncorrect),
        "0"
      );

      return parameters;
    }

    //! @brief converts an original training data set into an already classified version based on a given cutoff
    //! @param FEATURES original features of interest
    //! @param RESULTS corresponding results
    //! @return pre-processed data set
    util::ShPtr< storage::Vector< FeatureResultAndState> > ApproximatorDecisionTree::CreateDataSetReferences
    (
      const linal::MatrixConstInterface< float> &FEATURES,
      const linal::MatrixConstInterface< float> &RESULTS
    )
    {
      const size_t results_size( RESULTS.GetNumberCols());
      const size_t features_size( FEATURES.GetNumberCols());
      const size_t dataset_size( RESULTS.GetNumberRows());

      // setup the result state matrix
      m_ResultStates = linal::Matrix< size_t>( dataset_size, results_size);

      // create a preprocessed dataset from original dataset format
      util::ShPtr< storage::Vector< FeatureResultAndState> > sp_converted_dataset
      (
        new storage::Vector< FeatureResultAndState>( dataset_size)
      );
      storage::Vector< FeatureResultAndState> &converted_dataset( *sp_converted_dataset);
      storage::Vector< math::RunningAverage< float> > averages_above( results_size);
      storage::Vector< math::RunningAverage< float> > averages_below( results_size);

      // calculates the data set size once
      //!set up preprocessed data set for future use, contains bools for activities
      for( size_t dataset_counter( 0); dataset_counter < dataset_size; ++dataset_counter)
      {
        const float *results_row( RESULTS[ dataset_counter]);
        size_t *result_states_row( m_ResultStates[ dataset_counter]);
        // setup the row of m_ResultStates first
        for( size_t result_counter( 0); result_counter < results_size; ++result_counter)
        {
          result_states_row[ result_counter] = results_row[ result_counter] <= m_ActivityCutoff;
          if( result_states_row[ result_counter])
          {
            averages_below( result_counter) += results_row[ result_counter];
          }
          else
          {
            averages_above( result_counter) += results_row[ result_counter];
          }
        }

        converted_dataset( dataset_counter) =
          FeatureResultAndState
          (
            FeatureReference< float>(  features_size, FEATURES[ dataset_counter]),
            FeatureReference< float>(  results_size,  results_row),
            FeatureReference< size_t>( results_size,  m_ResultStates[ dataset_counter])
          );
      }

      m_AverageAboveCutoffSlopes = linal::Vector< float>( results_size, 0.0);
      m_AverageBelowCutoffSlopes = linal::Vector< float>( results_size, 0.0);
      for( size_t result_counter( 0); result_counter < results_size; ++result_counter)
      {
        m_AverageAboveCutoffSlopes( result_counter) = averages_above( result_counter).GetAverage() - m_ActivityCutoff;
        m_AverageBelowCutoffSlopes( result_counter) = m_ActivityCutoff - averages_below( result_counter).GetAverage();
      }

      // return preprocessed dataset
      return sp_converted_dataset;
    }

    //! @brief add a node, determines average activity, etc., sorts it into the list of nodes to expand
    //! @param NODE the node to add
    //! @param DATA the data that maps to that node
    void ApproximatorDecisionTree::AddExpandableNode
    (
      util::ShPtr< DecisionTree> NODE,
      const util::ShPtr< storage::Vector< FeatureResultAndState> > &DATA
    )
    {
      // for empty dataset, just reset the activity level of the tree
      if( DATA->IsEmpty())
      {
        NODE->SetActivity( storage::Vector< math::RunningAverage< float> >());
        return;
      }

      // determine # of results/outputs
      const size_t results_size( DATA->FirstElement().GetResult().GetSize());

      if
      (
        m_ObjectiveFunction->GetGoalType() != ObjectiveFunctionInterface::e_Classification
        && m_ObjectiveFunction->GetGoalType() != ObjectiveFunctionInterface::e_RankClassification
      )
      {
        // a regression (or other) objective function type.  Set the activity up to be the actual average
        // set up the activity to the raw average; this is the optimal value for regression-style objective functions

        storage::Vector< math::RunningAverage< float> > average_values( results_size);

        for
        (
          storage::Vector< FeatureResultAndState>::const_iterator itr( DATA->Begin()), itr_end( DATA->End());
          itr != itr_end;
          ++itr
        )
        {
          // get a reference on the result and its state
          const FeatureReference< float> &result_ref( itr->GetResult());
          for( size_t result_number( 0); result_number < results_size; ++result_number)
          {
            average_values( result_number) += result_ref( result_number);
          }
        }

        NODE->SetActivity( average_values);
      }
      else
      {
        // sum up the total in each state
        linal::Vector< size_t> state_counts( results_size, size_t( 0));
        const size_t number_results( DATA->GetSize());
        for
        (
          storage::Vector< FeatureResultAndState>::const_iterator itr( DATA->Begin()), itr_end( DATA->End());
          itr != itr_end;
          ++itr
        )
        {
          const FeatureReference< size_t> &state_ref( itr->GetResultState());
          for( size_t result_number( 0); result_number < results_size; ++result_number)
          {
            state_counts( result_number) += state_ref( result_number);
          }
        }

        // determine the results average value for categorical values such that the value reflects the purity of the node
        // That is, the there are exactly as many counts in state x as not in state x, use the cutoff
        // otherwise use the slopes and purity to give a value between cutoff and the average result value of all states
        // that far away from the cutoff according to the following function:
        // ave_state( x) = state_counts( x) / size
        // f(x) = cutoff + ( 1 - 2 * ave_state( x)) * m_AverageBelowCutoffSlopes( x) if ave_state( x) >= 0.5
        // else, f( x) = cutoff + ( 1 - 2 * ave_state( x)) * m_AverageAbvoveCutoffSlopes( x) if ave_state( x) < 0.5

        storage::Vector< math::RunningAverage< float> > activities( results_size);

        for( size_t result_number( 0); result_number < results_size; ++result_number)
        {
          const float ave_state( float( state_counts( result_number)) / float( number_results));
          float value( m_ActivityCutoff);
          if( ave_state >= 0.5)
          {
            value += ( 1.0 - 2.0 * ave_state) * m_AverageBelowCutoffSlopes( result_number);
          }
          else
          {
            value += ( 1.0 - 2.0 * ave_state) * m_AverageAboveCutoffSlopes( result_number);
          }
          activities( result_number).AddWeightedObservation( value, double( number_results));
        }

        // set the average activity for the node
        NODE->SetActivity( activities);
      }

      // determine the best partition of the data
      DtreeBinaryPartition node_partition( m_DataPartitioner->operator()( *DATA));

      // only add the partition if it was fully defined
      if
      (
        util::IsDefined( node_partition.GetFeatureIndex())
        && util::IsDefined( node_partition.GetSplitValue())
        && util::IsDefined( node_partition.GetSplitRating())
        && node_partition.GetInitialNumIncorrect() >= m_MinIncorrect
      )
      {
        // store the partition's rating for easy access
        const float partiton_rating( node_partition.GetData( m_PartitionScoreType));

        storage::Triplet
        <
          util::ShPtr< DecisionTree>,                             // the tree to be expanded
          util::ShPtr< storage::Vector< FeatureResultAndState> >, // references of feature results to expand
          DtreeBinaryPartition                                    // best binary partition for this dataset
        > new_node_to_expand( NODE, DATA, node_partition);

        // determine where to insert this node in the nodes to be partitioned by sorting it into the list based
        // on the split rating.  Always put the node with the highest rating last, so we can just always pop the last
        // member of the list off in operator()
        storage::List
        <
          storage::Triplet
          <
            util::ShPtr< DecisionTree>,                             // the tree to be expanded
            util::ShPtr< storage::Vector< FeatureResultAndState> >, // references of feature results to expand
            DtreeBinaryPartition                                    // best binary partition for this dataset
          >
        >::iterator itr_list( m_UnexpandedSubTrees.Begin()), itr_list_end( m_UnexpandedSubTrees.End());

        // determine a secondary partition rating as a tie-breaker
        const DtreeBinaryPartition::Data tie_breaker_score_type
        (
          m_PartitionScoreType == DtreeBinaryPartition::e_InitialNumIncorrect
          || m_PartitionScoreType == DtreeBinaryPartition::e_InitialIncorrectPlusFinalCorrect
          ? DtreeBinaryPartition::e_SplitRating
          : DtreeBinaryPartition::e_InitialNumIncorrect
        );

        const double tie_breaker_score( node_partition.GetData( tie_breaker_score_type));

        while
        (
          itr_list != itr_list_end
          && partiton_rating >= itr_list->Third().GetData( m_PartitionScoreType)
        )
        {
          if( partiton_rating == itr_list->Third().GetData( m_PartitionScoreType))
          {
            if( tie_breaker_score <= itr_list->Third().GetData( tie_breaker_score_type))
            {
              break;
            }
          }
          ++itr_list;
        }

        // insert the new node to expand
        m_UnexpandedSubTrees.InsertElement( itr_list, new_node_to_expand);
      }

    }

  } // namespace model
} // namespace bcl
