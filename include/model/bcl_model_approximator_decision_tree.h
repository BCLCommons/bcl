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

#ifndef BCL_MODEL_APPROXIMATOR_DECISION_TREE_H_
#define BCL_MODEL_APPROXIMATOR_DECISION_TREE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_decision_tree.h"
#include "bcl_model_dtree_data_partition_function_interface.h"
#include "bcl_model_feature_result_and_state.h"
#include "bcl_model_objective_function_wrapper.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorDecisionTree
    //! @brief is derived from ApproximatorBase and is part of the approximator framework.
    //!        It will train a decision tree accoding to a given objective function and training data.
    //!        DecisionTreeIterate will use a DtreeDataPartitionFunctionInterface to define which algorithm is used
    //!        for partitioning the training data in the learning process. No rescaling of training data is required.
    //!
    //! @see @link example_model_approximator_decision_tree.cpp @endlink
    //! @author lemmonwa, mendenjl
    //! @date Jul 03, 2009; July 01, 2010; Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorDecisionTree :
      public ApproximatorBase
    {

    private:

    //////////
    // data //
    //////////

      //! ShPtr to NodeSplitApproximater which does the splitting of nodes to be added to the decision tree
      util::Implementation< DtreeDataPartitionFunctionInterface> m_DataPartitioner;

      //! Value used to determine active or inactive
      float m_ActivityCutoff;

      //! Method by which to decide which node to split next (higher values are better)
      DtreeBinaryPartition::DataEnum m_PartitionScoreType;

      //! holds trees to be visited in future iterations
      storage::List
      <
        storage::Triplet
        <
          util::ShPtr< DecisionTree>,                             // the tree to be expanded
          util::ShPtr< storage::Vector< FeatureResultAndState> >, // references of feature results to expand
          DtreeBinaryPartition                                    // best binary partition for this dataset
        >
      > m_UnexpandedSubTrees;

      //! pointer to the current decision tree for evaluation purposes
      util::ShPtr< DecisionTree> m_DecisionTree;

      //! result states (all size_t's are actually binary)
      linal::Matrix< size_t> m_ResultStates;

      //! Minimum # of items incorrect to decide to split a node
      size_t m_MinIncorrect;

      //! Average of all values above the cutoff
      linal::Vector< float> m_AverageAboveCutoffSlopes;

      //! Average of all values below the cutoff
      linal::Vector< float> m_AverageBelowCutoffSlopes;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ApproximatorDecisionTree();

      //! @brief constructor from all necessary parameters
      //! @param DATA_PARTITIONER ShPtr to DtreeDataPartitionFunctionInterface
      //! @param OBJECTIVE_FUNCTION ShPtr to objective function
      //! @param ACTIVITY_CUTOFF the cutoff for binary classification
      //! @param TRAINING_DATA the dataset used to build the tree
      ApproximatorDecisionTree
      (
        const util::Implementation< DtreeDataPartitionFunctionInterface> &DATA_PARTITIONER,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
        const float ACTIVITY_CUTOFF,
        const util::ShPtr< descriptor::Dataset> &TRAINING_DATA
      );

      //! @brief clone function
      ApproximatorDecisionTree *Clone() const
      {
        return new ApproximatorDecisionTree( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns activity cutoff
      //! @return the activity cutoff
      float GetActivityCutoff() const
      {
        return m_ActivityCutoff;
      }

      //! @brief get training data set for a specific iterate in approximater framework
      //! @return dataset interface with training data
      const util::ShPtr< descriptor::Dataset> &GetTrainingData() const
      {
        return m_TrainingData;
      }

      //! @brief set training data set for a specific iterate in approximater framework
      //! @param DATA training data set
      void SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA);

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns ShPtr to NodeSplitter used
      //! @return ShPtr to NodeSplitter used
      const util::Implementation< DtreeDataPartitionFunctionInterface> &GetDataPartitionFunction() const
      {
        return m_DataPartitioner;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief Prune a decision tree model based on size & cutoff (if applicable)
      //! @param MODEL the decision tree to prune
      void Prune( DecisionTree &MODEL) const;

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const;

      //! @brief constructs children nodes, adds them to current node, and queues children for visiting
      //! @param CURRENT_NODE a ShPtr to the node to expand
      //! @param DATA data from the parent node
      //! @param BINARY_PARTITION partition to use
      void Expand
      (
        util::ShPtr< DecisionTree> &CURRENT_NODE,
        const util::ShPtr< storage::Vector< FeatureResultAndState> > &DATA,
        const DtreeBinaryPartition &BINARY_PARTITION
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief converts an original training data set into an already classified version based on a given cutoff
      //! @param FEATURES original features of interest
      //! @param RESULTS corresponding results
      //! @return pre-processed data set
      util::ShPtr< storage::Vector< FeatureResultAndState> > CreateDataSetReferences
      (
        const linal::MatrixConstInterface< float> &FEATURES,
        const linal::MatrixConstInterface< float> &RESULTS
      );

    private:

      //! @brief add a node, determines average activity, etc., sorts it into the list of nodes to expand
      //! @param NODE the node to add
      //! @param DATA the data that maps to that node
      void AddExpandableNode
      (
        util::ShPtr< DecisionTree> NODE,
        const util::ShPtr< storage::Vector< FeatureResultAndState> > &DATA
      );

    }; // DecisionTreeIterate

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_DECISION_TREE_H_
