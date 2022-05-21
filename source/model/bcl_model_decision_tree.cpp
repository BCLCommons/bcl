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
#include "model/bcl_model_decision_tree.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_feature_data_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DecisionTree::s_Instance
    (
      GetObjectInstances().AddInstance( new DecisionTree())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor for a decision tree
    DecisionTree::DecisionTree() :
      m_ModalActivity(),
      m_Count( 0),
      m_SplitIndex( util::GetUndefined< size_t>()),
      m_SplitValue( util::GetUndefined< float>())
    {
    }

    //! @brief constructor with parameters
    //! @param ACTIVITY average activity of this node
    DecisionTree::DecisionTree( const linal::Vector< float> &ACTIVITY) :
      m_ModalActivity( ACTIVITY),
      m_Count( 1),
      m_SplitIndex( util::GetUndefined< size_t>()),
      m_SplitValue( util::GetUndefined< float>())
    {
    }

    //! @brief clone this and any underlying objects
    DecisionTree *DecisionTree::DeepClone() const
    {
      DecisionTree *deep_clone( Clone());

      if( !IsLeaf())
      {
        deep_clone->ConfigureBranch
        (
          m_SplitIndex,
          m_SplitValue,
          util::ShPtr< DecisionTree>( m_LeftChild->DeepClone()),
          util::ShPtr< DecisionTree>( m_RightChild->DeepClone())
        );
      }
      return deep_clone;
    }

    //! @brief remove the leaves of this decision tree if both leaves have ave activities on the same side of cutoff
    //!        or if either leaf is has < SIZE_CUTOFF members
    //! @param ACTIVITY_CUTOFF the cutoff to use; use undefined if only filtering by leaf size
    //! @param SIZE_CUTOFF the size cutoff to use, use 0 if only filtering by activity
    //! @note if both criteria are specified, then either criterion can trigger the pruning
    void DecisionTree::Prune( const double &ACTIVITY_CUTOFF, const size_t &SIZE_CUTOFF)
    {
      if( IsLeaf())
      {
        // leaves cannot prune themselves, so just return
        return;
      }

      // prune child nodes first before testing whether they are leaves
      m_LeftChild->Prune( ACTIVITY_CUTOFF, SIZE_CUTOFF);
      m_RightChild->Prune( ACTIVITY_CUTOFF, SIZE_CUTOFF);

      // if either child is not a leaf, there is nothing more to do
      if( !m_LeftChild->IsLeaf() || !m_RightChild->IsLeaf())
      {
        return;
      }

      // if either child is too small, prune
      if( m_LeftChild->m_Count < SIZE_CUTOFF || m_RightChild->m_Count < SIZE_CUTOFF)
      {
        // if either of the children have too few nodes, make this node a leaf
        MakeLeaf();
      }
      else if( util::IsDefined( ACTIVITY_CUTOFF))
      {
        const linal::Vector< float> &activities_left( m_LeftChild->GetActivity());
        const linal::Vector< float> &activities_right( m_RightChild->GetActivity());
        for( size_t i( 0), n_outputs( GetNumberOutputs()); i < n_outputs; ++i)
        {
          // check whether the features are on opposite sides of the cutoff
          if( bool( activities_left( i) <= ACTIVITY_CUTOFF) != bool( activities_right( i) <= ACTIVITY_CUTOFF))
          {
            return;
          }
        }

        // all activities were on the same side of the cutoff, so both child nodes are equivalent to this node
        // and can be removed
        MakeLeaf();
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set the average activity of this node of the decision tree
    //! @param AVERAGES averages of each result
    void DecisionTree::SetActivity( const storage::Vector< math::RunningAverage< float> > &AVERAGES)
    {
      // set modal activity to the correct size
      m_ModalActivity = linal::Vector< float>( AVERAGES.GetSize(), util::GetUndefined< float>());

      if( AVERAGES.IsEmpty())
      {
        // empty dataset, just set m_Count to 0 and return
        m_Count = 0;
      }
      else
      {
        // set m_Count
        m_Count = AVERAGES.FirstElement().GetWeight();
      }

      // for each result, set the modal activity to the value from AVERAGES_ABOVE / AVERAGES_BELOW if there were more (Above/Below) the cutoff
      for
      (
        size_t result_number( 0), number_results( AVERAGES.GetSize());
        result_number < number_results;
        ++result_number
      )
      {
        m_ModalActivity( result_number) = AVERAGES( result_number).GetAverage();
      }
    }

    //! @brief get the size of the tree rooted at this node
    //! @return the size of the tree rooted at this node
    size_t DecisionTree::GetTreeSize() const
    {
      return 1 + ( IsLeaf() ? 0 : m_LeftChild->GetTreeSize() + m_RightChild->GetTreeSize());
    }

    //! @brief configure this node to be a branch node with the given parameters
    //! @param SPLIT_INDEX specifies which index of feature to split based off of
    //! @param SPLIT_VALUE specifies which value of selected feature; values <= this value will be sent to left child, > will be sent right
    //! @param DECISION_NODE_PAIR shptr to pair with node and their decision function
    //! @param LEFT_CHILD, RIGHT_CHILD the children of this branch
    void DecisionTree::ConfigureBranch
    (
      const size_t &SPLIT_INDEX,
      const float &SPLIT_VALUE,
      const util::ShPtr< DecisionTree> &LEFT_CHILD,
      const util::ShPtr< DecisionTree> &RIGHT_CHILD
    )
    {
      m_SplitIndex = SPLIT_INDEX;
      m_SplitValue = SPLIT_VALUE;
      m_LeftChild = LEFT_CHILD;
      m_RightChild = RIGHT_CHILD;
      BCL_Assert
      (
        m_LeftChild.IsDefined() == m_RightChild.IsDefined(),
        "Children of a decision tree node must be either both defined or undefined"
      );

      BCL_Assert
      (
        m_LeftChild != m_RightChild,
        "Children of a decision tree node must be either both defined or undefined"
      );
    }

    //! @brief get a vector of decision trees, ordered by depth in the sub-tree
    //! @return the subtree below this node
    storage::Vector< util::ShPtrVector< DecisionTree> > DecisionTree::GetSubtreeLevels() const
    {
      storage::Vector< util::ShPtrVector< DecisionTree> > forest;
      if( !IsLeaf())
      {
        // get the left and right forests too
        storage::Vector< util::ShPtrVector< DecisionTree> > left_forest( m_LeftChild->GetSubtreeLevels());
        storage::Vector< util::ShPtrVector< DecisionTree> > right_forest( m_RightChild->GetSubtreeLevels());

        // allocate enough memory to hold the new nodes
        forest.Resize( std::max( left_forest.GetSize(), right_forest.GetSize()) + 1);

        // add this node's children to the forest
        forest( 0).PushBack( m_LeftChild);
        forest( 0).PushBack( m_RightChild);

        // concatenate the resulting forests
        size_t depth( 1);
        for
        (
          storage::Vector< util::ShPtrVector< DecisionTree> >::const_iterator
          itr_forest( left_forest.Begin()), itr_forest_end( left_forest.End());
          itr_forest != itr_forest_end;
          ++itr_forest, ++depth
        )
        {
          forest( depth).Append( *itr_forest);
        }
        depth = 1;
        for
        (
          storage::Vector< util::ShPtrVector< DecisionTree> >::const_iterator
          itr_forest( right_forest.Begin()), itr_forest_end( right_forest.End());
          itr_forest != itr_forest_end;
          ++itr_forest, ++depth
        )
        {
          forest( depth).Append( *itr_forest);
        }
      }
      return forest;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURES not rescaled features
    //! @return predcited result vector using a model
    FeatureDataSet< float> DecisionTree::PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURES) const
    {
      linal::Matrix< float> results( FEATURES.GetNumberFeatures(), m_ModalActivity.GetSize());
      for
      (
        size_t feature_number( 0), number_features( FEATURES.GetNumberFeatures());
        feature_number < number_features;
        ++feature_number
      )
      {
        PredictFeature( FEATURES( feature_number), results[ feature_number]);
      }
      return FeatureDataSet< float>( results);
    }

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void DecisionTree::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      // no rescaling is ever performed by decision trees, so just remove any rescaling on the dataset
      FEATURE.DeScale();
    }

    //! @brief make a branch node into a leaf node
    void DecisionTree::MakeLeaf()
    {
      m_SplitIndex = util::GetUndefined< size_t>();
      m_SplitValue = util::GetUndefined< float>();
      m_LeftChild  = util::ShPtr< DecisionTree>();
      m_RightChild = util::ShPtr< DecisionTree>();
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURES normalized or rescaled features
    //! @return predcited results using a model
    FeatureDataSet< float> DecisionTree::operator()( const FeatureDataSetInterface< float> &FEATURES) const
    {
      return PredictWithoutRescaling( FEATURES);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DecisionTree::Read( std::istream &ISTREAM)
    {
      // read in members
      io::Serialize::Read( m_SplitIndex, ISTREAM);
      io::Serialize::Read( m_SplitValue, ISTREAM);
      io::Serialize::Read( m_Count, ISTREAM);
      io::Serialize::Read( m_ModalActivity, ISTREAM);

      // if this was a branch node, read in the children
      if( util::IsDefined( m_SplitIndex))
      {
        io::Serialize::Read( m_LeftChild, ISTREAM);
        io::Serialize::Read( m_RightChild, ISTREAM);
      }
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT indentation
    //! @return outputstream which was written to
    std::ostream &DecisionTree::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_SplitIndex, OSTREAM, INDENT) << ' ';
      io::Serialize::Write( m_SplitValue, OSTREAM) << ' ';
      io::Serialize::Write( m_Count, OSTREAM) << '\n';

      io::Serialize::Write( m_ModalActivity, OSTREAM, INDENT);

      // write the children if this was a branch node
      if( util::IsDefined( m_SplitIndex))
      {
        OSTREAM << '\n';
        io::Serialize::Write( m_LeftChild,  OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_RightChild, OSTREAM, INDENT);
      }

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DecisionTree::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "see http://en.wikipedia.org/wiki/Decision_tree");
      return parameters;
    }

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @param RESULT where to store the result
    void DecisionTree::PredictFeature
    (
      const FeatureReference< float> &FEATURE,
      float *RESULT
    ) const
    {
      // for branch nodes, pass the decision on to the correct branch
      if( !IsLeaf())
      {
        if( FEATURE( m_SplitIndex) <= m_SplitValue)
        {
          m_LeftChild->PredictFeature( FEATURE, RESULT);
        }
        else
        {
          m_RightChild->PredictFeature( FEATURE, RESULT);
        }
      }
      else
      {
        // for leaf nodes, just copy the result
        std::copy( m_ModalActivity.Begin(), m_ModalActivity.End(), RESULT);
      }
    }

  } // namespace model
} // namespace bcl
