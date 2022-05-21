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

#ifndef BCL_MODEL_DECISION_TREE_H_
#define BCL_MODEL_DECISION_TREE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_interface.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DecisionTree
    //! @brief is a decision tree implementation derived from model::Interface.
    //!        This class contains sub trees with their corresponding decision functions.
    //!        Leaf trees which do not contain subtrees are considered leaves and give back the ratio of activity
    //!        for the particular dataset partition the subtree operates on.
    //!
    //! @see @link example_model_decision_tree.cpp @endlink
    //! @author teixeipl, lemmonwa, butkiem1, mendenjl
    //! @date 7/04/2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DecisionTree :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! Most common activities; each result value in this array is the average of values above/below cutoff, whichever
      //! had more results assigned to it (in the case of equal #s of results, an average is made)
      linal::Vector< float> m_ModalActivity;
      size_t                m_Count; //!< number of training values assigned to this node

      // the following data members are only defined for branch nodes; for other nodes, they are undefined
      size_t m_SplitIndex; //!< for branch nodes, specifies which index of feature to split based off of, for leaf nodes, undefined
      float  m_SplitValue; //!< for branch nodes, specifies which value of selected feature; values <= this value will be sent to left child, > will be sent right

      util::ShPtr< DecisionTree> m_LeftChild; //!< child to send features to if feature( m_SplitIndex) <= m_Value
      util::ShPtr< DecisionTree> m_RightChild; //!< child to send features to if feature( m_SplitIndex) > m_Value

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor for a decision tree
      DecisionTree();

      //! @brief constructor with parameters
      //! @param AVERAGE_ACTIVITY average activity of this node
      explicit DecisionTree( const linal::Vector< float> &AVERAGE_ACTIVITY);

      //! @brief clone method
      DecisionTree *Clone() const
      {
        return new DecisionTree( *this);
      }

      //! @brief clone this and any underlying objects
      DecisionTree *DeepClone() const;

      //! @brief remove the leaves of this decision tree if both leaves have ave activities on the same side of cutoff
      //!        or if either leaf is has < SIZE_CUTOFF members
      //! @param ACTIVITY_CUTOFF the cutoff to use; use undefined if only filtering by leaf size
      //! @param SIZE_CUTOFF the size cutoff to use, use 0 if only filtering by activity
      void Prune( const double &ACTIVITY_CUTOFF, const size_t &SIZE_CUTOFF = 0);

    /////////////////
    // data access //
    /////////////////

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const
      {
        return m_ModalActivity.GetSize();
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_Name( "DecisionTree");
        return s_Name;
      }

      //! @brief get the split index of feature to split based off of (only defined for branch nodes)
      //! @return the split index of feature to split based off of (only defined for branch nodes)
      const size_t &GetSplitIndex() const
      {
        return m_SplitIndex;
      }

      //! @brief get the split value of feature to split based off of (only defined for branch nodes)
      //! @return the split value of feature to split based off of (only defined for branch nodes)
      const float &GetSplitValue() const
      {
        return m_SplitValue;
      }

      //! @brief get the average activity of training values assigned to this node
      //! @return the average activity of training values assigned to this node
      const linal::Vector< float> &GetActivity() const
      {
        return m_ModalActivity;
      }

      //! @brief set the average activity of this node of the decision tree
      //! @param AVERAGES averages of each result
      void SetActivity( const storage::Vector< math::RunningAverage< float> > &AVERAGES);

      //! @brief checks if this decision tree node is a leaf
      //! @return true if this decision tree node is a leaf
      bool IsLeaf() const
      {
        return !m_LeftChild.IsDefined();
      }

      //! @brief get the size of the tree rooted at this node
      //! @return the size of the tree rooted at this node
      size_t GetTreeSize() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief configure this node to be a branch node with the given parameters
      //! @param SPLIT_INDEX specifies which index of feature to split based off of
      //! @param SPLIT_VALUE specifies which value of selected feature; values <= this value will be sent to left child, > will be sent right
      //! @param LEFT_CHILD, RIGHT_CHILD the children of this branch
      void ConfigureBranch
      (
        const size_t &SPLIT_INDEX,
        const float &SPLIT_VALUE,
        const util::ShPtr< DecisionTree> &LEFT_CHILD,
        const util::ShPtr< DecisionTree> &RIGHT_CHILD
      );

      //! @brief get a vector of decision trees, ordered by depth in the sub-tree
      //! @return the subtree below this node
      storage::Vector< util::ShPtrVector< DecisionTree> > GetSubtreeLevels() const;

      //! @brief make a branch node into a leaf node
      void MakeLeaf();

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predcited result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURES normalized or rescaled features
      //! @return predcited results using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURES) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read decision tree from stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write decision tree to stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @param RESULT where to store the result
      void PredictFeature
      (
        const FeatureReference< float> &FEATURE,
        float *RESULT
      ) const;

    }; // class DecisionTree

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_DECISION_TREE_H_
