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
#include "model/bcl_model_objective_function_binary_operation.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionBinaryOperation::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionBinaryOperation()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new Evaluator
    ObjectiveFunctionBinaryOperation *ObjectiveFunctionBinaryOperation::Clone() const
    {
      return new ObjectiveFunctionBinaryOperation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionBinaryOperation::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      m_ObjectiveLeft->SetData( DATA, IDS);
      m_ObjectiveRight->SetData( DATA, IDS);
    }

    //! @brief determine what sign of the derivative of this objective function indicates improvement
    //! @return the sign of the derivative of this objective function indicates improvement
    opti::ImprovementType ObjectiveFunctionBinaryOperation::GetImprovementType() const
    {
      // ensure that all the members are defined
      BCL_Assert( m_ObjectiveLeft.IsDefined(), "LHS was undefined");
      BCL_Assert( m_ObjectiveRight.IsDefined(), "RHS was undefined");
      BCL_Assert( m_Op.IsDefined(), "Operation was undefined");

      // get the improvement type on each side
      const opti::ImprovementType lhs_type( m_ObjectiveLeft->GetImprovementType());
      const opti::ImprovementType rhs_type( m_ObjectiveRight->GetImprovementType());

      if( m_Op.GetAlias() == "+")
      {
        BCL_Assert
        (
          lhs_type == rhs_type,
          "Improvement type is ill defined: adding objectives with opposite signs of improvement"
        );
        return lhs_type;
      }
      else if( m_Op.GetAlias() == "-")
      {
        BCL_Assert
        (
          lhs_type != rhs_type,
          "Improvement type is ill defined: subtracting objectives with same signs of improvement"
        );
        return lhs_type;
      }
      else if( m_Op.GetAlias() == "*")
      {
        BCL_Assert
        (
          lhs_type == rhs_type,
          "Improvement type is ill defined: multiplying objectives with opposite signs of improvement"
        );
        return lhs_type;
      }
      else if( m_Op.GetAlias() == "/")
      {
        BCL_Assert
        (
          lhs_type != rhs_type,
          "Improvement type is ill defined: dividing objectives with same signs of improvement"
        );
        return lhs_type;
      }
      BCL_Assert
      (
        m_ObjectiveLeft->GetRankingParity() == m_ObjectiveRight->GetRankingParity(),
        "Objectives with opposite ranking parities cannot be combined!"
      );
      BCL_Exit( "Undefined improvement type for op: " + m_Op.GetAlias(), -1);
      return lhs_type;
    }

    //! @brief get the overall goal of the objective function
    //! @return the goal of the objective function
    ObjectiveFunctionInterface::Goal ObjectiveFunctionBinaryOperation::GetGoalType() const
    {
      BCL_Assert( m_ObjectiveLeft.IsDefined(), "LHS was undefined");
      BCL_Assert( m_ObjectiveRight.IsDefined(), "RHS was undefined");
      if( m_ObjectiveLeft->GetGoalType() == m_ObjectiveRight->GetGoalType())
      {
        return m_ObjectiveLeft->GetGoalType();
      }
      return e_Other;
    }

    //! @brief get the threshold, for classification type objectives
    //! @return the threshold, for classification type objectives
    float ObjectiveFunctionBinaryOperation::GetThreshold() const
    {
      if( m_ObjectiveLeft->GetThreshold() == m_ObjectiveRight->GetThreshold())
      {
        return m_ObjectiveLeft->GetThreshold();
      }
      return util::GetUndefined< float>();
    }

    //! @brief set the threshold
    //! @param THRESHOLD threshold that divides output classes
    void ObjectiveFunctionBinaryOperation::SetThreshold( const float &THRESHOLD)
    {
      m_ObjectiveLeft->SetThreshold( THRESHOLD);
      m_ObjectiveRight->SetThreshold( THRESHOLD);
    }

    //! @brief get the parity, for rank classification type objectives
    //! @return the parity, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
    //!       in prediction of values higher than the threshold, false means below threshold
    bool ObjectiveFunctionBinaryOperation::GetRankingParity() const
    {
      BCL_Assert
      (
        m_ObjectiveLeft->GetRankingParity() == m_ObjectiveRight->GetRankingParity(),
        "Objectives with opposite ranking parities cannot be combined!"
      );
      return m_ObjectiveLeft->GetRankingParity();
    }

    //! @brief get the desired hit rate (e.g. fraction of positive predictions)
    //! @return the desired hit rate, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
    //!       of values will be considered
    float ObjectiveFunctionBinaryOperation::GetDesiredHitRate() const
    {
      const float desired_left( m_ObjectiveLeft->GetDesiredHitRate());
      const float desired_right( m_ObjectiveRight->GetDesiredHitRate());
      if( !util::IsDefined( desired_left))
      {
        return desired_right;
      }
      else if( !util::IsDefined( desired_right))
      {
        return desired_left;
      }
      // both sides defined, take the average
      return 0.5 * ( desired_right + desired_left);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    linal::Matrix< char> ObjectiveFunctionBinaryOperation::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // number of data points in dataset
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());

      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == PREDICTED.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );
      BCL_Assert
      (
        EXPERIMENTAL.GetFeatureSize() == PREDICTED.GetFeatureSize(),
        "Experimental and predicted values have different result sizes"
      );

      linal::Matrix< char> prediction_classification
      (
        m_ObjectiveLeft->GetFeaturePredictionClassifications( EXPERIMENTAL, PREDICTED)
      );
      const linal::Matrix< char> prediction_classification_rhs
      (
        m_ObjectiveRight->GetFeaturePredictionClassifications( EXPERIMENTAL, PREDICTED)
      );

      // combine the matrices using the max operation, this prefers false predictions over true, and any prediction over none
      for( size_t counter( 0); counter < data_set_size; ++counter)
      {
        // iterate over each result
        for( size_t result_number( 0); result_number < result_size; ++result_number)
        {
          prediction_classification( counter, result_number) =
            std::max
            (
              prediction_classification( counter, result_number),
              prediction_classification_rhs( counter, result_number)
            );
        }
      }

      // return the classification
      return prediction_classification;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionBinaryOperation::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // compute the left hand side
      float value_left( m_ObjectiveLeft->operator ()( EXPERIMENTAL, PREDICTED));
      const float value_right( m_ObjectiveRight->operator ()( EXPERIMENTAL, PREDICTED));

      // carry out the operation
      m_Op->operator ()( value_left, value_right);

      // return the resulting value, now stored in value_left
      return value_left;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionBinaryOperation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Performs an arithmetic operation on the outputs of two objective functions"
      );

      parameters.AddInitializer
      (
        "lhs",
        "left hand side of the operation",
        io::Serialization::GetAgent( &m_ObjectiveLeft)
      );

      parameters.AddInitializer
      (
        "rhs",
        "right hand side of the operation",
        io::Serialization::GetAgent( &m_ObjectiveRight)
      );

      parameters.AddInitializer
      (
        "op",
        "operation to carry out on the lhs and rhs",
        io::Serialization::GetAgent( &m_Op)
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
