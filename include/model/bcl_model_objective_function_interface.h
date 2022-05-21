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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_INTERFACE_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "opti/bcl_opti_improvement_type.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionInterface
    //! @brief to unify all trainable objective functions for model::Interfaces in the approximator framework
    //!
    //! @see @link example_model_objective_function_interface.cpp @endlink
    //! @author mendenjl
    //! @date May 03, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionInterface :
      public math::BinaryFunctionInterfaceSerializable< FeatureDataSetInterface< float>, FeatureDataSetInterface< float>, float>
    {

    public:

    //////////
    // enum //
    //////////

      //! Goal type of the objective function
      enum Goal
      {
        e_Classification,     //!< Objective is based on an absolute threshold
        e_RankClassification, //!< Objective is based on ranking
        e_Regression,         //!< Objective is based on minimizing error
        e_Other,              //!< Objective is based on some other criteria
        s_NumberGoals
      };

      //! @brief Goal as string
      //! @param TYPE the goal
      //! @return the string for the goal
      static const std::string &GetGoalName( const Goal &TYPE);

      //! @brief opti::ImprovementType enum I/O helper
      typedef util::WrapperEnum< Goal, &GetGoalName, s_NumberGoals> GoalEnum;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Interface
      virtual ObjectiveFunctionInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief determine what sign of the derivative of this objective function indicates improvement
      //! @return the sign of the derivative of this objective function indicates improvement
      virtual opti::ImprovementType GetImprovementType() const = 0;

      //! @brief get the overall goal of the objective function
      //! @return the goal of the objective function
      virtual Goal GetGoalType() const = 0;

      //! @brief get the threshold, for classification type objectives
      //! @return the threshold, for classification type objectives
      virtual float GetThreshold() const;

      //! @brief set the threshold
      //! @param THRESHOLD threshold that divides output classes
      virtual void SetThreshold( const float &THRESHOLD)
      {
      }

      //! @brief get the parity, for rank classification type objectives
      //! @return the parity, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
      //!       in prediction of values higher than the threshold, false means below threshold
      virtual bool GetRankingParity() const;

      //! @brief get the desired hit rate (e.g. fraction of positive predictions)
      //! @return the desired hit rate, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
      //!       of values will be considered
      virtual float GetDesiredHitRate() const;

      //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
      //! @param DATA monitoring dataset results, non-scaled
      //! @param IDS ids; can be used by the objective function
      virtual void SetData
      (
        const FeatureDataSet< float> &DATA,
        const FeatureDataSet< char> &IDS = FeatureDataSet< char>()
      )
      {
      }

      //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
      //!        case indicates whether the prediction was true (upper) or false (lower)
      //! @param EXPERIMENTAL feature dataset with experimental values
      //! @param PREDICTED feature dataset with predicted values
      //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
      virtual linal::Matrix< char> GetFeaturePredictionClassifications
      (
        const FeatureDataSetInterface< float> &EXPERIMENTAL,
        const FeatureDataSetInterface< float> &PREDICTED
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
      //! @param EXPERIMENTAL feature dataset with experimental values
      //! @param PREDICTED feature dataset with predicted values
      //! @return objective function value based on the given data of experimental and predicted values
      virtual float operator()
      (
        const FeatureDataSetInterface< float> &EXPERIMENTAL,
        const FeatureDataSetInterface< float> &PREDICTED
      ) const = 0;

    }; // class ObjectiveFunctionInterface
  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_INTERFACE_H_
