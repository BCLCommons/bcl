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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_PARTIAL_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_PARTIAL_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_interface.h"
#include "math/bcl_math_range_set.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionPartial
    //! @brief Allows the core objective function to operate on a selected subset of each result, also allows weighting
    //!
    //! @see @link example_model_objective_function_partial.cpp @endlink
    //! @author mendenjl
    //! @date Jul 18, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionPartial :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! objective function that will be called and weighted
      util::Implementation< ObjectiveFunctionInterface> m_Objective;

      //! ranges of output data over which the objective function will be applied
      math::RangeSet< size_t> m_Ranges;

      //! weight for this partial objective function
      float m_Weight;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ObjectiveFunctionPartial
      ObjectiveFunctionPartial *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_Name( "PartialObjective");
        return s_Name;
      }

      //! @brief determine what sign of the derivative of this objective function indicates improvement
      //! @return the sign of the derivative of this objective function indicates improvement
      opti::ImprovementType GetImprovementType() const;

      //! @brief get the overall goal of the objective function
      //! @return the goal of the objective function
      Goal GetGoalType() const;

      //! @brief get the threshold, for classification type objectives
      //! @return the threshold, for classification type objectives
      float GetThreshold() const
      {
        return m_Objective->GetThreshold();
      }

      //! @brief set the threshold
      //! @param THRESHOLD threshold that divides output classes
      void SetThreshold( const float &THRESHOLD)
      {
        m_Objective->SetThreshold( THRESHOLD);
      }

      //! @brief get the parity, for rank classification type objectives
      //! @return the parity, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
      //!       in prediction of values higher than the threshold, false means below threshold
      bool GetRankingParity() const
      {
        return m_Objective->GetRankingParity();
      }

      //! @brief get the desired hit rate (e.g. fraction of positive predictions)
      //! @return the desired hit rate, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
      //!       of values will be considered
      float GetDesiredHitRate() const
      {
        return m_Objective->GetDesiredHitRate();
      }

      //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
      //! @param DATA monitoring dataset results, non-scaled
      //! @param IDS ids; can be used by the objective function
      void SetData
      (
        const FeatureDataSet< float> &DATA,
        const FeatureDataSet< char> &IDS = FeatureDataSet< char>()
      );

    ////////////////
    // operations //
    ////////////////

      //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
      //!        case indicates whether the prediction was true (upper) or false (lower)
      //! @param EXPERIMENTAL feature dataset with experimental values
      //! @param PREDICTED feature dataset with predicted values
      //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
      linal::Matrix< char> GetFeaturePredictionClassifications
      (
        const FeatureDataSetInterface< float> &EXPERIMENTAL,
        const FeatureDataSetInterface< float> &PREDICTED
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
      //! @param EXPERIMENTAL feature dataset with experimental values
      //! @param PREDICTED feature dataset with predicted values
      //! @return objective function value based on the given data of experimental and predicted values
      float operator()
      (
        const FeatureDataSetInterface< float> &EXPERIMENTAL,
        const FeatureDataSetInterface< float> &PREDICTED
      ) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ObjectiveFunctionPartial

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_PARTIAL_H_
