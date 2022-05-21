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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_ENRICHMENT_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_ENRICHMENT_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionEnrichment
    //! @brief The objective function for monitoring a model::Interface
    //! @details This class determines the prediction error of the model::Interface model on the monitoring data set.
    //! The model::Interface improves as long as the error goes down.
    //!
    //! @see @link example_model_objective_function_enrichment.cpp @endlink
    //! @author butkiem1
    //! @date Jan 8, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionEnrichment :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief cutoff to determining the area under the curve
      float m_Cutoff;

      //! @brief oversampling factor in monitoring dataset
      float m_OverSamplingFactor;

      //! @brief enrichment cutoff
      float m_FalsePositiveRateCutoff;

      //! flag for sorting list of predicted and experimental values in roc curve
      bool m_PositivesAboveThreshold;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionEnrichment();

      //! @brief constructor from monitor data
      ObjectiveFunctionEnrichment
      (
        const float CUTOFF,
        const float OVERSAMPLING_FACTOR,
        const bool POSITIVES_ABOVE_THRESHOLD
      );

      //! @brief Clone function
      //! @return pointer to new ObjectiveFunctionEnrichment
      ObjectiveFunctionEnrichment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_Name( "Enrichment");
        return s_Name;
      }

      //! @brief determine what sign of the derivative of this objective function indicates improvement
      //! @return the sign of the derivative of this objective function indicates improvement
      opti::ImprovementType GetImprovementType() const
      {
        return opti::e_LargerEqualIsBetter;
      }

      //! @brief get the overall goal of the objective function
      //! @return the goal of the objective function
      Goal GetGoalType() const
      {
        return e_RankClassification;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief set cutoff to for monitoring dataset evaluation
      //! @param CUTOFF cutoff used in ROC curve determination
      void SetCutoff( const float &CUTOFF)
      {
        m_Cutoff = CUTOFF;
      }

      //! @brief get the threshold, for classification type objectives
      //! @return the threshold, for classification type objectives
      float GetThreshold() const
      {
        return m_Cutoff;
      }

      //! @brief set the threshold
      //! @param THRESHOLD threshold that divides output classes
      void SetThreshold( const float &THRESHOLD)
      {
        m_Cutoff = THRESHOLD;
      }

      //! @brief get the parity, for rank classification type objectives
      //! @return the parity, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
      //!       in prediction of values higher than the threshold, false means below threshold
      bool GetRankingParity() const
      {
        return m_PositivesAboveThreshold;
      }

      //! @brief set oversampling factor for monitoring dataset evaluation
      //! @param FACTOR oversampling factor used in ROC curve determination
      void SetOverSamplingFactor( const float &FACTOR)
      {
        BCL_Assert
        (
          util::IsDefined( FACTOR) && FACTOR > 0.0,
          GetClassIdentifier() + " needs a defined oversampling factor > 0"
        );
        m_OverSamplingFactor = FACTOR;
      }

      //! @brief get oversampling factor for model::Interface
      //! @return oversampling factor used to construct monitoring dataset
      const float &GetOverSamplingFactor() const
      {
        return m_OverSamplingFactor;
      }

      //! @brief get the desired hit rate (e.g. fraction of positive predictions)
      //! @return the desired hit rate, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
      //!       of values will be considered
      float GetDesiredHitRate() const
      {
        return m_FalsePositiveRateCutoff;
      }

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

    }; // class ObjectiveFunctionEnrichment

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_ENRICHMENT_H_
