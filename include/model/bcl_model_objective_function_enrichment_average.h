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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_ENRICHMENT_AVERAGE_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_ENRICHMENT_AVERAGE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

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
    //! @class ObjectiveFunctionEnrichmentAverage
    //! @brief An objective function for monitoring a model::Interface using an enrichment average
    //! @details This class determines the average enrichment of the model::Interface model on the monitoring data set.
    //! @see @link example_model_objective_function_enrichment_average.cpp @endlink
    //!
    //! @author butkiem1
    //! @date Oct 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionEnrichmentAverage :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! @brief cutoff to determining the area under the curve
      float m_Cutoff;

      //! @brief max cutoff of enrichment sampling
      float m_EnrichmentCutOffMax;

      //! @brief stepsize for evaluating enrichment at cutoff
      float m_EnrichmentCutOffStepSize;

      //! flag for sorting list of predicted and experimental values in roc curve
      bool m_PositivesAboveThreshold;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////
    // data //
    //////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionEnrichmentAverage() :
        m_Cutoff( 0.0),
        m_EnrichmentCutOffMax( 0.0),
        m_EnrichmentCutOffStepSize( 0.0),
        m_PositivesAboveThreshold( true)
      {
      }

      //! @brief constructor from monitor data
      //! @param CUTOFF cutoff to determine roc curve that is used for enrichment calculation
      //! @param ENRICHMENT_CUTOFF_MAX max cutoff value for that enrichment will be sampled
      //! @param ENRICHMENT_CUTOFF_STEPSIZE step size of cutoffs for sampling enrichment values
      //! @param DATA dataset on which objective function is evaluated
      ObjectiveFunctionEnrichmentAverage
      (
        const float CUTOFF,
        const float ENRICHMENT_CUTOFF_MAX,
        const float ENRICHMENT_CUTOFF_STEPSIZE,
        const bool POSITIVES_ABOVE_THRESHOLD
      ) :
        m_Cutoff( CUTOFF),
        m_EnrichmentCutOffMax( ENRICHMENT_CUTOFF_MAX),
        m_EnrichmentCutOffStepSize( ENRICHMENT_CUTOFF_STEPSIZE),
        m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD)
      {
      }

      //! @brief Clone function
      //! @return pointer to new ObjectiveFunctionEnrichmentAverage
      ObjectiveFunctionEnrichmentAverage *Clone() const;

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
        static const std::string s_Name( "EnrichmentAverage");
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

      //! @brief get the threshold, for classification type objectives
      //! @return the threshold, for classification type objectives
      float GetThreshold() const
      {
        return m_Cutoff;
      }

      //! @brief get the parity, for rank classification type objectives
      //! @return the parity, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
      //!       in prediction of values higher than the threshold, false means below threshold
      bool GetRankingParity() const
      {
        return m_PositivesAboveThreshold;
      }

      //! @brief get the desired hit rate (e.g. fraction of positive predictions)
      //! @return the desired hit rate, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
      //!       of values will be considered
      float GetDesiredHitRate() const
      {
        return m_EnrichmentCutOffMax;
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

      //! @brief set the threshold
      //! @param THRESHOLD threshold that divides output classes
      void SetThreshold( const float &THRESHOLD)
      {
        m_Cutoff = THRESHOLD;
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

    }; // class ObjectiveFunctionEnrichmentAverage

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_ENRICHMENT_AVERAGE_H_
