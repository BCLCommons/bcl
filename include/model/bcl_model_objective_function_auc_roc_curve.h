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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_AUC_ROC_CURVE_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_AUC_ROC_CURVE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_interface.h"
#include "math/bcl_math_polynomial.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionAucRocCurve
    //! @brief The objective function for monitoring a model::Interface
    //! @details This class determines the prediction error based on area under the calculated roc curve
    //! of the model::Interface model on the monitoring data set.
    //! The model::Interface improves as long as the auc get larger.
    //!
    //! @see @link example_model_objective_function_auc_roc_curve.cpp @endlink
    //! @author butkiem1
    //! @date Sep 19, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionAucRocCurve :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! cutoff to determining the area under the curve
      float m_Cutoff;

      //! weighting function to prioritize regions of auc
      math::Polynomial m_AucWeightingFunction;

      //! flag for sorting list of predicted and experimental values in roc curve
      bool m_PositivesAboveThreshold;

      //! flag for calculating the AUC when x axis is on a log10 scale
      bool m_XAxisLogScale;

      //! fraction predicted cutoff range on x axis of roc curve
      math::Range< double> m_FPRCutoffRange;

      //! Number of data points; needed to establish lower bounds for m_FPRCutoffRange
      size_t m_NumberFeatures;

      //! expected value; what would a perfectly naive model get
      double m_ExpectedResult;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionAucRocCurve() :
        m_Cutoff( 0.0),
        m_AucWeightingFunction( math::Polynomial::MakeFromCoefficients( linal::Vector< float>( 1, 1.0))),
        m_PositivesAboveThreshold( true),
        m_XAxisLogScale( false),
        m_FPRCutoffRange( 0.0, 1.0),
        m_NumberFeatures( 0)
      {
      }

      //! @brief constructor from monitor data with no weighting of area under roc curve
      ObjectiveFunctionAucRocCurve
      (
        const float CUTOFF,
        const bool POSITIVES_ABOVE_THRESHOLD
      ) :
        m_Cutoff( CUTOFF),
        m_AucWeightingFunction( math::Polynomial::MakeFromCoefficients( linal::Vector< float>( 1, 1.0))),
        m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD),
        m_XAxisLogScale( false),
        m_FPRCutoffRange( 0.0, 1.0),
        m_NumberFeatures( 0)
      {
      }

      //! @brief constructor from monitor data and weighting function
      ObjectiveFunctionAucRocCurve
      (
        const float CUTOFF,
        const math::Polynomial &WEIGHTING_FUNCTION,
        const bool POSITIVES_ABOVE_THRESHOLD
      ) :
        m_Cutoff( CUTOFF),
        m_AucWeightingFunction( WEIGHTING_FUNCTION),
        m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD),
        m_XAxisLogScale( false),
        m_FPRCutoffRange( 0.0, 1.0),
        m_NumberFeatures( 0)
      {
      }

      //! @brief Clone function
      //! @return pointer to new ObjectiveFunctionAucRocCurve
      ObjectiveFunctionAucRocCurve *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

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

      //! @brief set the threshold
      //! @param THRESHOLD threshold that divides output classes
      void SetThreshold( const float &THRESHOLD)
      {
        m_Cutoff = THRESHOLD;
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

      //! @brief set weighting function for auc of roc curve
      //! @param WEIGHTING_FUNCTION weighting function to weight auc of roc curve
      void SetWeightingFunction( const math::Polynomial &WEIGHTING_FUNCTION)
      {
        m_AucWeightingFunction = WEIGHTING_FUNCTION;
      }

      //! @brief get weighting function for auc of roc curve
      const math::Polynomial &GetWeightingFunction() const
      {
        return m_AucWeightingFunction;
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
      float GetDesiredHitRate() const;

      //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
      //! @param DATA monitoring dataset results, non-scaled
      //! @param IDS ids; can be used by the objective function
      void SetData
      (
        const FeatureDataSet< float> &DATA,
        const FeatureDataSet< char> &IDS = FeatureDataSet< char>()
      );

    //////////////////////
    // input and output //
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

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class ObjectiveFunctionAucRocCurve

  } // namespace model

} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_AUC_ROC_CURVE_H_
