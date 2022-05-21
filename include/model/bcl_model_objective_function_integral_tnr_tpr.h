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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_INTEGRAL_TNR_TPR_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_INTEGRAL_TNR_TPR_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

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
    //! @class ObjectiveFunctionIntegralTnrTpr
    //! @brief objective function that evaluates precision in comparison to fraction predicted positives
    //! @details The integral of Precision (PPV) versus fraction predicted positives (FPP=(TP+FP)/(P+N))
    //!          is calculated until a given cutoff is reached. This is similar to a classical roc curve plot only that
    //!          the x-coordinate is PPV and the y-coordinate is FPP.
    //!          (in roc curve x-coordinate is false positive rate FPR and y-coordinate is TPR)
    //!
    //! @see @link example_model_objective_function_integral_tnr_tpr.cpp @endlink
    //! @author butkiem1
    //! @date Aug 06, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionIntegralTnrTpr :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! potency cutoff to distinguish actives from inactives
      float m_Cutoff;

      //! flag for sorting list of predicted and experimental values in roc curve
      bool m_PositivesAboveThreshold;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionIntegralTnrTpr();

      //! @brief constructor from monitor data
      ObjectiveFunctionIntegralTnrTpr
      (
        const float CUTOFF,
        const bool POSITIVES_ABOVE_THRESHOLD
      );

      //! @brief Clone function
      //! @return pointer to new ObjectiveFunctionIntegralTnrTpr
      ObjectiveFunctionIntegralTnrTpr *Clone() const;

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
        static const std::string s_Name( "TnrTprAuc");
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

    }; // class ObjectiveFunctionIntegralTnrTpr

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_OBJECTIVE_FUNCTION_INTEGRAL_TNR_TPR_H_
