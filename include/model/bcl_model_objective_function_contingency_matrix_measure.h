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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_CONTINGENCY_MATRIX_MEASURE_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_CONTINGENCY_MATRIX_MEASURE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_interface.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "math/bcl_math_contingency_matrix.h"
#include "util/bcl_util_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionContingencyMatrixMeasure
    //! @brief Calculates the fraction of predictions that were correct using the given cutoff
    //!
    //! @see @link example_model_objective_function_contingency_matrix_measure.cpp @endlink
    //! @author mendenjl
    //! @date Nov 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionContingencyMatrixMeasure :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! activity cutoff to destinguish between active and inactive compounds
      float m_ActivityCutoff;

      //! measure of interest
      util::Implementation< util::FunctionInterfaceSerializable< math::ContingencyMatrix, double> > m_Measure;

      //! parity
      bool m_PositivesAboveThreshold;

      //! Whether to optimize the activity cutoff for the predicted values
      //! This is really only useful with contingency matrix measures that make use of all values in the contingency matrix,
      //! for example, accuracy, information gain, matthews correlation coefficient, and (possibly) enrichment
      bool m_OptimizePredictedCutoff;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionContingencyMatrixMeasure();

      //! @brief constructor taking all necessary parameters
      //! @param ACTIVITY_CUTOFF the activity cutoff
      //! @param MEASURE the measure to compute
      ObjectiveFunctionContingencyMatrixMeasure
      (
        const float &ACTIVITY_CUTOFF,
        const math::ContingencyMatrixMeasures &MEASURE,
        const bool &POSITIVES_ABOVE_THRESHOLD = false,
        const bool &OPTIMIZE_CUTOFF = false
      );

      //! @brief Clone function
      //! @return pointer to new Evaluator
      ObjectiveFunctionContingencyMatrixMeasure *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_Name( "ContingencyMatrixMeasure");
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
        return m_OptimizePredictedCutoff ? e_RankClassification : e_Classification;
      }

      //! @brief get the threshold, for classification type objectives
      //! @return the threshold, for classification type objectives
      float GetThreshold() const
      {
        return m_ActivityCutoff;
      }

      //! @brief set the threshold
      //! @param THRESHOLD threshold that divides output classes
      void SetThreshold( const float &THRESHOLD)
      {
        m_ActivityCutoff = THRESHOLD;
      }

      //! @brief get the parity, for rank classification type objectives
      //! @return the parity, for rank classification type objectives
      //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
      //!       in prediction of values higher than the threshold, false means below threshold
      bool GetRankingParity() const
      {
        return m_PositivesAboveThreshold;
      }

    ////////////////
    // operators //
    ////////////////

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

    };

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_OBJECTIVE_FUNCTION_CONTINGENCY_MATRIX_MEASURE_H_
