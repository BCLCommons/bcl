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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_CATEGORICAL_MAX_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_CATEGORICAL_MAX_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

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
    //! @class ObjectiveFunctionCategoricalMax
    //! @brief determines the score of a Decision Tree.
    //!
    //! @see @link example_model_objective_function_categorical_max.cpp @endlink
    //! @author mendenjl
    //! @date Jul 18, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionCategoricalMax :
      public ObjectiveFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      linal::Vector< size_t> m_ClassBoundaries;   //!< End index for each classification boundary
      linal::Matrix< size_t> m_ResultsBestIndex;  //!< Index of each classication for each feature
      mutable linal::Matrix< size_t> m_PredictedClasses; //!< Storage for predicted class matrix

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new Evaluator
      ObjectiveFunctionCategoricalMax *Clone() const;

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
        static const std::string s_Name( "CategoricalMax");
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
        return 0.5;
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

    ///////////////
    // operators //
    ///////////////

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

      //! @brief get the class boundaries
      //! @return indices of the class boundaries
      const linal::Vector< size_t> &GetClassBoundaries() const
      {
        return m_ClassBoundaries;
      }

      //! @brief get the results classes
      //! @return class for each result class
      const linal::Matrix< size_t> &GetActualClasses() const
      {
        return m_ResultsBestIndex;
      }

      //! @brief get the predicted classes (from the most recent call to operator())
      //! @return the predicted classes (from the most recent call to operator())
      const linal::Matrix< size_t> &GetPredictedClasses() const
      {
        return m_PredictedClasses;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_OBJECTIVE_FUNCTION_CATEGORICAL_MAX_H_
