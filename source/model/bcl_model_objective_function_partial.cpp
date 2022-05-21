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
#include "model/bcl_model_objective_function_partial.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionPartial::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionPartial()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! copy constructor
    ObjectiveFunctionPartial *ObjectiveFunctionPartial::Clone() const
    {
      return new ObjectiveFunctionPartial( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectiveFunctionPartial::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief determine what sign of the derivative of this objective function indicates improvement
    //! @return the sign of the derivative of this objective function indicates improvement
    opti::ImprovementType ObjectiveFunctionPartial::GetImprovementType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Weight >= 0.0
               ? m_Objective->GetImprovementType()
               : m_Objective->GetImprovementType() == opti::e_LargerEqualIsBetter
                 ? opti::e_SmallerEqualIsBetter
                 : opti::e_LargerEqualIsBetter;
      }
      return opti::e_LargerEqualIsBetter;
    }

    //! @brief get the overall goal of the objective function
    //! @return the goal of the objective function
    ObjectiveFunctionInterface::Goal ObjectiveFunctionPartial::GetGoalType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Objective->GetGoalType();
      }
      return e_Other;
    }

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    linal::Matrix< char> ObjectiveFunctionPartial::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      const size_t dataset_size( EXPERIMENTAL.GetNumberFeatures());
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());
      linal::Matrix< char> prediction_classification( dataset_size, result_size, '\0');

      // filter the undesired outputs from the ranges
      const FeatureDataSet< float> experimental_reduced( EXPERIMENTAL.GetMatrix(), m_Ranges);
      const FeatureDataSet< float> predicted_reduced( PREDICTED.GetMatrix(), m_Ranges);

      // get the classifications from the internal objective function
      linal::Matrix< char> internal_prediction_classification
      (
        m_Objective->GetFeaturePredictionClassifications( EXPERIMENTAL, PREDICTED)
      );

      // determine the corresponding index of each value in the range
      storage::Vector< size_t> mapping;
      for( size_t j( 0); j < result_size; ++j)
      {
        if( m_Ranges.IsWithin( j))
        {
          mapping.PushBack( j);
        }
      }

      BCL_Assert( mapping.GetSize() == internal_prediction_classification.GetNumberCols(), "Mapping is broken");

      const size_t mapping_size( mapping.GetSize());

      // copy the values from the internal objective function into the returned matrix at the appropriate positions
      for( size_t i( 0); i < dataset_size; ++i)
      {
        for( size_t j( 0); j < mapping_size; ++j)
        {
          prediction_classification( i, mapping( j)) = internal_prediction_classification( i, j);
        }
      }

      return prediction_classification;
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionPartial::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      const FeatureDataSet< float> experimental_reduced( DATA.GetMatrix(), m_Ranges);

      m_Objective->SetData( experimental_reduced, IDS);
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionPartial::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // filter the undesired outputs from the ranges
      const FeatureDataSet< float> experimental_reduced( EXPERIMENTAL.GetMatrix(), m_Ranges);
      const FeatureDataSet< float> predicted_reduced( PREDICTED.GetMatrix(), m_Ranges);

      return m_Weight * m_Objective->operator ()( experimental_reduced, predicted_reduced);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionPartial::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "A weighted meta-objective function that uses only a subset of the result data"
      );

      parameters.AddInitializer
      (
        "weight",
        "amount by which to weight the objective function's output",
        io::Serialization::GetAgent( &m_Weight),
        "1.0"
      );

      parameters.AddInitializer
      (
        "function",
        "core objective function to use",
        io::Serialization::GetAgent( &m_Objective)
      );
      parameters.AddInitializer
      (
        "outputs",
        "ranges of chunks to load, e.g. outputs=\"[ 0, 5)+(7,10)\"",
        io::Serialization::GetAgent( &m_Ranges),
        "[0]"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
