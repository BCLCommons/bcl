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
#include "model/bcl_model_objective_function_bootstrap.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_data_set_select_columns.h"
#include "model/bcl_model_feature_label_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionBootstrap::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionBootstrap()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! copy constructor
    ObjectiveFunctionBootstrap *ObjectiveFunctionBootstrap::Clone() const
    {
      return new ObjectiveFunctionBootstrap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief default constructor
    ObjectiveFunctionBootstrap::ObjectiveFunctionBootstrap() :
      m_NumberBootstraps( 2000),
      m_ConfidenceInterval( 0.95)
    {
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectiveFunctionBootstrap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief determine what sign of the derivative of this objective function indicates improvement
    //! @return the sign of the derivative of this objective function indicates improvement
    opti::ImprovementType ObjectiveFunctionBootstrap::GetImprovementType() const
    {
      if( m_Objective.IsDefined())
      {
        return m_Objective->GetImprovementType();
      }
      return opti::e_LargerIsBetter;
    }

    //! @brief get the overall goal of the objective function
    //! @return the goal of the objective function
    ObjectiveFunctionInterface::Goal ObjectiveFunctionBootstrap::GetGoalType() const
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
    linal::Matrix< char> ObjectiveFunctionBootstrap::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      return m_Objective->GetFeaturePredictionClassifications( EXPERIMENTAL, PREDICTED);
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionBootstrap::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      m_RealIds = IDS;
      m_Objective->SetData( DATA, m_RealIds);
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionBootstrap::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      const size_t dataset_size( EXPERIMENTAL.GetNumberFeatures());
      FeatureDataSet< float> experimental_sample( EXPERIMENTAL);
      FeatureDataSet< float> predicted_sample( PREDICTED);
      FeatureDataSet< char> ids_sample( dataset_size, m_RealIds.GetFeatureSize(), char( 0));

      math::RunningAverageSD< float> ave_val;
      storage::Vector< float> objective_output;
      objective_output.AllocateMemory( m_NumberBootstraps);
      util::Implementation< ObjectiveFunctionInterface> obj_copy( m_Objective);
      for( size_t run_number( 0); run_number < m_NumberBootstraps; ++run_number)
      {
        for( size_t row( 0); row < dataset_size; ++row)
        {
          const size_t selected_row( random::GetGlobalRandom().Random( dataset_size - 1));
          experimental_sample.GetRawMatrix().ReplaceRow( row, EXPERIMENTAL.GetMatrix().GetRow( selected_row));
          predicted_sample.GetRawMatrix().ReplaceRow( row, PREDICTED.GetMatrix().GetRow( selected_row));
          ids_sample.GetRawMatrix().ReplaceRow( row, m_RealIds.GetMatrix().GetRow( selected_row));
        }
        obj_copy->SetData( experimental_sample, ids_sample);
        const float result( m_Objective->operator ()( experimental_sample, predicted_sample));
        ave_val += result;
        objective_output.PushBack( result);
        util::GetLogger().LogStatus
        (
          util::Format()( run_number) + " / " + util::Format()( m_NumberBootstraps) + " bootstraps finished "
          + util::Format().FFP( 2)( 100.0 * float( run_number + 1) / float( m_NumberBootstraps)) + "% complete"
        );
      }
      const size_t conf_low( ( 1.0 - m_ConfidenceInterval) * m_NumberBootstraps);
      const size_t conf_hi( m_NumberBootstraps - conf_low);
      objective_output.Sort( std::less< float>());

      BCL_MessageCrt
      (
        "ObjFunction: " + m_Objective->GetString()
        + " Ave: " + util::Format()( ave_val.GetAverage())
        + " SD: " + util::Format()( ave_val.GetStandardDeviation())
        + " " + util::Format().FFP( 3)( 100.0 * m_ConfidenceInterval)
        + "% confidence interval: "
        + util::Format()( objective_output( conf_low)) + " - " + util::Format()( objective_output( conf_hi))
      );

      return ave_val.GetAverage();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionBootstrap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes bootstrap (sampling with replacement) computation for the objective function, reports average value "
        "confidence intervals"
      );

      parameters.AddInitializer
      (
        "repeats",
        "number of times to repeat the bootstrap to obtain a confidence interval",
        io::Serialization::GetAgent( &m_NumberBootstraps),
        "2000"
      );

      parameters.AddInitializer
      (
        "function",
        "core objective function to use",
        io::Serialization::GetAgent( &m_Objective)
      );

      parameters.AddInitializer
      (
        "confidence interval",
        "confidence interval to report upper and lower bounds for the objective function (95%)",
        io::Serialization::GetAgentWithRange( &m_ConfidenceInterval, 0.5, 1.0),
        "0.95"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
