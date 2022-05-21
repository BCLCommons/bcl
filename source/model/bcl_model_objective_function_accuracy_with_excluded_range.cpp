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
#include "model/bcl_model_objective_function_accuracy_with_excluded_range.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionAccuracyWithExcludedRange::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionAccuracyWithExcludedRange()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionAccuracyWithExcludedRange::ObjectiveFunctionAccuracyWithExcludedRange() :
      m_ActivityCutoff( 0.0),
      m_ActivityCutoffLow( 0.0),
      m_ActivityCutoffHigh( 0.0)
    {
    }

    //! @brief constructor
    ObjectiveFunctionAccuracyWithExcludedRange::ObjectiveFunctionAccuracyWithExcludedRange( const float &ACTIVITY_CUTOFF) :
      m_ActivityCutoff( ACTIVITY_CUTOFF),
      m_ActivityCutoffLow( ACTIVITY_CUTOFF),
      m_ActivityCutoffHigh( ACTIVITY_CUTOFF)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Evaluator
    ObjectiveFunctionAccuracyWithExcludedRange *ObjectiveFunctionAccuracyWithExcludedRange::Clone() const
    {
      return new ObjectiveFunctionAccuracyWithExcludedRange( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ObjectiveFunctionAccuracyWithExcludedRange::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    linal::Matrix< char> ObjectiveFunctionAccuracyWithExcludedRange::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // number of data points in dataset
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());

      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == PREDICTED.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );
      BCL_Assert
      (
        EXPERIMENTAL.GetFeatureSize() == PREDICTED.GetFeatureSize(),
        "Experimental and predicted values have different result sizes"
      );

      linal::Matrix< char> prediction_classification( data_set_size, result_size, '\0');

      // iterate over all experimental and predicted values and check for accuracy
      for( size_t counter( 0); counter < data_set_size; ++counter)
      {
        // iterate over each result
        for( size_t result_number( 0); result_number < result_size; ++result_number)
        {
          const float res( EXPERIMENTAL( counter)( result_number));
          if( res < m_ActivityCutoffLow || res > m_ActivityCutoffHigh)
          {
            // check if experimental was above or below the cutoff
            const bool experimental( res < m_ActivityCutoff);
            // check if predicted was above or below the cutoff
            const bool prediction( PREDICTED( counter)( result_number) < m_ActivityCutoff);

            // if experimental and prediction belong to same group, the prediction was correct
            prediction_classification( counter, result_number) = char( experimental == prediction ? 'P' : 'n');
          }
        }
      }

      // return accuracy
      return prediction_classification;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionAccuracyWithExcludedRange::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // number of data points in dataset
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());

      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == PREDICTED.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );
      BCL_Assert
      (
        EXPERIMENTAL.GetFeatureSize() == PREDICTED.GetFeatureSize(),
        "Experimental and predicted values have different result sizes"
      );

      // count the number of predictions that were on the same side of the cutoff as the experimental predictions
      size_t accurate_prediction_count( 0);

      // count the number of predictions that were outside the specified range
      size_t count_outside_range( 0);

      // iterate over all experimental and predicted values and check for accuracy
      for( size_t counter( 0); counter < data_set_size; ++counter)
      {
        // iterate over each result
        for( size_t result_number( 0); result_number < result_size; ++result_number)
        {
          const float res( EXPERIMENTAL( counter)( result_number));
          if( res < m_ActivityCutoffLow || res > m_ActivityCutoffHigh)
          {
            // check if experimental was above or below the cutoff
            const bool experimental( EXPERIMENTAL( counter)( result_number) < m_ActivityCutoff);
            // check if predicted was above or below the cutoff
            const bool prediction( PREDICTED( counter)( result_number) < m_ActivityCutoff);

            // if experimental and prediction belong to same group, the prediction was correct
            if( experimental == prediction)
            {
              ++accurate_prediction_count;
            }
            ++count_outside_range;
          }
        }
      }

      // return accuracy
      return float( accurate_prediction_count) / float( count_outside_range) / float( result_size);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionAccuracyWithExcludedRange::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the fraction of predictions that were correct using the given cutoff"
      );

      parameters.AddInitializer
      (
        "cutoff",
        "result value that separates correct from incorrect results",
        io::Serialization::GetAgent( &m_ActivityCutoff),
        "0"
      );
      parameters.AddInitializer
      (
        "cutoff low",
        "result value that separates correct from incorrect results",
        io::Serialization::GetAgent( &m_ActivityCutoffLow),
        "0"
      );
      parameters.AddInitializer
      (
        "cutoff high",
        "result value that separates correct from incorrect results",
        io::Serialization::GetAgent( &m_ActivityCutoffHigh),
        "0"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
