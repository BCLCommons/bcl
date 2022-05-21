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
#include "model/bcl_model_objective_function_contingency_matrix_measure.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionContingencyMatrixMeasure::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionContingencyMatrixMeasure()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionContingencyMatrixMeasure::ObjectiveFunctionContingencyMatrixMeasure() :
      m_ActivityCutoff( 0.0),
      m_PositivesAboveThreshold( false),
      m_OptimizePredictedCutoff( false)
    {
    }

    //! @brief constructor taking all necessary parameters
    //! @param ACTIVITY_CUTOFF the activity cutoff
    //! @param MEASURE the measure to compute
    ObjectiveFunctionContingencyMatrixMeasure::ObjectiveFunctionContingencyMatrixMeasure
    (
      const float &ACTIVITY_CUTOFF,
      const math::ContingencyMatrixMeasures &MEASURE,
      const bool &POSITIVES_ABOVE_THRESHOLD,
      const bool &OPTIMIZE_CUTOFF
    ) :
      m_ActivityCutoff( ACTIVITY_CUTOFF),
      m_Measure( MEASURE),
      m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD),
      m_OptimizePredictedCutoff( OPTIMIZE_CUTOFF)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Evaluator
    ObjectiveFunctionContingencyMatrixMeasure *ObjectiveFunctionContingencyMatrixMeasure::Clone() const
    {
      return new ObjectiveFunctionContingencyMatrixMeasure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ObjectiveFunctionContingencyMatrixMeasure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionContingencyMatrixMeasure::operator()
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

      // compute the average of the measure over all outputs
      math::RunningAverage< float> average_measure;

      // generate an output message to inform the user of progress
      std::ostringstream output_info;

      // iterate over each result
      for( size_t result_number( 0); result_number < result_size; ++result_number)
      {
        if( !m_OptimizePredictedCutoff)
        {
          // count the number of predictions of each type in TP/FP/TN/FN
          size_t tp( 0), tn( 0), fp( 0), fn( 0);
          // iterate over all experimental and predicted values and check for accuracy
          for( size_t counter( 0); counter < data_set_size; ++counter)
          {
            if( util::IsDefined( EXPERIMENTAL( counter)( result_number)))
            {
              // check if experimental was above or below the cutoff
              const bool experimental( EXPERIMENTAL( counter)( result_number) >= m_ActivityCutoff);
              // check if predicted was above or below the cutoff
              const bool prediction( PREDICTED( counter)( result_number) >= m_ActivityCutoff);

              // if experimental and prediction belong to same group, the prediction was correct
              if( experimental == prediction)
              {
                ++( experimental == m_PositivesAboveThreshold ? tp : tn);
              }
              else
              {
                ++( prediction == m_PositivesAboveThreshold ? fp : fn);
              }
            }
          }
          const math::ContingencyMatrix matrix( tp, fp, fn, tn);
          if( matrix.GetTotal())
          {
            const float measure( ( *m_Measure)( matrix));
            average_measure += measure;
            output_info << "\tResult " << result_number << " TP/TN/FP/FN: "
                        << tp << ' ' << tn << ' ' << fp << ' ' << fn << ' '
                        << m_Measure.GetAlias() << ": " << measure;
          }
        }
        else
        {
          // find the optimal value for the measure
          storage::List< storage::Pair< double, double> > results_list;
          // iterate over all experimental and predicted values and check for accuracy
          for( size_t counter( 0); counter < data_set_size; ++counter)
          {
            if( util::IsDefined( EXPERIMENTAL( counter)( result_number)))
            {
              results_list.PushBack
              (
                storage::Pair< double, double>
                (
                  PREDICTED( counter)( result_number),
                  EXPERIMENTAL( counter)( result_number)
                )
              );
            }
          }
          // create a roc curve
          math::ROCCurve curve( results_list, m_ActivityCutoff, m_PositivesAboveThreshold);

          if( !curve.GetSortedCounts().IsEmpty())
          {
            // find the optimal value for the given contingency matrix measure
            // first, get the last element of the curve, which is needed for computing the
            const math::ROCCurve::Point &roc_end_point( curve.GetSortedCounts().LastElement());

            std::pair< storage::Vector< math::ROCCurve::Point>::const_iterator, double>
              best_cutoff_and_metric_value( curve.GetMaxima( *m_Measure));
            const math::ContingencyMatrix best_contingency_matrix
            (
              best_cutoff_and_metric_value.first->GetContingencyMatrix( roc_end_point)
            );
            const float best_metric_value( best_cutoff_and_metric_value.second);
            const float best_cutoff( best_cutoff_and_metric_value.first->GetCutoff());

            // output the best contingency matrix info
            average_measure += best_metric_value;
            output_info << "\tResult " << result_number << " TP/TN/FP/FN: "
                        << best_contingency_matrix.GetNumberTruePositives() << ' '
                        << best_contingency_matrix.GetNumberTrueNegatives() << ' '
                        << best_contingency_matrix.GetNumberFalsePositives() << ' '
                        << best_contingency_matrix.GetNumberFalseNegatives() << ' '
                        << m_Measure.GetAlias() << ": " << best_metric_value
                        << " @cutoff " << best_cutoff;
          }
        }
      }

      BCL_MessageStd( output_info.str());

      // return accuracy
      return average_measure.GetAverage();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionContingencyMatrixMeasure::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates any contingency matrix measure at a specific cutoff. "
        "Multiple output columns will result in the average of the measure being reported."
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
        "measure",
        "contingency matrix measure to consider for each result",
        io::Serialization::GetAgent( &m_Measure)
      );
      parameters.AddInitializer
      (
        "parity",
        "true if positives are above the cutoff",
        io::Serialization::GetAgent( &m_PositivesAboveThreshold)
      );
      parameters.AddInitializer
      (
        "adjustable cutoff",
        "Useful when the predicted cutoff is allowed to be different than the experimental cutoff. This is only "
        "appropriate for measures that consider all aspects of the contingency matrix; e.g. FPR is always 0 at the "
        " start of the ROC curve. In this case, the predicted cutoff will be optimized",
        io::Serialization::GetAgent( &m_OptimizePredictedCutoff),
        "False"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
