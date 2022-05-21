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
#include "model/bcl_model_objective_function_information_gain_ratio.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionInformationGainRatio::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionInformationGainRatio()
      )
    );

    //! copy constructor
    ObjectiveFunctionInformationGainRatio *ObjectiveFunctionInformationGainRatio::Clone() const
    {
      return new ObjectiveFunctionInformationGainRatio( *this);
    }

    //! @brief evaluate for a given dataset with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionInformationGainRatio::operator()
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
        "Experimental and predicted values have different result sizes, experimental size: "
        + util::Format()( EXPERIMENTAL.GetFeatureSize()) + ", predicted: " + util::Format()( PREDICTED.GetFeatureSize())
      );

      // sum the enrichments for each result
      math::RunningAverage< float> avg_best_contingency_matrix_measure;

      for( size_t result_number( 0); result_number < result_size; ++result_number)
      {
        // list of pairs with pred and exp values for ROC Curve
        storage::List< storage::Pair< double, double> > values_predicted_experimental;

        // iterate over all experimental and predicted values and check for accuracy
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          values_predicted_experimental.PushBack
          (
            storage::Pair< double, double>( PREDICTED( counter)( result_number), EXPERIMENTAL( counter)( result_number))
          );
        }

        // create roc curve according to cutoff
        math::ROCCurve roc_curve( values_predicted_experimental, m_Cutoff, m_PositivesAboveThreshold);

        storage::Vector< math::ROCCurve::Point>::const_iterator itr_sorted_counts( roc_curve.GetSortedCounts().Begin());
        storage::Vector< math::ROCCurve::Point>::const_iterator itr_sorted_counts_end( roc_curve.GetSortedCounts().End());

        // maximum information gain ratio seen so far
        float max_information_gain_ratio( 0.0);

        // keep track of best performing contingency matrix
        math::ContingencyMatrix best_contingency_matrix;

        // iterate over all contingency matrix cutoffs
        for( ; itr_sorted_counts != itr_sorted_counts_end; ++itr_sorted_counts)
        {
          if( itr_sorted_counts->GetNumberPredictedPositives() <= roc_curve.GetNumberActualPositives())
          {
            continue;
          }

          const math::ContingencyMatrix contingency_matrix
          (
            itr_sorted_counts->GetNumberTruePositives(),
            itr_sorted_counts->GetNumberFalsePositives(),
            roc_curve.GetNumberActualPositives() - itr_sorted_counts->GetNumberTruePositives(),
            roc_curve.GetNumberActualNegatives() - itr_sorted_counts->GetNumberFalsePositives()
          );

          const float information_gain_ratio( contingency_matrix.GetInformationGainRatio());

          // check if the current max information gain ratio is reached
          if( max_information_gain_ratio < information_gain_ratio)
          {
            // save the best information gain ratio seen so far
            max_information_gain_ratio = information_gain_ratio;
            // save the best performing contingency matrix
            best_contingency_matrix = contingency_matrix;
          }
        }

        BCL_MessageStd
        (
          "best PPV: " + util::Format()( m_ContingencyMatrixMeasure( best_contingency_matrix))
          + " @infogainratio " + util::Format()( max_information_gain_ratio)
        );

        // add the averaged information gain ratio as objective function score
        avg_best_contingency_matrix_measure += m_ContingencyMatrixMeasure( best_contingency_matrix);
      }

      // return averaged information gain ratio as objective function score
      return avg_best_contingency_matrix_measure.GetAverage();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionInformationGainRatio::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the Information Gain Ratio as an objective function value"
      );
      parameters.AddInitializer
      (
        "measure",
        "contingency matrix measure that will be returned",
        io::Serialization::GetAgent( &m_ContingencyMatrixMeasure),
        "InformationGainRatio"
      );
      parameters.AddInitializer
      (
        "cutoff",
        "result value that separates correct from incorrect results",
        io::Serialization::GetAgent( &m_Cutoff),
        "0"
      );
      parameters.AddInitializer
      (
        "parity",
        "specifies actives are above or below cutoff, 0 - below, 1 - above",
        io::Serialization::GetAgent( &m_PositivesAboveThreshold),
        "0"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
