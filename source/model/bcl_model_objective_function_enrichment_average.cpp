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
#include "model/bcl_model_objective_function_enrichment_average.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionEnrichmentAverage::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionEnrichmentAverage()
      )
    );

    //! copy constructor
    ObjectiveFunctionEnrichmentAverage *ObjectiveFunctionEnrichmentAverage::Clone() const
    {
      return new ObjectiveFunctionEnrichmentAverage( *this);
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionEnrichmentAverage::operator()
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
      float sum_enrichment_averages( 0);

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

        // average of enrichment values
        math::RunningAverage< float> enrichment_average;

        // calculate enrichment by sampling of different cutoffs
        for( float cutoff( m_EnrichmentCutOffStepSize); cutoff <= m_EnrichmentCutOffMax; cutoff += m_EnrichmentCutOffStepSize)
        {
          const math::ContingencyMatrix matrix( roc_curve.ContingencyMatrixFraction( cutoff));
          const float enrichment( matrix.GetEnrichment());

          // if the enrichment cutoff step size is too small, then there may be 0 predicted actives/inactives, thus making
          // the enrichment undefined.  Therefore, only add defined enrichments
          if( util::IsDefined( enrichment))
          {
            BCL_MessageDbg
            (
              "#cutoff: " + util::Format()( cutoff)
              + " enrichment: " + util::Format()( enrichment) +
              " FP/TP/FN/TN: "
              + util::Format()( matrix.GetNumberFalsePositives()) + "/"
              + util::Format()( matrix.GetNumberTruePositives()) + "/"
              + util::Format()( matrix.GetNumberFalseNegatives()) + "/"
              + util::Format()( matrix.GetNumberTrueNegatives()) + "/"
            );

            // sum enrichments
            enrichment_average += enrichment;
          }
        }

        // add the averaged enrichment as objective function score
        sum_enrichment_averages += enrichment_average.GetAverage();
      }

      // return averaged enrichment as objective function score
      return sum_enrichment_averages / float( result_size);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionEnrichmentAverage::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the average enrichment over a range of cutoffs"
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
        "enrichment max",
        "maximum cutoff to use for calculating enrichment",
        io::Serialization::GetAgentWithRange( &m_EnrichmentCutOffMax, 0.0, 1.0),
        "1.0"
      );
      parameters.AddInitializer
      (
        "step size",
        "step size for enrichment cutoff, smaller step sizes improve accuracy at minor cost in speed",
        io::Serialization::GetAgentWithRange( &m_EnrichmentCutOffStepSize, 0.0, 1.0),
        "0.1"
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
