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
#include "model/bcl_model_objective_function_enrichment.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_roc_curve.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionEnrichment::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionEnrichment()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionEnrichment::ObjectiveFunctionEnrichment() :
      m_Cutoff( 0.0),
      m_OverSamplingFactor( 1.0),
      m_PositivesAboveThreshold( true)
    {
    }

    //! @brief constructor from monitor data
    ObjectiveFunctionEnrichment::ObjectiveFunctionEnrichment
    (
      const float CUTOFF,
      const float OVERSAMPLING_FACTOR,
      const bool POSITIVES_ABOVE_THRESHOLD
    ) :
      m_Cutoff( CUTOFF),
      m_OverSamplingFactor( OVERSAMPLING_FACTOR),
      m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD)
    {
    }

    //! copy constructor
    ObjectiveFunctionEnrichment *ObjectiveFunctionEnrichment::Clone() const
    {
      return new ObjectiveFunctionEnrichment( *this);
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionEnrichment::operator()
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

      // sum the enrichments for each result
      float sum_enrichment( 0);

      for( size_t result_number( 0); result_number < result_size; ++result_number)
      {
        // list of pairs with pred and exp values for ROC Curve
        storage::List< storage::Pair< double, double> > values_predicted_experimental;

        // iterate over all experimental and predicted values and check for accuracy
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          values_predicted_experimental.PushBack
          (
            storage::Pair< double, double>
            (
              PREDICTED( counter)( result_number),
              EXPERIMENTAL( counter)( result_number)
            )
          );
        }

        // create roc curve according to cutoff
        math::ROCCurve roc_curve( values_predicted_experimental, m_Cutoff, m_PositivesAboveThreshold);

        // write out contingency matrix for top fraction of predicted in roc curve
        BCL_MessageDbg
        (
          "FPR cutoff: "
          + util::Format()( m_FalsePositiveRateCutoff)
          + util::Format()( " oversampling factor: ")
          + util::Format()( m_OverSamplingFactor)
          + " ConMatx: "
          + util::Format()( roc_curve.ContingencyMatrixFraction( m_FalsePositiveRateCutoff))
        );

        // sum enrichment of top fraction according to oversampling factor
        sum_enrichment += roc_curve.ContingencyMatrixFraction( m_FalsePositiveRateCutoff).GetOversampledEnrichment( m_OverSamplingFactor);
      }

      // return the average enrichment
      return sum_enrichment / float( result_size);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionEnrichment::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the enrichment = Precision / Positive Sampling Rate; Precision = #True-Positives / (#True-Positives + #False-Positives), Positive Sampling Rate = #Positives / Sample-Size"
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
        "oversampling factor",
        "factor how many times a group of compounds is repeated",
        io::Serialization::GetAgent( &m_OverSamplingFactor),
        "1"
      );
      parameters.AddInitializer
      (
        "FPR cutoff",
        "false positive rate (FPR) cutoff for which enrichment should be calculated",
        io::Serialization::GetAgent( &m_FalsePositiveRateCutoff),
        "0.03"
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
