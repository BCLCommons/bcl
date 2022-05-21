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
#include "model/bcl_model_objective_function_integral_tnr_tpr.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_roc_curve.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionIntegralTnrTpr::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionIntegralTnrTpr()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionIntegralTnrTpr::ObjectiveFunctionIntegralTnrTpr() :
      m_Cutoff(),
      m_PositivesAboveThreshold()
    {
    }

    //! @brief constructor from monitor data
    //! @param CUTOFF cutoff to determining the area under the curve
    //! @param FRACTION_PRED_POS_CUTOFF cuoff for fraction predicted positives
    //! @param POSITIVES_ABOVE_THRESHOLD flag for sorting list of predicted and experimental values in roc curve
    ObjectiveFunctionIntegralTnrTpr::ObjectiveFunctionIntegralTnrTpr
    (
      const float CUTOFF,
      const bool POSITIVES_ABOVE_THRESHOLD
    ) :
      m_Cutoff( CUTOFF),
      m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ObjectiveFunctionIntegralTnrTpr
    ObjectiveFunctionIntegralTnrTpr *
    ObjectiveFunctionIntegralTnrTpr::Clone() const
    {
      return new ObjectiveFunctionIntegralTnrTpr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionIntegralTnrTpr::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionIntegralTnrTpr::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      // number of data points in dataset
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());

      // sum the precisions and precision idel for each result
      float sum_precision( 0);

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

        sum_precision += roc_curve.Integral
          (
            &math::ContingencyMatrix::GetTrueNegativeRate,           // x-coordinate
            &math::ContingencyMatrix::GetTruePositiveRate,           // y-coordinate
            math::Range< double>(0,1),                               // x-coordinate cutoff range
            false                                                    // x-axis log10 scaling on
          );
      }

      // return the average precision
      return sum_precision / float( result_size);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionIntegralTnrTpr::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates integral of plotting true negative rate (TNR) vs true positive rate (TPR); "
        "TNR =  #True-Negatives / (#False-Positives + #True-Negatives), "
        "TPR =  #True-Positives / (#True-Positives + #False-Negatives)"
      );

      parameters.AddInitializer
      (
        "cutoff",
        "potency cutoff that separates actives from inactives",
        io::Serialization::GetAgent( &m_Cutoff)
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
