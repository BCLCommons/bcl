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
#include "model/bcl_model_objective_function_integral_precision_fraction_predicted.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_allowed.h"
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
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionIntegralPrecisionFractionPredicted::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionIntegralPrecisionFractionPredicted()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionIntegralPrecisionFractionPredicted::ObjectiveFunctionIntegralPrecisionFractionPredicted() :
      m_Cutoff(),
      m_FPPCutoffRange(),
      m_PositivesAboveThreshold(),
      m_DatasetSize( 0)
    {
    }

    //! @brief constructor from monitor data
    //! @param CUTOFF cutoff to determining the area under the curve
    //! @param FRACTION_PRED_POS_CUTOFF cuoff for fraction predicted positives
    //! @param POSITIVES_ABOVE_THRESHOLD flag for sorting list of predicted and experimental values in roc curve
    ObjectiveFunctionIntegralPrecisionFractionPredicted::ObjectiveFunctionIntegralPrecisionFractionPredicted
    (
      const float CUTOFF,
      const math::Range< double> &FRACTION_PRED_POS_CUTOFF,
      const bool POSITIVES_ABOVE_THRESHOLD
    ) :
      m_Cutoff( CUTOFF),
      m_FPPCutoffRange( FRACTION_PRED_POS_CUTOFF),
      m_PositivesAboveThreshold( POSITIVES_ABOVE_THRESHOLD),
      m_DatasetSize( 0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ObjectiveFunctionIntegralPrecisionFractionPredicted
    ObjectiveFunctionIntegralPrecisionFractionPredicted *
    ObjectiveFunctionIntegralPrecisionFractionPredicted::Clone() const
    {
      return new ObjectiveFunctionIntegralPrecisionFractionPredicted( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the desired hit rate (e.g. fraction of positive predictions)
    //! @return the desired hit rate, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
    //!       of values will be considered
    float ObjectiveFunctionIntegralPrecisionFractionPredicted::GetDesiredHitRate() const
    {
      // compute average exponent:
      const double average_exponent
      (
        m_FPPCutoffRange.GetMin() == double( 0)
        ? std::log10( m_FPPCutoffRange.GetMax()) + std::log10( double( m_DatasetSize - 1))
        : std::log10( m_FPPCutoffRange.GetMax()) - std::log10( m_FPPCutoffRange.GetMin())
      );
      // compute average value
      return pow( double( 10), average_exponent);
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionIntegralPrecisionFractionPredicted::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      m_DatasetSize = DATA.GetNumberFeatures();

      if( m_FPPCutoffType == "compound")
      {
        m_FPPCutoffRange = math::Range< double>
        (
          m_FPPCompoundCutoffRange.GetMin() / m_DatasetSize,
          m_FPPCompoundCutoffRange.GetMax() / m_DatasetSize
        );

        BCL_MessageDbg
        (
          "Converting compound range into x axis pct range. New range: "
          + util::Format()( m_FPPCutoffRange)
        );
      }
    }

    //! @brief calculate integral of ideal PPV vs FPP curve. Necessary for normalization purposes of a calculated curve.
    //! @param ROC_CURVE calculated roc curve of interest
    //! @param X_AXIS_RANGE integration range on x axis (FPP)
    //! @return integral of ideal PPV vs FPP curve
    float ObjectiveFunctionIntegralPrecisionFractionPredicted::CalculateIdealIntegral
    (
      const math::ROCCurve ROC_CURVE,
      const math::Range< double> X_AXIS_RANGE
    ) const
    {
      const size_t total_results( ROC_CURVE.GetNumberResults());
      const double num_actual_positives( ROC_CURVE.GetNumberActualPositives());
      const double step_size( 1.0 / double( total_results - 1));
      double prev_progress( 0);

      bool initialize( true);

      // sum the precisions and precision idel for each result
      float sum_precision_ideal( 0);

      // iterate over sorted counts of roc curve and plot ideal curve for PPV vs FractionPositivePredicted
      for
      (
        storage::Vector< math::ROCCurve::Point>::const_iterator
          itr( ROC_CURVE.GetSortedCounts().Begin()), itr_end( ROC_CURVE.GetSortedCounts().End());
        itr != itr_end;
        ++itr
      )
      {
        const size_t counts( itr->GetNumberPredictedPositives());
        const double progress( step_size * counts);

        if( X_AXIS_RANGE.GetMin() < progress && X_AXIS_RANGE.GetMax() >= progress)
        {
          if( initialize)
          {
            prev_progress = progress;
            initialize = false;
            continue;
          }
          sum_precision_ideal += std::log10( progress / prev_progress) * std::min( 1.0, num_actual_positives / ( counts + 1));
          prev_progress = progress;
        }
      }

      sum_precision_ideal /=
          ( X_AXIS_RANGE.GetMin() == double( 0))
          ? std::log10( X_AXIS_RANGE.GetMax()) - std::log10( step_size)
          : std::log10( X_AXIS_RANGE.GetMax()) - std::log10( X_AXIS_RANGE.GetMin());

      return sum_precision_ideal;
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionIntegralPrecisionFractionPredicted::operator()
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

      // sum the precisions and precision idel for each result
      float sum_precision( 0), sum_precision_ideal( 0);

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

        sum_precision +=
          roc_curve.Integral
          (
            &math::ContingencyMatrix::GetFractionPredictedPositives, // x-coordinate
            &math::ContingencyMatrix::GetPrecision,                  // y-coordinate
            m_FPPCutoffRange,                                        // x-coordinate cutoff range
            true                                                     // x-axis log10 scaling on
          );

        // calculate integral of corresponding ideal ppv vs fpp curve
        sum_precision_ideal = CalculateIdealIntegral( roc_curve, m_FPPCutoffRange);

        if( sum_precision_ideal > 0)
        {
          sum_precision /= sum_precision_ideal;
        }
        else
        {
          BCL_MessageStd( "Ideal FPPvsPPV integral is 0!");
          sum_precision = float( 0);
        }
      }

      // return the average precision
      return sum_precision / float( result_size);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionIntegralPrecisionFractionPredicted::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates integral of plot Precision vs Fraction Predicted Positives (FPP); "
        "Precision = #True-Positives / (#True-Positives + #False-Positives), "
        "FPP = ( #True-Positives + #False-Positives) / (#Positives + #Negatives)"
      );

      parameters.AddInitializer
      (
        "cutoff",
        "potency cutoff that separates actives from inactives",
        io::Serialization::GetAgent( &m_Cutoff)
      );
      parameters.AddInitializer
      (
        "cutoff_type",
        "x_axis_range determined by a compound range (compound) or a percent range on x axis (fpp_percent)",
        io::Serialization::GetAgentWithCheck
        (
          &m_FPPCutoffType,
          command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "compound", "fpp_percent"))
        )
      );
      parameters.AddInitializer
      (
        "x_axis_range",
        "Fraction Predicted Positives (FPP) cutoff range on x-axis for which integral should be calculated",
        io::Serialization::GetAgent( &m_FPPCutoffRange),
        "x_axis_range=\"[0.0 , 1.0]\""
      );
      parameters.AddInitializer
      (
        "x_axis_compound_range",
        "Fraction Predicted Positives (FPP) cutoff compound range on x-axis for which integral should be calculated"
        "the compound range translates into a x_axis_range between 0 and 1, it super-seeds x_axis_range",
        io::Serialization::GetAgent( &m_FPPCompoundCutoffRange),
        "x_axis_compound_range=\"[25 , 250]\""
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
