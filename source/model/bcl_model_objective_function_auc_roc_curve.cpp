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
#include "model/bcl_model_objective_function_auc_roc_curve.h"

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
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionAucRocCurve::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionAucRocCurve()
      )
    );

    //! copy constructor
    ObjectiveFunctionAucRocCurve *ObjectiveFunctionAucRocCurve::Clone() const
    {
      return new ObjectiveFunctionAucRocCurve( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ObjectiveFunctionAucRocCurve::GetAlias() const
    {
      static const std::string s_Name( "AucRocCurve");
      return s_Name;
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionAucRocCurve::operator()
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
      BCL_Assert
      (
        m_FPRCutoffRange.GetMin() < m_FPRCutoffRange.GetMax(),
        "Min fpr must be < max fpr for integration range"
      );

      // sum the integrals for each result
      float sum_integral( 0);

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
        // weighted auc under roc curve
        double local_integral
        (
          roc_curve.Integral
          (
            m_AucWeightingFunction,
            &math::ContingencyMatrix::GetFalsePositiveRate,
            &math::ContingencyMatrix::GetTruePositiveRate,
            m_FPRCutoffRange,
            m_XAxisLogScale
          )
        );
        sum_integral += local_integral;
        BCL_MessageStd
        (
          "weighted auc for result # " + util::Format()( result_number) + ": " + util::Format()( local_integral) + "; "
          + util::Format()( local_integral / m_ExpectedResult) + " x the naive predictor result"
        );
      }

      // compute the average
      sum_integral /= float( result_size);

      BCL_MessageStd
      (
        "average weighted auc: " + util::Format()( sum_integral) + "; "
        + util::Format()( sum_integral / m_ExpectedResult) + " x the naive predictor result"
      );

      // return average integral of roc curve
      return sum_integral;
    }

    //! @brief get the desired hit rate (e.g. fraction of positive predictions)
    //! @return the desired hit rate, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
    //!       of values will be considered
    float ObjectiveFunctionAucRocCurve::GetDesiredHitRate() const
    {
      BCL_Assert( m_NumberFeatures || !m_XAxisLogScale, "SetData needs to be called before GetDesiredHitRate");

      if( !m_XAxisLogScale)
      {
        // check for default parameters for Min/Max; in which case the default behavior of dynamically selecting
        // based on actual features is likely best
        if( m_FPRCutoffRange.GetMin() <= float( 0.0) && m_FPRCutoffRange.GetMax() >= float( 1.0))
        {
          return util::GetUndefined< float>();
        }
        // return the midpoint of the range
        return m_FPRCutoffRange.GetMiddle();
      }

      // determine adjusted min point
      const float log_fpr_min( std::log( std::max( m_FPRCutoffRange.GetMin(), 1.0 / double( m_NumberFeatures))));
      const float log_fpr_max( std::log( std::max( m_FPRCutoffRange.GetMax(), 1.0 / double( m_NumberFeatures))));

      // return log-mid point
      return std::exp( 0.5 * ( log_fpr_min + log_fpr_max));
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionAucRocCurve::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      m_NumberFeatures = DATA.GetNumberFeatures();

      // list of pairs to compute the expected value
      storage::List< storage::Pair< double, bool> > values_fake;
      values_fake.PushBack( storage::Pair< double, bool>( -1.0, false));
      values_fake.PushBack( storage::Pair< double, bool>( 1.0, true));
      // compute expected value
      math::ROCCurve roc_curve( values_fake);
      m_ExpectedResult =
        roc_curve.Integral
        (
          m_AucWeightingFunction,
          &math::ContingencyMatrix::GetFalsePositiveRate,
          &math::ContingencyMatrix::GetTruePositiveRate,
          m_FPRCutoffRange,
          m_XAxisLogScale
        );
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionAucRocCurve::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the area under the receiver-operator curve"
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
        "x_axis_log",
        "calculate AUC while x axis is on log scale",
        io::Serialization::GetAgent( &m_XAxisLogScale),
        "0"
      );
      parameters.AddInitializer
      (
        "polynomial",
        "coefficient of the polynomial, in increasing degree, e.g. polynomial(1, 2) = 1+2x",
        io::Serialization::GetAgent( &m_AucWeightingFunction),
        "(1)"
      );
      parameters.AddInitializer
      (
        "parity",
        "specifies actives are above or below cutoff, 0 - below, 1 - above",
        io::Serialization::GetAgent( &m_PositivesAboveThreshold),
        "0"
      );
      parameters.AddInitializer
      (
        "min fpr",
        "minimum fpr to begin integration at. "
        "Values smaller than 1/(P+N) will automatically be rounded up to 1/(P+N) if log scaling is set",
        io::Serialization::GetAgentWithRange( &m_FPRCutoffRange.GetMin(), 0.0, 1.0),
        "0"
      );
      parameters.AddInitializer
      (
        "max fpr",
        "maximum fpr to end integration at. ",
        io::Serialization::GetAgentWithRange( &m_FPRCutoffRange.GetMax(), 0.0, 1.0),
        "1"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
