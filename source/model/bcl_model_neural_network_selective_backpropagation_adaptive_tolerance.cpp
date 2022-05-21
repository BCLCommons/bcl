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
#include "model/bcl_model_neural_network_selective_backpropagation_adaptive_tolerance.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::s_Instance
    (
      util::Enumerated< NeuralNetworkSelectiveBackpropagationInterface>::AddInstance
      (
        new NeuralNetworkSelectiveBackpropagationAdaptiveTolerance()
      )
    );

    //! @brief default constructor
    NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::NeuralNetworkSelectiveBackpropagationAdaptiveTolerance() :
      m_ResultStartID( 0),
      m_ResultEndID( 1000),
      m_MaxTolerance( 0.0),
      m_MinTolerance( 0.0),
      m_InitialTolerance( 0.0),
      m_ToleranceMinusError( 1.0),
      m_AdaptationRate( 0.5),
      m_ReduceErrorByTolerance( false),
      m_PureClassification( false),
      m_ResultsSize( 1000),
      m_NumberThreads( 0)
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkSelectiveBackpropagationAdaptiveTolerance copied from this instance
    NeuralNetworkSelectiveBackpropagationAdaptiveTolerance *NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::Clone() const
    {
      return new NeuralNetworkSelectiveBackpropagationAdaptiveTolerance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::GetAlias() const
    {
      static const std::string s_name( "AdaptiveTolerant");
      return s_name;
    }

    //! @brief Initialize this object with the rescaled dataset and other important parameters
    //! @param RESCALED_DATA the already-rescaled training dataset
    //! @param OBJECTIVE the actual objective function
    //! @param THREADS # of threads that may be accessing this object at once
    void NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::Initialize
    (
      const descriptor::Dataset &RESCALED_DATA,
      const ObjectiveFunctionInterface &OBJECTIVE,
      const size_t &NUMBER_THREADS
    )
    {
      // initialization
      if( m_ResultEndID > RESCALED_DATA.GetResultSize())
      {
        m_ResultEndID = RESCALED_DATA.GetResultSize();
      }
      m_ResultsSize = m_ResultEndID - m_ResultStartID;
      m_NumberThreads = NUMBER_THREADS;
      BCL_Assert( m_MinTolerance <= m_MaxTolerance, "Min tolerance must be less than max tolerance");
      m_RescaledTolerance = linal::Vector< float>( m_ResultsSize, std::max( m_InitialTolerance, m_MinTolerance));
      m_RescaledMaxTolerance = linal::Vector< float>( m_ResultsSize, m_MaxTolerance);
      m_RescaledMinTolerance = linal::Vector< float>( m_ResultsSize, m_MinTolerance);
      m_ThreadResultAverageError.Reset();
      m_ThreadResultAverageError.Resize
      (
        m_NumberThreads,
        storage::Vector< math::RunningAverage< float> >( m_ResultsSize)
      );
      m_ScaledCutoffs = linal::Vector< float>( m_ResultsSize, float( 0.0));
      m_Results = RESCALED_DATA.GetResultsPtr();

      const util::SiPtr< const RescaleFeatureDataSet> results_scaling( RESCALED_DATA.GetResultsPtr()->GetScaling());
      const float rescale_to_range_width( results_scaling->GetRange().GetWidth());
      for( size_t result( 0); result < m_ResultsSize; ++result)
      {
        const float rescaled_range_ratio( results_scaling->GetRescaleRanges()( result).GetWidth() / rescale_to_range_width);
        if( m_PureClassification)
        {
          m_ScaledCutoffs( result) = results_scaling->RescaleValue( result, OBJECTIVE.GetThreshold());
        }
        m_RescaledTolerance( result) /= rescaled_range_ratio;
        m_RescaledMaxTolerance( result) /= rescaled_range_ratio;
        m_RescaledMinTolerance( result) /= rescaled_range_ratio;
      }
    } // Initialize

    //! @brief select whether or not to backpropagate the current feature, and edit the error vector, if desired
    //! @param PREDICTION the prediction that was made by the neural network
    //! @param ERROR Reference to the error vector, (should already have been set to RESULT - PREDICTION)
    //! @param FEATURE_ID id of this feature in the dataset
    //! @param THREAD_ID id of this thread
    //! @return true if the feature should be backpropagated
    bool NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::ShouldBackpropagate
    (
      const linal::VectorConstInterface< float> &PREDICTION,
      linal::VectorInterface< float> &ERROR,
      const size_t &FEATURE_ID,
      const size_t &THREAD_ID
    )
    {
      bool backprop_should_continue( false);
      // pointer to the error vector for this thread
      storage::Vector< math::RunningAverage< float> >::iterator
        itr_error_ave( m_ThreadResultAverageError( THREAD_ID).Begin());
      linal::Vector< float>::const_iterator itr_tol( m_RescaledTolerance.Begin()),
                                            itr_result( ( *m_Results)[ FEATURE_ID] + m_ResultStartID);
      for
      (
        size_t result_offset( 0);
        result_offset < m_ResultsSize;
        ++result_offset, ++itr_error_ave, ++itr_tol, ++itr_result
      )
      {
        // get the actual index
        const size_t result_id( result_offset + m_ResultStartID);
        float &er( ERROR( result_id));

        if( m_PureClassification)
        {
          // for pure classification problems, a predicted value that is farther from the cutoff (and in the same direction)
          // is not an error, and hence should not normally be backpropagated
          const float pred( PREDICTION( result_id));
          const float cutoff( m_ScaledCutoffs( result_offset));
          const int result_over_cutoff( *itr_result > cutoff);
          const int pred_over_cutoff( pred > cutoff);
          const int pred_over_result( pred >= *itr_result);
          const int sum_over( result_over_cutoff + pred_over_cutoff + pred_over_result);
          if( sum_over == 0 || sum_over == 3)
          {
            // either prediction > result > cutoff or
            // prediction < result < cutoff
            er = 0.0;
          }
        }
        const float abs_er( math::Absolute( er));
        const float tolerance( *itr_tol);
        *itr_error_ave += abs_er;
        if( abs_er <= tolerance)
        {
          er = 0.0;
        }
        else
        {
          // absolute error was beyond the tolerance; backpropagate
          backprop_should_continue = true;
          if( m_ReduceErrorByTolerance)
          {
            if( er > 0.0)
            {
              er -= tolerance;
            }
            else
            {
              er += tolerance;
            }
          }
        }
      }

      // return true if any part of ERROR remains to be backpropagated
      return backprop_should_continue;
    }

    //! @brief finalize the current round; occurs only after all threads were already joined
    void NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::FinalizeRound()
    {
      // accumulate results
      for( size_t thread_number( 1); thread_number < m_NumberThreads; ++thread_number)
      {
        for( size_t i( 0); i < m_ResultsSize; ++i)
        {
          m_ThreadResultAverageError( 0)( i) += m_ThreadResultAverageError( thread_number)( i);
        }
      }

      // update tolerance
      for( size_t i( 0); i < m_ResultsSize; ++i)
      {
        const float average_error( m_ThreadResultAverageError( 0)( i).GetAverage());
        BCL_MessageVrb
        (
          "Average error for result " + util::Format()( i)
          + ": " + util::Format()( average_error)
        );
        const float target_tolerance( average_error + m_ToleranceMinusError);
        const float current_tolerance( m_RescaledTolerance( i));
        const float new_tolerance_unbounded
        (
          target_tolerance * m_AdaptationRate + current_tolerance * ( 1.0 - m_AdaptationRate)
        );
        m_RescaledTolerance( i) = std::max( m_RescaledMinTolerance( i), std::min( m_RescaledMaxTolerance( i), new_tolerance_unbounded));
      }

      // reset thread-specific arrays
      for( size_t thread_number( 0); thread_number < m_NumberThreads; ++thread_number)
      {
        for( size_t i( 0); i < m_ResultsSize; ++i)
        {
          m_ThreadResultAverageError( thread_number)( i).Reset();
        }
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkSelectiveBackpropagationAdaptiveTolerance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Always back-propagating features that are mispredicted beyond the tolerance. "
        "If the predicted value has less than the current tolerance of error, backpropagation will not occur. "
        "This implementation adapts the tolerance based on the average error (per result column) each round."
      );
      parameters.AddInitializer
      (
        "begin",
        "result columns to use this triager for",
        io::Serialization::GetAgent( &m_ResultStartID),
        "0"
      );
      parameters.AddInitializer
      (
        "end",
        "1+max result column to consider",
        io::Serialization::GetAgent( &m_ResultEndID),
        "1000"
      );
      parameters.AddInitializer
      (
        "max tolerance",
        "Tolerance will not be allowed to go above this value.  This will also be the initial tolerance value",
        io::Serialization::GetAgent( &m_MaxTolerance)
      );
      parameters.AddInitializer
      (
        "min tolerance",
        "Tolerance will not be allowed to go below this value",
        io::Serialization::GetAgent( &m_MinTolerance),
        "0"
      );
      parameters.AddInitializer
      (
        "initial tolerance",
        "Initial tolerance value.  (If the min tolerance is larger than this value, it will be used instead)",
        io::Serialization::GetAgent( &m_InitialTolerance),
        "0"
      );
      parameters.AddInitializer
      (
        "adaptation rate",
        "Tolerance for each result will move at this rate to the current target tolerance (error multiplier * average error for this result)",
        io::Serialization::GetAgentWithRange( &m_AdaptationRate, 0.0, 1.0),
        "0.5"
      );
      parameters.AddInitializer
      (
        "tolerance offset from error",
        "this number is added to the average error per column to get the target tolerance",
        io::Serialization::GetAgent( &m_ToleranceMinusError),
        "0"
      );
      parameters.AddInitializer
      (
        "reduce error by tolerance",
        "Whether to reduce the error for backpropagated features by the tolerance when the feature is backpropagated",
        io::Serialization::GetAgent( &m_ReduceErrorByTolerance),
        "True"
      );
      parameters.AddInitializer
      (
        "pure classification",
        "True if result values should be pretended to be binary, depending on their side of the cutoff",
        io::Serialization::GetAgent( &m_PureClassification),
        "False"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
