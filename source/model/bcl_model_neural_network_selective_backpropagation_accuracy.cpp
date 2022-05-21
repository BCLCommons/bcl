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
#include "model/bcl_model_neural_network_selective_backpropagation_accuracy.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkSelectiveBackpropagationAccuracy::s_Instance
    (
      util::Enumerated< NeuralNetworkSelectiveBackpropagationInterface>::AddInstance
      (
        new NeuralNetworkSelectiveBackpropagationAccuracy()
      )
    );

    //! @brief default constructor
    NeuralNetworkSelectiveBackpropagationAccuracy::NeuralNetworkSelectiveBackpropagationAccuracy() :
      m_ResultStartID( 0),
      m_ResultEndID( 1000),
      m_IgnoreTruePredictionsBelow( 0.1),
      m_IgnoreTruePredictionsAbove( 0.9),
      m_IgnoreFalsePredictionsBelow ( 0.0),
      m_IgnoreFalsePredictionsAbove ( 1.0),
      m_PureClassification( false),
      m_BalanceError( false),
      m_FlatErrorFunction( false),
      m_ResultsSize( 1000),
      m_NumberThreads( 0)
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkSelectiveBackpropagationAccuracy copied from this instance
    NeuralNetworkSelectiveBackpropagationAccuracy *NeuralNetworkSelectiveBackpropagationAccuracy::Clone() const
    {
      return new NeuralNetworkSelectiveBackpropagationAccuracy( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkSelectiveBackpropagationAccuracy::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkSelectiveBackpropagationAccuracy::GetAlias() const
    {
      static const std::string s_name( "Accuracy");
      return s_name;
    }

    //! @brief Initialize this object with the rescaled dataset and other important parameters
    //! @param RESCALED_DATA the already-rescaled training dataset
    //! @param OBJECTIVE the actual objective function
    //! @param THREADS # of threads that may be accessing this object at once
    void NeuralNetworkSelectiveBackpropagationAccuracy::Initialize
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

      // set up sizes of all internally-held members
      m_ScaledCutoffs = linal::Vector< float>( m_ResultsSize);

      m_IncorrectAboveCutoffTally.Resize( m_NumberThreads);
      m_IncorrectAboveCutoffTally.SetAllElements( linal::Vector< size_t>( m_ResultsSize, size_t( 0)));
      m_IncorrectBelowCutoffTally = m_CorrectAboveCutoffTally = m_CorrectBelowCutoffTally
          = m_IncorrectAboveCutoffTally;

      // regression task: Compute standard deviation over the selected columns
      const float cutoff( OBJECTIVE.GetThreshold());
      m_ScaledCutoffs = linal::Vector< float>( m_ResultsSize, float( 0.0));
      m_Results = RESCALED_DATA.GetResultsPtr();

      const util::SiPtr< const RescaleFeatureDataSet> results_scaling( RESCALED_DATA.GetResultsPtr()->GetScaling());
      for( size_t result( 0); result < m_ResultsSize; ++result)
      {
        m_ScaledCutoffs( result) = results_scaling->RescaleValue( result, cutoff);
      }
      m_LastRoundInaccuraciesAboveCutoff = m_LastRoundInaccuraciesBelowCutoff
        = linal::Vector< float>( m_ResultsSize, float( 0.0));
      m_LastRoundRelativeErrorAboveCutoff = m_LastRoundRelativeErrorBelowCutoff
        = linal::Vector< float>( m_ResultsSize, float( 1.0));
      m_RoundNumber = 1;
    } // Initialize

    //! @brief select whether or not to backpropagate the current feature, and edit the error vector, if desired
    //! @param PREDICTION the prediction that was made by the neural network
    //! @param ERROR Reference to the error vector, (should already have been set to RESULT - PREDICTION)
    //! @param FEATURE_ID id of this feature in the dataset
    //! @param THREAD_ID id of this thread
    //! @return true if the feature should be backpropagated
    bool NeuralNetworkSelectiveBackpropagationAccuracy::ShouldBackpropagate
    (
      const linal::VectorConstInterface< float> &PREDICTION,
      linal::VectorInterface< float> &ERROR,
      const size_t &FEATURE_ID,
      const size_t &THREAD_ID
    )
    {
      bool backprop_should_continue( false);
      linal::VectorConstReference< float> result_row
      (
        m_ResultsSize,
        m_Results->GetMatrix()[ FEATURE_ID] + m_ResultStartID
      );
      for( size_t result_offset( 0); result_offset < m_ResultsSize; ++result_offset)
      {
        // get the actual index
        const size_t result( result_offset + m_ResultStartID);

        const float experimental( result_row( result_offset));
        const bool exp_was_above_cutoff( experimental > m_ScaledCutoffs( result_offset));
        const float prediction( PREDICTION( result));
        const float cutoff
        (
          m_PureClassification
          ? ( exp_was_above_cutoff ? std::min( float( 0.5), m_IgnoreFalsePredictionsAbove) : std::max( float( 0.5), m_IgnoreFalsePredictionsBelow))
          : m_ScaledCutoffs( result_offset)
        );
        const bool predicted_was_above_cutoff( prediction > cutoff);

        float &er( ERROR( result));

        const float inaccuracy_last_round
        (
          ( exp_was_above_cutoff ? m_LastRoundInaccuraciesAboveCutoff : m_LastRoundInaccuraciesBelowCutoff)( result_offset)
        );

        float distance_to_saturation( 0.0), cutoff_distance_to_saturation( 0);

        // tally correct / incorrect counts based on cutoff side
        if( exp_was_above_cutoff == predicted_was_above_cutoff)
        {
          if( predicted_was_above_cutoff)
          {
            ++m_CorrectAboveCutoffTally( THREAD_ID)( result_offset);
            distance_to_saturation = m_IgnoreTruePredictionsAbove - prediction;
            cutoff_distance_to_saturation = m_IgnoreTruePredictionsAbove - cutoff;
          }
          else
          {
            ++m_CorrectBelowCutoffTally( THREAD_ID)( result_offset);
            distance_to_saturation = prediction - m_IgnoreTruePredictionsBelow;
            cutoff_distance_to_saturation = cutoff - m_IgnoreTruePredictionsBelow;
          }
        }
        else
        {
          if( predicted_was_above_cutoff)
          {
            ++m_IncorrectBelowCutoffTally( THREAD_ID)( result_offset);
            distance_to_saturation = m_IgnoreFalsePredictionsAbove - prediction;
          }
          else
          {
            ++m_IncorrectAboveCutoffTally( THREAD_ID)( result_offset);
            distance_to_saturation = prediction - m_IgnoreFalsePredictionsBelow;
          }
        }

        const bool last_round_was_mostly_right( inaccuracy_last_round < 0.5);
        const bool on_correct_side( exp_was_above_cutoff == predicted_was_above_cutoff);
        if( ( on_correct_side || last_round_was_mostly_right) && distance_to_saturation <= 0.0)
        {
          // prediction beyond saturation range, and in the correct direction, or at least half the values were in the
          // correct direction last time, so assume features that yield saturated values cannot reasonably be predicted
          // from the descriptors

          // saturated result; suppose results are scaled to 0.1 - 0.9,
          // so if the predicted value was on the correct side of
          // the cutoff and is either between 0.0 - 0.1 or 0.9 - 1.0, then we're in the saturation range of the transfer
          // function and it results in an ill-conditioned network in this case to backpropagate an error term that
          // tries to pin the result exactly at 0.1 or 0.9.
          er = 0.0;
          continue;
        }
        if( er == 0.0)
        {
          continue;
        }
        backprop_should_continue = true;
        if( m_FlatErrorFunction || !on_correct_side)
        {
          er = exp_was_above_cutoff ? 1.0 : -1.0;
        }
        else if( m_PureClassification)
        {
          er = ( er > 0.0 ? 1.0 : -1.0) * math::Absolute( distance_to_saturation) / cutoff_distance_to_saturation;
        }
        if( m_BalanceError)
        {
          er *= ( exp_was_above_cutoff ? m_LastRoundRelativeErrorAboveCutoff : m_LastRoundRelativeErrorBelowCutoff)( result_offset);
        }
      }

      // return true if any part of ERROR remains to be backpropagated
      return backprop_should_continue;
    }

    //! @brief finalize the current round; occurs only after all threads were already joined
    void NeuralNetworkSelectiveBackpropagationAccuracy::FinalizeRound()
    {
      // accumulate results
      for( size_t thread_number( 1); thread_number < m_NumberThreads; ++thread_number)
      {
        m_IncorrectAboveCutoffTally( 0) += m_IncorrectAboveCutoffTally( thread_number);
        m_IncorrectBelowCutoffTally( 0) += m_IncorrectBelowCutoffTally( thread_number);
        m_CorrectAboveCutoffTally( 0)   += m_CorrectAboveCutoffTally( thread_number);
        m_CorrectBelowCutoffTally( 0)   += m_CorrectBelowCutoffTally( thread_number);
      }
      for( size_t i( 0); i < m_ResultsSize; ++i)
      {
        m_LastRoundInaccuraciesAboveCutoff( i) = float( m_IncorrectAboveCutoffTally( 0)( i)) / float( m_IncorrectAboveCutoffTally( 0)( i) + m_CorrectAboveCutoffTally( 0)( i));
        m_LastRoundInaccuraciesBelowCutoff( i) = float( m_IncorrectBelowCutoffTally( 0)( i)) / float( m_IncorrectBelowCutoffTally( 0)( i) + m_CorrectBelowCutoffTally( 0)( i));
        if( m_BalanceError)
        {
          if( m_LastRoundInaccuraciesAboveCutoff( i) > m_LastRoundInaccuraciesBelowCutoff( i))
          {
            m_LastRoundRelativeErrorAboveCutoff( i) *= 1.05;
            m_LastRoundRelativeErrorBelowCutoff( i) /= 1.05;
          }
          else
          {
            m_LastRoundRelativeErrorAboveCutoff( i) /= 1.05;
            m_LastRoundRelativeErrorBelowCutoff( i) *= 1.05;
          }
        }
        BCL_MessageStd
        (
          "Round #" + util::Format()( m_RoundNumber) + " inaccuracy below cutoff for result " + util::Format()( i)
          + " = " + util::Format()( m_LastRoundInaccuraciesBelowCutoff( i))
          + " above cutoff = " + util::Format()( m_LastRoundInaccuraciesAboveCutoff( i)) +
          (
            m_BalanceError
            ? " error ratio above cutoff: " + util::Format()( m_LastRoundRelativeErrorAboveCutoff( i))
              + " error ratio below cutoff: " + util::Format()( m_LastRoundRelativeErrorBelowCutoff( i))
            : std::string()
          )
        );
      }
      for( size_t thread_number( 0); thread_number < m_NumberThreads; ++thread_number)
      {
        m_IncorrectAboveCutoffTally( thread_number) = size_t( 0);
        m_IncorrectBelowCutoffTally( thread_number) = size_t( 0);
        m_CorrectAboveCutoffTally( thread_number) = size_t( 0);
        m_CorrectBelowCutoffTally( thread_number) = size_t( 0);
      }
      ++m_RoundNumber;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkSelectiveBackpropagationAccuracy::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Always back-propagating features that are mispredicted on either side of the cutoff. "
        "Features predicted on the currect side of the cutoff, but closer to the cutoff than the actual result, are "
        "backpropagated with a user-specified probability"
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
        "pure classification",
        "True if result values should be pretended to be binary, depending on their side of the cutoff",
        io::Serialization::GetAgent( &m_PureClassification),
        "False"
      );
      parameters.AddInitializer
      (
        "balance error",
        "True to weight error values so as to try to balance inaccuracies on each side of the cutoff. This is done via"
        "the use of multipliers that are adjusted in an iterative fashion",
        io::Serialization::GetAgent( &m_BalanceError),
        "False"
      );
      parameters.AddInitializer
      (
        "tolerance above",
        "Rescaled boundary above which accurate predictions are not backpropagated",
        io::Serialization::GetAgent( &m_IgnoreTruePredictionsAbove),
        "0.9"
      );
      parameters.AddInitializer
      (
        "tolerance below",
        "Rescaled boundary below which accurate predictions are not backpropagated",
        io::Serialization::GetAgent( &m_IgnoreTruePredictionsBelow),
        "0.1"
      );
      parameters.AddInitializer
      (
        "ignore false predictions below",
        "Do not backpropagate false predictions that are below this value",
        io::Serialization::GetAgent( &m_IgnoreFalsePredictionsBelow),
        "0.0"
      );
      parameters.AddInitializer
      (
        "ignore false predictions above",
        "Ddo not backpropagate false predictions that are above this value",
        io::Serialization::GetAgent( &m_IgnoreFalsePredictionsAbove),
        "1.0"
      );
      parameters.AddInitializer
      (
        "flat error function",
        "If true, all error values outside the saturation region are set to +/- 1.",
        io::Serialization::GetAgent( &m_FlatErrorFunction),
        "False"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
