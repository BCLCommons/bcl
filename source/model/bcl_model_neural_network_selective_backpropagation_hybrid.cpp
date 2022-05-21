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
#include "model/bcl_model_neural_network_selective_backpropagation_hybrid.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkSelectiveBackpropagationHybrid::s_Instance
    (
      util::Enumerated< NeuralNetworkSelectiveBackpropagationInterface>::AddInstance
      (
        new NeuralNetworkSelectiveBackpropagationHybrid()
      )
    );

    //! @brief default constructor
    NeuralNetworkSelectiveBackpropagationHybrid::NeuralNetworkSelectiveBackpropagationHybrid() :
      m_PureClassification( false),
      m_ResultStartID( 0),
      m_ResultEndID( 1000),
      m_EnrichmentCutoff( 0.0),
      m_StabilityRatio( 2),
      m_SensitivityWeight( 0.5),
      m_FalsePositiveBonus( 0),
      m_FalseNegativeBonus( 0),
      m_StartWithAll( false),
      m_ResultsSize( 1000),
      m_DatasetSize( 1),
      m_NumberThreads( 0),
      m_Parity( true),
      m_RoundNumber( 0)
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkSelectiveBackpropagationHybrid copied from this instance
    NeuralNetworkSelectiveBackpropagationHybrid *NeuralNetworkSelectiveBackpropagationHybrid::Clone() const
    {
      return new NeuralNetworkSelectiveBackpropagationHybrid( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkSelectiveBackpropagationHybrid::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkSelectiveBackpropagationHybrid::GetAlias() const
    {
      static const std::string s_name( "Hybrid");
      return s_name;
    }

    //! @brief Initialize this object with the rescaled dataset and other important parameters
    //! @param RESCALED_DATA the already-rescaled training dataset
    //! @param OBJECTIVE the actual objective function
    //! @param THREADS # of threads that may be accessing this object at once
    void NeuralNetworkSelectiveBackpropagationHybrid::Initialize
    (
      const descriptor::Dataset &RESCALED_DATA,
      const ObjectiveFunctionInterface &OBJECTIVE,
      const size_t &NUMBER_THREADS
    )
    {
      // initialization

      // bound m_ResultEndID by the results size
      if( m_ResultEndID > RESCALED_DATA.GetResultSize())
      {
        m_ResultEndID = RESCALED_DATA.GetResultSize();
      }
      // compute internal results size
      m_ResultsSize = m_ResultEndID - m_ResultStartID;

      m_RoundNumber = 1;

      // store required data
      m_NumberThreads = NUMBER_THREADS;
      m_Results = RESCALED_DATA.GetResultsPtr()->GetMatrix();
      m_PredictedValue = m_Results;
      m_DatasetSize = RESCALED_DATA.GetResultsPtr()->GetNumberFeatures();
      m_Parity = OBJECTIVE.GetRankingParity();

      // set up sizes of all internally-held members
      m_ScaledCutoffs = linal::Vector< float>( m_ResultsSize);
      m_MaxDelta = m_ScaledCutoffs;
      m_IncorrectTallyAboveCutoff.Resize( m_NumberThreads);
      m_IncorrectTallyAboveCutoff.SetAllElements( linal::Vector< size_t>( m_ResultsSize, size_t( 0)));
      m_IncorrectTallyBelowCutoff = m_IncorrectTallyAboveCutoff;
      m_BackpropagatedTallyAboveCutoff = m_BackpropagatedTallyBelowCutoff = m_IncorrectTallyAboveCutoff;
      m_BelowCutoffCount = m_AboveCutoffCount = linal::Vector< size_t>( m_ResultsSize, size_t( 0));
      m_BackPropPreference = linal::Matrix< float>( m_DatasetSize, m_ResultsSize, float( 1.0));
      m_ChangeLastRound = linal::Matrix< float>( m_DatasetSize, m_ResultsSize, float( 1.0));
      m_AverageChange = linal::Matrix< float>( m_DatasetSize, m_ResultsSize, float( 0.0));
      m_ThreadMaxDelta.Resize( m_NumberThreads);
      m_ThreadMaxDelta.SetAllElements( storage::Vector< float>( m_ResultsSize, float( 0.0)));
      m_ThresholdBPHigh = m_ThresholdBPLow = m_SensitivityScalingHigh = m_SensitivityScalingLow = linal::Vector< float>( m_ResultsSize, 0.5);

      // get scaling and cutoff information
      const float cutoff( OBJECTIVE.GetThreshold());
      const util::SiPtr< const RescaleFeatureDataSet> results_scaling( RESCALED_DATA.GetResultsPtr()->GetScaling());

      // take data about the dataset
      for( size_t result( 0); result < m_ResultsSize; ++result)
      {
        // get the scaled cutoff
        m_ScaledCutoffs( result) = results_scaling->RescaleValue( result, cutoff);

        // compute average distance from cutoff above and below the cutoff
        math::RunningAverage< float> ave_high, ave_low;
        float max_delta_high( 0), max_delta_low( 0);

        // count the number of times results are seen above cutoff, also find the max/min of the dataset
        for( size_t i( 0); i < m_DatasetSize; ++i)
        {
          const float result_value( m_Results( i, result));
          if( result_value >= m_ScaledCutoffs( result))
          {
            ++m_AboveCutoffCount( result);
            max_delta_high = std::max( result_value, max_delta_high);
          }
          else
          {
            ++m_BelowCutoffCount( result);
            max_delta_low = std::min( result_value, max_delta_low);
          }
        }

        // Initialize the max delta with the min/max results seen in the dataset
        if( m_AboveCutoffCount( result))
        {
          m_MaxDelta( result) = max_delta_high - m_ScaledCutoffs( result);
        }
        if( m_BelowCutoffCount( result))
        {
          m_MaxDelta( result) = std::max( m_MaxDelta( result), m_ScaledCutoffs( result) - max_delta_low);
        }

        // Print off counts above and below the cutoff
        BCL_MessageStd
        (
          "For result # " + util::Format()( result + m_ResultStartID)
          + " #Above cutoff: " + util::Format()( m_AboveCutoffCount( result))
          + " #Below: " + util::Format()( m_BelowCutoffCount( result))
        );

        // determine initial backpropagation probabilities so as to balance the number of positives and negatives
        // backpropagated
        float bp_prob_high( 1.0);
        float bp_prob_low( 1.0);

        // handle datasets that have overrepresentation factors greater than 3/2.  Datasets with smaller
        // overrepresentation factors typically do not require balancing on the first round
        if( m_AboveCutoffCount( result) > m_BelowCutoffCount( result) * 3 / 2)
        {
          bp_prob_high
            = std::min
              (
                float( 1.0),
                float( m_StabilityRatio * m_BelowCutoffCount( result))
                / float( std::max( m_AboveCutoffCount( result), size_t( 1)))
              );
        }
        else if( m_AboveCutoffCount( result) * 3 / 2 < m_BelowCutoffCount( result))
        {
          bp_prob_low
            = std::min
              (
                float( 1.0),
                float( m_StabilityRatio * m_AboveCutoffCount( result))
                / float( std::max( m_BelowCutoffCount( result), size_t( 1)))
              );
        }

        // set initial backpropagation propensities deterministically, in accordance with the backpropagation
        // probabilities necessary to balance the solution.  The first instance of the class is always backpropagated
        float current_bp_low( 0.0), current_bp_high( 0.0);
        for( size_t i( 0); i < m_DatasetSize; ++i)
        {
          if( m_Results( i, result) >= m_ScaledCutoffs( result))
          {
            // positives
            current_bp_high += bp_prob_high;
            m_BackPropPreference( i, result) = current_bp_high;
            if( current_bp_high >= 0.0)
            {
              current_bp_high -= 1.0;
            }
          }
          else
          {
            // negatives
            current_bp_low += bp_prob_low;
            m_BackPropPreference( i, result) = current_bp_low;
            if( current_bp_low >= 0.0)
            {
              current_bp_low -= 1.0;
            }
          }
        }

        // avoid numerical issues that arise if there was only one class represented in the data
        if( !m_AboveCutoffCount( result))
        {
          ++m_AboveCutoffCount( result);
          --m_BelowCutoffCount( result);
        }
        else if( !m_BelowCutoffCount( result))
        {
          --m_AboveCutoffCount( result);
          ++m_BelowCutoffCount( result);
        }
      }
      m_EffectiveCutoffs = m_ScaledCutoffs;
      if( m_PureClassification)
      {
        m_EffectiveCutoffs = float( 0.5);
      }
    } // Initialize

    //! @brief select whether or not to backpropagate the current feature, and edit the error vector, if desired
    //! @param PREDICTION the prediction that was made by the neural network
    //! @param ERROR Reference to the error vector, (should already have been set to RESULT - PREDICTION)
    //! @param FEATURE_ID id of this feature in the dataset
    //! @param THREAD_ID id of this thread
    //! @return true if the feature should be backpropagated
    bool NeuralNetworkSelectiveBackpropagationHybrid::ShouldBackpropagate
    (
      const linal::VectorConstInterface< float> &PREDICTION,
      linal::VectorInterface< float> &ERROR,
      const size_t &FEATURE_ID,
      const size_t &THREAD_ID
    )
    {
      // get a reference to the result and prediction row, already offset by m_ResultStartID
      linal::VectorConstReference< float> result_row( m_ResultsSize, m_Results[ FEATURE_ID] + m_ResultStartID);
      linal::VectorConstReference< float> prediction_row( m_ResultsSize, PREDICTION.Begin() + m_ResultStartID);

      // get references to the tallying arrays for this thread
      linal::Vector< size_t> &incorrect_above_tally( m_IncorrectTallyAboveCutoff( THREAD_ID));
      linal::Vector< size_t> &incorrect_below_tally( m_IncorrectTallyBelowCutoff( THREAD_ID));
      linal::Vector< size_t> &bp_above_tally( m_BackpropagatedTallyAboveCutoff( THREAD_ID));
      linal::Vector< size_t> &bp_below_tally( m_BackpropagatedTallyBelowCutoff( THREAD_ID));
      storage::Vector< float> &max_delta( m_ThreadMaxDelta( THREAD_ID));

      bool backprop_should_continue( false);

      // for each result
      for( size_t result_id( 0); result_id < m_ResultsSize; ++result_id)
      {
        // get the actual result, prediction, and cutoff
        const float actual_result_value( result_row( result_id));
        const float prediction( prediction_row( result_id));
        const float cutoff( m_ScaledCutoffs( result_id));
        const float effective_cutoff( m_EffectiveCutoffs( result_id));
        float &last_round_prediction( m_PredictedValue( FEATURE_ID, result_id));

        // determine whether the value was above cutoff
        const bool is_high( actual_result_value >= cutoff);

        // compute the difference between last rounds prediction and this rounds'
        const float delta_prediction( math::Absolute( prediction - last_round_prediction));
        const float delta_cutoff( math::Absolute( prediction - effective_cutoff));

        // keep track of how much the prediction has changed over the last round
        m_ChangeLastRound( FEATURE_ID, result_id) = delta_prediction;

        // update the threads maximum delta by this point
        max_delta( result_id) = std::max( max_delta( result_id), delta_prediction);

        // determine whether this example was predicted incorrectly
        const bool is_predicted_wrong( is_high != ( prediction >= effective_cutoff));

        // determine whether this exemplar must be backpropagated
        //bool backprop_this_column( m_BackPropPreference( FEATURE_ID, result_id) >= 0.0);
        bool backprop_this_column( m_RoundNumber == 1 && m_BackPropPreference( FEATURE_ID, result_id) >= 0);
        if( is_high)
        {
          const float bonus( is_predicted_wrong && m_Parity ? m_FalseNegativeBonus : 0.0);
          const float sensitivity_addition
          (
            m_SensitivityWeight * m_SensitivityScalingHigh( result_id)
            * ( m_AverageChange( FEATURE_ID, result_id) + delta_prediction) / 2.0
          );
          const float importance( delta_cutoff - bonus - sensitivity_addition);
          m_BackPropPreference( FEATURE_ID, result_id) = importance;
          if( m_RoundNumber > 1 || m_StartWithAll)
          {
            backprop_this_column = importance <= m_ThresholdBPHigh( result_id);
          }
        }
        else
        {
          const float bonus( is_predicted_wrong && !m_Parity ? m_FalsePositiveBonus : 0.0);
          const float sensitivity_addition
          (
            m_SensitivityWeight * m_SensitivityScalingLow( result_id)
            * ( m_AverageChange( FEATURE_ID, result_id) + delta_prediction) / 2.0
          );
          const float importance( delta_cutoff - bonus - sensitivity_addition);
          m_BackPropPreference( FEATURE_ID, result_id) = importance;
          if( m_RoundNumber > 1 || m_StartWithAll)
          {
            backprop_this_column = importance <= m_ThresholdBPLow( result_id);
          }
        }

        // determine whether the predicted value is on the correct side of the cutoff
        // note that for pure classification tasks, a cutoff of 0.5 is imposed.  This is only the proper behavior though
        // if the objective function is a rank classifier though
        if( is_predicted_wrong)
        {
          // update incorrect tally
          ++( is_high ? incorrect_above_tally : incorrect_below_tally)( result_id);
        }

        // update last round prediction
        last_round_prediction = prediction;

        // handle the case that this column will be backpropagated
        if( backprop_this_column)
        {
          // feature must be backpropagated
          backprop_should_continue = true;

          // update backpropagated tally
          ++( is_high ? bp_above_tally : bp_below_tally)( result_id);

          // modify the backpropagated error for pure classification tasks so that all exemplars have a result of 0 or 1
          if( m_PureClassification)
          {
            ERROR( result_id + m_ResultStartID) = ( is_high ? 1.0 : 0.0) - prediction;
          }
        }
        else
        {
          // this column is not to be backpropagated; so set the error to 0 (since it is possible other columns may be
          // backpropagated, but it has already been decided not to backpropagate this column)
          ERROR( result_id + m_ResultStartID) = 0.0;
        }
      }

      // return true if any part of ERROR remains to be backpropagated
      return backprop_should_continue;
    }

    //! @brief finalize the current round; occurs only after all threads were already joined
    void NeuralNetworkSelectiveBackpropagationHybrid::FinalizeRound()
    {
      // only keep the average change over the last 10 rounds ( 10 * 9 / 2 = 45)
      if( m_ChangeAverager.GetWeight() > 10.0)
      {
        m_ChangeAverager.SetWeight( 9.0);
      }

      // update average change
      m_ChangeAverager += m_ChangeLastRound;

      // increment round number
      m_RoundNumber += 1;

      // accumulate results
      for( size_t thread_number( 1); thread_number < m_NumberThreads; ++thread_number)
      {
        m_IncorrectTallyAboveCutoff( 0) += m_IncorrectTallyAboveCutoff( thread_number);
        m_IncorrectTallyBelowCutoff( 0) += m_IncorrectTallyBelowCutoff( thread_number);
        m_BackpropagatedTallyAboveCutoff( 0) += m_BackpropagatedTallyAboveCutoff( thread_number);
        m_BackpropagatedTallyBelowCutoff( 0) += m_BackpropagatedTallyBelowCutoff( thread_number);
        for( size_t i( 0); i < m_ResultsSize; ++i)
        {
          m_ThreadMaxDelta( 0)( i) = std::max( m_ThreadMaxDelta( 0)( i), m_ThreadMaxDelta( thread_number)( i));
        }
      }

      const linal::Matrix< float> &avg_change( m_ChangeAverager.GetAverage());
      m_AverageChange = avg_change;
      for( size_t i( 0); i < m_ResultsSize; ++i)
      {
        // store the max deltas
        m_MaxDelta( i) = m_ThreadMaxDelta( 0)( i);

        // compute the previous round's inaccuracy above / below cutoff, to display for the user
        const float inaccuracy_above( float( m_IncorrectTallyAboveCutoff( 0)( i)) / float( m_AboveCutoffCount( i)));
        const float inaccuracy_below( float( m_IncorrectTallyBelowCutoff( 0)( i)) / float( m_BelowCutoffCount( i)));

        // compute the absolute minimum and maximum number of exemplars that should be backpropagated based on the user's
        // given cutoffs
        const float max_above( m_AboveCutoffCount( i));
        const float max_below( m_BelowCutoffCount( i));

        // save the # inaccurate above and below the cutoff last turn
        const float n_inaccurate_above( m_IncorrectTallyAboveCutoff( 0)( i));
        const float n_inaccurate_below( m_IncorrectTallyBelowCutoff( 0)( i));

        // expected_high is the number of features that must be backpropagated next turn that are above the cutoff
        // it is calculated as 1/2 the number that were backpropagated last turn + the optimal # that should be calculated
        // to optimize enrichment (if enrichment cutoff was set) or to ensure balancing
        // Carrying 1/2 the number that were backpropagated last turn forward stabilizes the calculation against, e.g.
        // predicting all the positives correct one turn, then all the negatives correct the following turn
        float expected_high( 0.75 * m_BackpropagatedTallyAboveCutoff( 0)( i));
        float expected_low( 0.75 * m_BackpropagatedTallyBelowCutoff( 0)( i));
        if( m_EnrichmentCutoff)
        {
          // enrichment cutoff set; optimize for enrichment

          // compute the current contingency matrix for the training data using the experimental cutoff
          float false_positives( 0), false_negatives( 0);

          // Enrichment at a constant cutoff (e.g. TP+FP = Constant) = TP/C = 1-FP/C
          // taking the derivative with respect to TP and FP, it should be clear that the change in enrichment is identical
          // in value (but opposite in sign) for exchanging a TP for an FP.  Given that we currently have TP=x, and FP=y,
          // there are FN chances to increase enrichment by backpropagating a positive, and min(FN,FP) chances
          // to increase the value by BPing an FP.
          float chances_to_increase_above( 0), chances_to_increase_below( 0);
          if( m_Parity)
          {
            false_positives = n_inaccurate_below;
            false_negatives = n_inaccurate_above;
            chances_to_increase_above = false_negatives + 1.0;
            chances_to_increase_below = std::min( false_negatives, false_positives) + 1.0;
            if( m_EnrichmentCutoff * m_DatasetSize > max_above && ( false_positives + false_negatives) < m_EnrichmentCutoff * m_DatasetSize)
            {
              chances_to_increase_below = 1;
            }
          }
          else
          {
            false_positives = n_inaccurate_above;
            false_negatives = n_inaccurate_below;
            chances_to_increase_above = std::min( false_negatives, false_positives) + 1.0;
            chances_to_increase_below = false_negatives + 1.0;
            if( m_EnrichmentCutoff * m_DatasetSize > max_below && ( false_positives + false_negatives) < m_EnrichmentCutoff * m_DatasetSize)
            {
              chances_to_increase_above = 1;
            }
          }

          // if there are more incorrect results than the user-specified enrichment cutoff would allow,
          // then do not optimize enrichment,
          // otherwise this leads to an instability where the optimum is reached by predicting everything positive
          // Instead, just backpropagate some factor of the number inaccurate on either side of the cutoff
          if( m_EnrichmentCutoff * m_DatasetSize < ( false_negatives + false_positives))
          {
            chances_to_increase_above = n_inaccurate_above * m_StabilityRatio;
            chances_to_increase_below = n_inaccurate_below * m_StabilityRatio;
          }
          else
          {
            chances_to_increase_above = std::max( chances_to_increase_above, float( 1));
            chances_to_increase_below = std::max( chances_to_increase_below, float( 1));

            // compute the effective stability ratio
            const float effective_stability_ratio
            (
              std::min
              (
                std::min( m_StabilityRatio, max_above / chances_to_increase_above),
                max_below / chances_to_increase_below
              )
            );

            // optimize directly for enrichment
            chances_to_increase_above *= effective_stability_ratio;
            chances_to_increase_below *= effective_stability_ratio;
          }

          // bound by the user-specified limits
          chances_to_increase_above = std::min( max_above, std::max( chances_to_increase_above, float( 1)));
          chances_to_increase_below = std::min( max_below, std::max( chances_to_increase_below, float( 1)));

          // add the 1/4 the optimal value to the expected value to provide stability
          expected_high += 0.25 * chances_to_increase_above;
          expected_low += 0.25 * chances_to_increase_below;
        }
        else
        {
          float desired_high_steady_state( n_inaccurate_above * m_StabilityRatio);
          float desired_low_steady_state( n_inaccurate_below * m_StabilityRatio);

          // bound by the user-specified limits
          desired_high_steady_state = std::min( max_above, std::max( desired_high_steady_state, float( 0)));
          desired_low_steady_state = std::min( max_below, std::max( desired_low_steady_state, float( 0)));

          expected_high += 0.25 * desired_high_steady_state;
          expected_low += 0.25 * desired_low_steady_state;
        }

        // sort the next round's backpropagation preferences above and below the cutoff
        // Compute m_BackPropPreference( j, i) -= x with x such that the number of exemplars with
        // m_BackPropPreference( j, i) >= 0 is the desired number of exemplars for that class (expected above or below)

        // also compute the average distance to cutoff above and below cutoff and the average change above and below cutoff
        // to compute the sensitivity scaling

        math::RunningAverage< float> ave_dist_to_cutoff_high, ave_dist_to_cutoff_low;
        math::RunningAverage< float> ave_change_high, ave_change_low;

        storage::Vector< float> high_values, low_values;
        high_values.AllocateMemory( m_AboveCutoffCount( i));
        low_values.AllocateMemory( m_BelowCutoffCount( i));
        for( size_t j( 0); j < m_DatasetSize; ++j)
        {
          const float result( m_Results( j, i));
          const float dist_to_cutoff( math::Absolute( m_PredictedValue( j, i) - m_EffectiveCutoffs( i)));
          if( result >= m_ScaledCutoffs( i))
          {
            high_values.PushBack( m_BackPropPreference( j, i));
            if( m_BackPropPreference( j, i) <= m_ThresholdBPHigh( i))
            {
              ave_dist_to_cutoff_high += dist_to_cutoff;
              ave_change_high += m_AverageChange( j, i);
            }
          }
          else
          {
            low_values.PushBack( m_BackPropPreference( j, i));
            if( m_BackPropPreference( j, i) <= m_ThresholdBPLow( i))
            {
              ave_dist_to_cutoff_low += dist_to_cutoff;
              ave_change_low += m_AverageChange( j, i);
            }
          }
        }

        m_SensitivityScalingHigh( i) = std::max( float( 0.001), ave_dist_to_cutoff_high.GetAverage()) / std::max( float( 0.001), ave_change_high.GetAverage());
        m_SensitivityScalingLow( i) = std::max( float( 0.001), ave_dist_to_cutoff_low.GetAverage()) / std::max( float( 0.001), ave_change_low.GetAverage());

        // sort the backpropagation preferences above and below this results' cutoff
        high_values.Sort( std::less< float>());
        low_values.Sort( std::less< float>());

        // choose the cutoff for exemplars below cutoff such that there are expected_low number of values >= 0
        const float cutoff_low
        (
          low_values.IsEmpty()
          ? float( 0.5)
          : low_values( std::min( std::max( size_t( expected_low + 0.5), size_t( 1)) - 1, low_values.GetSize() - 1))
        );

        // choose the cutoff for exemplars abpve cutoff such that there are expected_high number of values >= 0
        const float cutoff_high
        (
          high_values.IsEmpty()
          ? float( 0.5)
          : high_values( std::min( std::max( size_t( expected_high + 0.5), size_t( 1)) - 1, high_values.GetSize() - 1))
        );

        m_ThresholdBPHigh( i) = cutoff_high;
        m_ThresholdBPLow( i) = cutoff_low;

        BCL_MessageVrb
        (
          "Last round inaccuracy on result # " + util::Format()( i + m_ResultStartID)
          + " above cutoff = " + util::Format()( inaccuracy_above)
          + " below cutoff = " + util::Format()( inaccuracy_below)
          + " # Incorrect above/below: " + util::Format()( n_inaccurate_above)
          + " " + util::Format()( n_inaccurate_below)
          + " # BP above / below: " + util::Format()( m_BackpropagatedTallyAboveCutoff( 0)( i))
          + " " + util::Format()( m_BackpropagatedTallyBelowCutoff( 0)( i))
          + " Cutoff Low: " + util::Format()( m_ThresholdBPLow( i))
          + " Cutoff high: " + util::Format()( m_ThresholdBPHigh( i))
        );
      }

      // reset thread-specific arrays
      for( size_t thread_number( 0); thread_number < m_NumberThreads; ++thread_number)
      {
        m_IncorrectTallyAboveCutoff( thread_number) = size_t( 0);
        m_IncorrectTallyBelowCutoff( thread_number) = size_t( 0);
        m_BackpropagatedTallyAboveCutoff( thread_number) = size_t( 0);
        m_BackpropagatedTallyBelowCutoff( thread_number) = size_t( 0);
        m_ThreadMaxDelta( thread_number).SetAllElements( 0.0);
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkSelectiveBackpropagationHybrid::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Backpropagates preferrentially according to sensitivity, distance from cutoff, "
        "classification status (correct/incorrect), and cutoff side. "
        "See full explanation and advice about parameters at "
        "https://structbio.vanderbilt.edu:8443/display/MeilerLab/BclANNParameterSelection"
      );
      parameters.AddInitializer
      (
        "begin",
        "result columns to consider",
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
        "enrichment max",
        "flag to deliberately optimize for enrichment with the given cutoff, leave 0 to focus on balancing rather than enrichment",
        io::Serialization::GetAgentWithRange( &m_EnrichmentCutoff, float( 0.0), float( 1.0)),
        "0.0"
      );
      parameters.AddInitializer
      (
        "pure classification",
        "flag; set to ignore the result values when backpropagating, instead, backprop 0 or 1 (for above cutoff)",
        io::Serialization::GetAgent( &m_PureClassification),
        "False"
      );
      parameters.AddInitializer
      (
        "stability",
        "Ratio of FP+FN backpropagated due to incorrectness to those backpropped due to rapid change. "
        "Use 1 if the data are very clean and representative of the overall space; larger if not",
        io::Serialization::GetAgent( &m_StabilityRatio),
        "2"
      );
      parameters.AddInitializer
      (
        "sensitivity weight",
        "Weighting of sensitivity. The only other weight term (which is 1-sensitivity weight) accounts for outliers and"
        " false predictions about the cutoff",
        io::Serialization::GetAgent( &m_SensitivityWeight),
        "0.5"
      );
      parameters.AddInitializer
      (
        "fn bonus",
        "Give a bonus to all false negatives BP tendancy",
        io::Serialization::GetAgent( &m_FalseNegativeBonus),
        "0"
      );
      parameters.AddInitializer
      (
        "fp bonus",
        "Give a bonus to all false positives BP tendancy",
        io::Serialization::GetAgent( &m_FalsePositiveBonus),
        "0"
      );
      parameters.AddInitializer
      (
        "start with all",
        "BP all features in the first round",
        io::Serialization::GetAgent( &m_StartWithAll),
        "False"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
