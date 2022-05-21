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
#include "model/bcl_model_neural_network_selective_backpropagation_leading_sequence.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_data_set_select_columns.h"
#include "model/bcl_model_feature_label_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkSelectiveBackpropagationLeadingSequence::s_Instance
    (
      util::Enumerated< NeuralNetworkSelectiveBackpropagationInterface>::AddInstance
      (
        new NeuralNetworkSelectiveBackpropagationLeadingSequence()
      )
    );

    //! @brief default constructor
    NeuralNetworkSelectiveBackpropagationLeadingSequence::NeuralNetworkSelectiveBackpropagationLeadingSequence() :
      m_ResultStartID( 0),
      m_ResultEndID( 1000),
      m_ResultsSize( 1000),
      m_DatasetSize( 1),
      m_BackPropagationRound( false)
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkSelectiveBackpropagationLeadingSequence copied from this instance
    NeuralNetworkSelectiveBackpropagationLeadingSequence *NeuralNetworkSelectiveBackpropagationLeadingSequence::Clone() const
    {
      return new NeuralNetworkSelectiveBackpropagationLeadingSequence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkSelectiveBackpropagationLeadingSequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkSelectiveBackpropagationLeadingSequence::GetAlias() const
    {
      static const std::string s_name( "LeadingSequence");
      return s_name;
    }

    //! @brief Initialize this object with the rescaled dataset and other important parameters
    //! @param RESCALED_DATA the already-rescaled training dataset
    //! @param OBJECTIVE the actual objective function
    //! @param THREADS # of threads that may be accessing this object at once
    void NeuralNetworkSelectiveBackpropagationLeadingSequence::Initialize
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

      m_DatasetSize = RESCALED_DATA.GetResultsPtr()->GetNumberFeatures();
      m_MappingIds.Reset();
      m_GroupId = storage::Vector< size_t>( m_DatasetSize);

      DataSetSelectColumns sequence_id_selector;
      // get the feature labels so that sequences can be identified
      const FeatureDataSet< char> &ids( *RESCALED_DATA.GetIdsPtr());

      if( !m_DbIdLabel.GetValue().empty() || m_DbIdLabel.GetNumberArguments())
      {
        const FeatureLabelSet model_feature_labels( *ids.GetFeatureLabelSet());
        sequence_id_selector =
          DataSetSelectColumns
          (
            model_feature_labels.GetSize(),
            model_feature_labels.GetPropertyIndices
            (
              util::ObjectDataLabel( m_DbIdLabel.GetValue(), m_DbIdLabel.GetArguments())
            )
          );
      }
      else
      {
        // no labels were specified, use them all
        sequence_id_selector =
          DataSetSelectColumns( ids.GetFeatureSize(), storage::CreateIndexVector( ids.GetFeatureSize()));
      }

      BCL_Assert( sequence_id_selector.GetOutputFeatureSize(), "Could not find " + m_DbIdLabel.ToString() + " in IDs!");

      // locate all the sequence boundaries and ids
      std::string this_sequence_id( sequence_id_selector.GetOutputFeatureSize(), ' ');
      storage::Map< std::string, size_t> map_uniqueid_index;
      size_t n_groups( 0);
      for( size_t i( 0); i < m_DatasetSize; ++i)
      {
        sequence_id_selector( ids( i), &this_sequence_id[ 0]);
        if( !map_uniqueid_index.Has( this_sequence_id))
        {
          m_GroupId( i) = map_uniqueid_index[ this_sequence_id] = n_groups++;
        }
        else
        {
          m_GroupId( i) = map_uniqueid_index[ this_sequence_id];
        }
      }
      m_MappingIds.Resize( n_groups);
      for( size_t i( 0); i < m_DatasetSize; ++i)
      {
        m_MappingIds( m_GroupId( i)).PushBack( i);
      }

      m_LastRoundHighestIndex = m_CurrentRoundHighestIndex =
        storage::Vector< linal::Vector< size_t> >
        (
          n_groups,
          linal::Vector< size_t>( m_ResultsSize, util::GetUndefined< size_t>())
        );

      m_LastRoundHighestPrediction = m_CurrentRoundHighestPrediction =
        storage::Vector< linal::Vector< float> >
        (
          n_groups,
          linal::Vector< float>( m_ResultsSize, math::GetLowestBoundedValue< float>())
        );
    } // Initialize

    //! @brief select whether or not to backpropagate the current feature, and edit the error vector, if desired
    //! @param PREDICTION the prediction that was made by the neural network
    //! @param ERROR Reference to the error vector, (should already have been set to RESULT - PREDICTION)
    //! @param FEATURE_ID id of this feature in the dataset
    //! @param THREAD*_ID id of this thread
    //! @return true if the feature should be backpropagated
    bool NeuralNetworkSelectiveBackpropagationLeadingSequence::ShouldBackpropagate
    (
      const linal::VectorConstInterface< float> &PREDICTION,
      linal::VectorInterface< float> &ERROR,
      const size_t &FEATURE_ID,
      const size_t &THREAD_ID
    )
    {
      size_t unique_id( m_GroupId( FEATURE_ID));
      if( m_BackPropagationRound)
      {
        // get a reference to the result and prediction row, already offset by m_ResultStartID
        linal::VectorConstReference< float> prediction_row( m_ResultsSize, PREDICTION.Begin() + m_ResultStartID);

        bool backprop_should_continue( false);
        for( size_t result_id( 0); result_id < m_ResultsSize; ++result_id)
        {
          const size_t &last_round_max_index( m_CurrentRoundHighestIndex( unique_id)( result_id));
          if( PREDICTION( result_id) > m_CurrentRoundHighestPrediction( unique_id)( result_id))
          {
            m_CurrentRoundHighestPrediction( unique_id)( result_id) = PREDICTION( result_id);
          }
          ERROR( result_id) *=
            PREDICTION( result_id)
            /
            std::max
            (
              float( 1.0e-10),
              m_CurrentRoundHighestPrediction( unique_id)( result_id)
            );
          ERROR( result_id) = std::min( 1.0f, std::max( ERROR( result_id), -1.0f));
        }
        return true;
      }

      // get a reference to the result and prediction row, already offset by m_ResultStartID
      linal::VectorConstReference< float> prediction_row( m_ResultsSize, PREDICTION.Begin() + m_ResultStartID);

      // for each result
      for( size_t result_id( 0); result_id < m_ResultsSize; ++result_id)
      {
        // get the actual result, prediction, and cutoff
        const float prediction( prediction_row( result_id));

        float &this_round_prediction( m_CurrentRoundHighestPrediction( unique_id)( result_id));

        size_t &current_round_max_index( m_CurrentRoundHighestIndex( unique_id)( result_id));

        if( prediction > this_round_prediction)
        {
          current_round_max_index = FEATURE_ID;
          this_round_prediction = prediction;
        }
      }
      // return true if any part of ERROR remains to be backpropagated
      return false;
    }

    //! @brief finalize the current round; occurs only after all threads were already joined
    void NeuralNetworkSelectiveBackpropagationLeadingSequence::FinalizeRound()
    {
      if( m_BackPropagationRound)
      {
        m_CurrentRoundHighestIndex.SetAllElements
        (
          linal::Vector< size_t>( m_ResultsSize, util::GetUndefined< size_t>())
        );

        m_CurrentRoundHighestPrediction.SetAllElements
        (
          linal::Vector< float>( m_ResultsSize, math::GetLowestBoundedValue< float>())
        );
      }
      m_BackPropagationRound = false;
    }

    //! @brief finalize the current round; occurs only after all threads were already joined
    void NeuralNetworkSelectiveBackpropagationLeadingSequence::FinalizeConformation()
    {
      if( m_BackPropagationRound)
      {
        m_CurrentRoundHighestIndex.SetAllElements
        (
          linal::Vector< size_t>( m_ResultsSize, util::GetUndefined< size_t>())
        );

        m_CurrentRoundHighestPrediction.SetAllElements
        (
          linal::Vector< float>( m_ResultsSize, math::GetLowestBoundedValue< float>())
        );
      }
      m_BackPropagationRound = !m_BackPropagationRound;
    }

    //! @brief function to calculate constitution mapping when given the order vector of dataset presentation
    //! @param ORDER vector containing order in which neural network is presented with training dataset
    void NeuralNetworkSelectiveBackpropagationLeadingSequence::ConstitutionMapping
    (
      storage::Vector< size_t> &ORDER
    )
    {
      storage::Vector< size_t> order_vector;
      order_vector.AllocateMemory( ORDER.GetSize());
      for
      (
        storage::Vector< size_t>::const_iterator itr( ORDER.Begin()), itr_end( ORDER.End());
          itr != itr_end;
        ++itr
      )
      {
        size_t cur_index( *itr);
        size_t unique_id( m_GroupId( cur_index));

        if( m_MappingIds( unique_id).FirstElement() == cur_index)
        {
            order_vector.Append( m_MappingIds( unique_id));
        }
      }
      ORDER = order_vector;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkSelectiveBackpropagationLeadingSequence::GetSerializer() const
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
        "label",
        "label to use",
        io::Serialization::GetAgent( &m_DbIdLabel),
        ""
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
