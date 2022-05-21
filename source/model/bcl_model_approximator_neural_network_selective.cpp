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
#include "model/bcl_model_approximator_neural_network_selective.h"

// includes from bcl - sorted alphabetically
#include "iterate/bcl_iterate_reflecting.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_statistics.h"
#include "model/bcl_model_neural_network_selective_backpropagation_default.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ApproximatorNeuralNetworkSelective::s_IterateInstance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorNeuralNetworkSelective( false))
    );
    const util::SiPtr< const util::ObjectInterface> ApproximatorNeuralNetworkSelective::s_PretrainInstance
    (
      util::Enumerated< PretrainNeuralNetworkInterface>::AddInstance( new ApproximatorNeuralNetworkSelective( true))
    );

    //! @brief InputLayerDropoutType as string
    //! @param TYPE the type
    //! @return the string for the type
    const std::string &ApproximatorNeuralNetworkSelective::GetInputLayerDropoutTypeName( const InputLayerDropoutType &TYPE)
    {
      static const std::string s_names[] =
      {
        "Zero",
        "Noise",
        "CopySingleRandomFeature",
        "CopyEachRandomFeature",
        "CopySingleRandomPeerFeature",
        "CopyEachRandomPeerFeature",
        "CopySingleBalancedFeature",
        "CopyEachBalancedFeature",
        GetStaticClassName< ApproximatorNeuralNetworkSelective::InputLayerDropoutType>()
      };
      return s_names[ size_t( TYPE)];
    }

    //! @brief default constructor
    ApproximatorNeuralNetworkSelective::ApproximatorNeuralNetworkSelective( const bool &PRETRAIN) :
      m_UpdateEveryNthFeature( 0),
      m_IterationsPerRMSDMessage( 1),
      m_Shuffle( false),
      m_ConnectionDensity( 1.0),
      m_AlignCutoff( false),
      m_RescaleOutputDynamicRange( true),
      m_DataSetRangePosition( 0),
      m_TransferFunction(),
      m_DataSelector( NeuralNetworkSelectiveBackpropagationDefault()),
      m_IsPretrainer( PRETRAIN),
      m_NumIterations( 0),
      m_Balance( false),
      m_BalanceMaxRepeatedFeatures( 1000000),
      m_BalanceMaxOversampling( 1.0),
      m_NoiseZScore( 0.0),
      m_InputDropoutType( e_Zero),
      m_RescaleType( RescaleFeatureDataSet::e_AveStd)
    {
    }

    //! @brief constructor from training data, transfer, rescale, and objective functions
    //! @param TRAINING_DATA data to train the NeuralNetwork on
    //! @param UPDATE_EVERY_NTH_FEATURE how often the weights get updated
    //! @param ARCHITECTURE the # of neurons in each hidden layer of the network
    //! @param TRANSFER_FUNCTION ShPtr to the transfer function between input and output of each neuron
    //! @param OBJECTIVE_FUNCTION ShPtr to objective function
    //! @param WEIGHT_UPDATE_FUNCTION method by which to update the weights
    //! @param ITERATIONS_PER_RMSD_REPORT # iterations per report of the rmsd
    ApproximatorNeuralNetworkSelective::ApproximatorNeuralNetworkSelective
    (
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t UPDATE_EVERY_NTH_FEATURE,
      const storage::Vector< size_t> &ARCHITECTURE,
      const util::Implementation< TransferFunctionInterface> &TRANSFER_FUNCTION,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
      const util::Implementation< NeuralNetworkUpdateWeightsInterface> &WEIGHT_UPDATE_FUNCTION,
      const size_t &ITERATIONS_PER_RMSD_REPORT,
      const RescaleFeatureDataSet::TypeEnum RESCALE_TYPE
    ) :
      ApproximatorBase( OBJECTIVE_FUNCTION),
      m_UpdateEveryNthFeature( UPDATE_EVERY_NTH_FEATURE),
      m_IterationsPerRMSDMessage( ITERATIONS_PER_RMSD_REPORT),
      m_Shuffle( false),
      m_ConnectionDensity( 1.0),
      m_AlignCutoff( false),
      m_RescaleOutputDynamicRange( true),
      m_DataSetRangePosition( 0),
      m_HiddenArchitecture( ARCHITECTURE),
      m_TransferFunction( TRANSFER_FUNCTION),
      m_WeightUpdateType( WEIGHT_UPDATE_FUNCTION),
      m_DataSelector( NeuralNetworkSelectiveBackpropagationDefault()),
      m_IsPretrainer( false),
      m_NumIterations( 0),
      m_Balance( false),
      m_BalanceMaxRepeatedFeatures( 1000000),
      m_BalanceMaxOversampling( 1.0),
      m_NoiseZScore( 0.0),
      m_InputDropoutType( e_Zero),
      m_RescaleType( RESCALE_TYPE)
    {
      // set and rescale training data set
      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorNeuralNetworkSelective copied from this instance
    ApproximatorNeuralNetworkSelective *ApproximatorNeuralNetworkSelective::Clone() const
    {
      return new ApproximatorNeuralNetworkSelective( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorNeuralNetworkSelective::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorNeuralNetworkSelective::GetAlias() const
    {
      static const std::string s_Name( "NeuralNetworkSelective");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorNeuralNetworkSelective::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      m_TrainingData = DATA;
      if( !DATA->GetFeaturesPtr()->IsRescaled())
      {
        DATA->GetFeatures().Rescale( NeuralNetwork::s_DefaultInputRange, m_RescaleType);
      }
      if( !DATA->GetResultsPtr()->IsRescaled())
      {
        DATA->GetResults().Rescale
        (
          m_RescaleOutputDynamicRange
          ? m_TransferFunction->GetDynamicOutputRange()
          : m_TransferFunction->GetOutputRange()
        );
      }
      m_RescaleOutputLastRound = DATA->GetResultsPtr()->GetScaling();

      m_Order.Resize( m_TrainingData->GetSize());
      for( size_t i( 0), n_features( m_TrainingData->GetSize()); i < n_features; ++i)
      {
        m_Order( i) = i;
      }
      if( m_Balance)
      {
        BalanceFeatures( m_BalanceMaxRepeatedFeatures, m_BalanceMaxOversampling);
      }
      if( int( m_InputDropoutType) >= int( s_InputLayerDropoutFirstClassBasedMethod))
      {
        if( m_Dropout.IsEmpty() || !m_Dropout( 0))
        {
          BCL_MessageCrt
          (
            m_InputDropoutType.GetString() + " has no effect because no neurons are to be dropped in the input layer"
          );
        }
        if( !util::IsDefined( m_ObjectiveFunction->GetThreshold()))
        {
          const std::string &old_input_dropout_type( m_InputDropoutType.GetString());
          m_InputDropoutType
            = (
                m_InputDropoutType == e_CopyEachRandomPeerFeature || m_InputDropoutType == e_CopyEachBalancedFeature
                ? e_CopyEachRandomFeature
                : e_CopySingleRandomFeature
               );
          BCL_MessageCrt
          (
            old_input_dropout_type
            + " has no effect because the objective function does not have a cutoff, switching to "
            + m_InputDropoutType.GetString()
          );
        }
        DetermineResultClasses();
      }

      BCL_MessageStd
      (
        "Setting up training data with " + util::Format()( m_TrainingData->GetSize()) + " points"
      );

      // initialize the dropout, if necessary
      InitializeDropout();

      // if the current instance is NOT used through IterateInterfaceFromFile
      // then initialize the Iterate properly
      if( !IsTrainingContinued())
      {
        if( m_Pretrainer.IsDefined())
        {
          m_Pretrainer->SetIds( ApproximatorBase::GetIdCode());
          m_Pretrainer->SetFeatures( ApproximatorBase::GetFeatureCode());
          m_Pretrainer->SetResults( ApproximatorBase::GetResultCode());

          BCL_MessageStd( "Pretraining with method: " + m_Pretrainer->GetString());
          util::ShPtr< NeuralNetwork> initial_network( m_Pretrainer->PretrainNetwork( DATA, m_ObjectiveFunction));
          // ensure that a model was pretrained
          if( !initial_network.IsDefined())
          {
            BCL_Exit( "Pretraining failed!", -1);
          }

          // copy the bias and weight from the model
          m_Bias = initial_network->GetBias();
          m_Weight = initial_network->GetWeight();

          // copy the hidden architecture from the model
          linal::Vector< size_t> old_hidden_architecture( m_HiddenArchitecture);
          m_HiddenArchitecture = initial_network->GetArchitecture();
          // remove the first and last elements, which are the input and output layer sizes, respectively
          m_HiddenArchitecture.RemoveElements( 0, 1);
          m_HiddenArchitecture.RemoveElements( m_HiddenArchitecture.GetSize() - 1, 1);

          // need to adjust weights if performing dropout in this layer
          if( m_NumberToDrop.Sum())
          {
            // must multiply m_Weight( i) by (1-m_Dropout( i))
            for( size_t layer( 0), n_layers( m_Weight.GetSize()); layer < n_layers; ++layer)
            {
              if( !m_NumberToDrop( layer))
              {
                continue;
              }
              const size_t neurons_this_layer( m_Weight( layer).GetNumberCols());
              const float kept_fraction( float( neurons_this_layer) / float( neurons_this_layer - m_NumberToDrop( layer)));
              m_Weight( layer) *= kept_fraction;
            }
          }

          // make sure that either no hidden architecture was specified or that it agreed with the architecture implicit
          // in the model
          if
          (
            !old_hidden_architecture.IsEmpty()
            && linal::Vector< size_t>( m_HiddenArchitecture) != old_hidden_architecture
          )
          {
            BCL_MessageCrt
            (
              "Ignoring hidden architecture given in iterate; using implicit architecture from initial network file"
            );
          }
        }
      }

      // set the data ranges up for threading (if applicable)
      SetupDataSetRanges();

      m_DataSelector->Initialize( *m_TrainingData, *m_ObjectiveFunction->GetImplementation(), m_NumberThreads);
      m_DataSelector->ConstitutionMapping( m_Order);

      // if the current instance is NOT used through IterateInterfaceFromFile
      // then initialize the Iterate properly
      if( !IsTrainingContinued())
      {
        // if the initial weights and bias were not read in from a file, they need to be initialized
        if( !m_Pretrainer.IsDefined())
        {
          // set the architecture (m_Bias, m_Weights, etc) up using m_HiddenArchitecture
          SetupArchitecture();
        }
        else
        {
          // destroy the pretrainer to free up memory
          m_Pretrainer = m_Pretrainer.GetLabel();
        }
        // weight updaters need to be initialized regardless
        InitializeWeightUpdaters();
      }

      // update neurons for dropout
      UpdateDroppedNeurons();

    } // SetTrainingData

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorNeuralNetworkSelective::GetCurrentModel() const
    {
      if( m_NumberToDrop.Sum())
      {
        // must multiply m_Weight( i) by (1-m_Dropout( i)).  To avoid unnecessary creation/destruction of matrices,
        // use m_SlopesWeight
        for( size_t layer( 0), n_layers( m_Weight.GetSize()); layer < n_layers; ++layer)
        {
          m_SlopesWeight( 0)( layer) = m_Weight( layer);
          if( m_NumberToDrop( layer))
          {
            const float n_dropped( m_ChosenDropped( layer).GetSize());
            if( layer || m_InputDropoutType == e_Noise || m_InputDropoutType == e_Zero || m_InputDropoutType == e_CopyEachRandomFeature)
            {
              const float kept_fraction( 1.0 - n_dropped / float( m_Weight( layer).GetNumberCols()));
              m_SlopesWeight( 0)( layer) *= kept_fraction;
            }
//            else if( m_InputDropoutType == e_CopyEachBalancedFeature || m_InputDropoutType == e_CopySingleBalancedFeature)
//            {
//              const float n_classes( m_PeerFeatures.GetSize());
//              const float kept_fraction( 1.0 - n_dropped / float( m_Weight( layer).GetNumberCols() * n_classes));
//              m_SlopesWeight( 0)( layer) *= kept_fraction;
//            }
          }
        }
        return util::ShPtr< Interface>
               (
                 new NeuralNetwork
                 (
                   GetRescaleFeatureDataSet(),
                   m_RescaleOutputLastRound,
                   m_Bias,
                   m_SlopesWeight( 0),
                   m_TransferFunction
                 )
               );
      }
      // make a new model out of the current data members
      return util::ShPtr< Interface>
      (
        new NeuralNetwork
        (
          GetRescaleFeatureDataSet(),
          m_RescaleOutputLastRound,
          m_Bias,
          m_Weight,
          m_TransferFunction
        )
      );
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorNeuralNetworkSelective::GetCurrentApproximation() const
    {
      util::ShPtr< Interface> model( GetCurrentModel());
      return
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>
          (
            model,
            m_ObjectiveFunction->operator ()( model)
          )
        );
    }

    //! @brief the main operation, pretrains a neural network
    //! @param DATA the data for use in pretraining
    //! @param OBJECTIVE ShPtr to the objective function for the network
    util::ShPtr< NeuralNetwork> ApproximatorNeuralNetworkSelective::PretrainNetwork
    (
      util::ShPtr< descriptor::Dataset> &DATA,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
    )
    {
      this->SetObjectiveFunction( OBJECTIVE);
      this->SetTrainingData( DATA);
      while( this->GetTracker().GetIteration() < m_NumIterations)
      {
        this->Next();
      }
      return GetCurrentModel();
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorNeuralNetworkSelective::Next()
    {
      // handle shuffling
      if( m_Shuffle)
      {
        m_Order.Shuffle();
      }

      util::ShPtrVector< sched::JobInterface> all_jobs( m_NumberThreads);

      const size_t group_id( 0);

      // store the id of each thread for use by the jobs and reset all results tallying objects.
      storage::Vector< size_t> thread_ids( m_NumberThreads);
      for( size_t thread_number( 0); thread_number < m_NumberThreads; ++thread_number)
      {
        thread_ids( thread_number) = thread_number;
        m_RMSDError( thread_number) = float( 0.0);
      }

      // run through the data set once, updating weights after every epoch (m_UpdateEveryNthFeature features, data set size by default)
      // launch jobs to train the NN with the features for this epoch
      for( size_t thread_number( 0); thread_number < m_NumberThreads; ++thread_number)
      {
        all_jobs( thread_number)
          = util::ShPtr< sched::JobInterface>
            (
              new sched::UnaryFunctionJobWithData
              <
                const size_t,
                void,
                ApproximatorNeuralNetworkSelective
              >
              (
                group_id,
                *this,
                &ApproximatorNeuralNetworkSelective::TrainThread,
                thread_ids( thread_number),
                sched::JobInterface::e_READY,
                NULL
              )
            );
      }

      m_DataSetRangePosition = 0;
      size_t n_updated = 0;
      while( m_DataSetRangePosition < m_DataSetRanges.GetSize())
      {
        if( m_NumberThreads > size_t( 1))
        {
          for( size_t thread_number( 0); thread_number < m_NumberThreads; thread_number++)
          {
            sched::GetScheduler().SubmitJob( all_jobs( thread_number));
          }

          sched::GetScheduler().Join( all_jobs( 0));
        }
        else
        {
          // just run the job directly to avoid the scheduler's overhead
          all_jobs( 0)->Run();
        }

        size_t n_updated_this_time( 0);
        n_updated_this_time += m_NumberFeaturesUpdated( 0);
        m_NumberFeaturesUpdated( 0) = 0;
        for( size_t thread_number( 1); thread_number < m_NumberThreads; thread_number++)
        {
          sched::GetScheduler().Join( all_jobs( thread_number));
          n_updated_this_time += m_NumberFeaturesUpdated( thread_number);

          if( m_NumberFeaturesUpdated( thread_number))
          {
            // accumulate the slopes from the newly joined jobs with the slopes from the earlier jobs
            for
            (
              size_t hidden_layer_number( 0), number_hidden_layers( m_HiddenArchitecture.GetSize());
              hidden_layer_number < number_hidden_layers;
              ++hidden_layer_number
            )
            {
              m_SlopesBias( 0)( hidden_layer_number) += m_SlopesBias( thread_number)( hidden_layer_number);
              m_SlopesWeight( 0)( hidden_layer_number) += m_SlopesWeight( thread_number)( hidden_layer_number);
            }
            m_NumberFeaturesUpdated( thread_number) = 0;
          }
        }

        if( n_updated_this_time)
        {
          UpdateWeights();
          n_updated += n_updated_this_time;
        }
        else
        {
          UpdateDroppedNeurons();
        }

        // move on the the next epoch
        m_DataSetRangePosition += m_NumberThreads;
      }

      if( ( this->GetTracker().GetIteration() % m_IterationsPerRMSDMessage) == size_t( 0))
      {
        for( size_t thread_number( 1); thread_number < m_NumberThreads; thread_number++)
        {
          m_RMSDError( 0) += m_RMSDError( thread_number);
        }

        m_RMSDError( 0) /= m_TrainingData->GetResultSize();
        m_RMSDError( 0) /= m_TrainingData->GetSize();
        m_RMSDError( 0) = std::max( m_RMSDError( 0), float( 0.0));
        BCL_MessageStd
        (
          "Training relative RMSD for iteration " + util::Format()( this->GetTracker().GetIteration())
          + ( n_updated != m_Order.GetSize() ? ", " + util::Format()( n_updated) + " features trained" : std::string())
          + ": " + util::Format()( m_RMSDError( 0))
        );
      }

      m_NumberFeaturesUpdated = 0;

      m_DataSelector->FinalizeRound();

      m_DataSetRangePosition = 0;

      // get current model
      util::ShPtr< NeuralNetwork> current_model( GetCurrentModel());
      current_model->SetRescaleOutput( GetRescaleResultDataSet());

      // determine inaccuracy counts
      if( m_IterationWeightUpdateType.IsDefined())
      {
        linal::Matrix< float> old_weight( m_Weight( 0));
        ( *m_IterationWeightUpdaters( 0))( m_Weight( 0));
        float total_change( 0), total_weight( 0);
        for
        (
          const float *itr_old( old_weight.Begin()), *itr_old_end( old_weight.End()), *itr_new( m_Weight( 0).Begin());
          itr_old != itr_old_end;
          ++itr_old, ++itr_new
        )
        {
          total_change += math::Absolute( *itr_old - *itr_new);
          total_weight += math::Absolute( *itr_new);
        }
        BCL_MessageStd
        (
          "Average weight change: " + util::Format()( total_change / float( m_Weight( 0).GetNumberOfElements()))
          + " Average weight now : " + util::Format()( total_weight / float( m_Weight( 0).GetNumberOfElements()))
        );
      }

      if( m_IsPretrainer)
      {
        this->GetTracker().Track
        (
          util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
          (
            new storage::Pair< util::ShPtr< Interface>, float>
            (
              util::ShPtr< Interface>(),
              float( 0.0)
            )
          )
        );
        // act as a pretrainer for the specified number of steps, no need to create a model
        return;
      }

      // compute results on monitoring dataset
      FeatureDataSet< float> predicted_results( m_ObjectiveFunction->Predict( util::ShPtr< Interface>( current_model)));

      const float objective_result( m_ObjectiveFunction->Evaluate( predicted_results));
      if( m_AlignCutoff)
      {
        // optimize the rescaling function
        m_RescaleOutputLastRound = m_ObjectiveFunction->OptimizeRescalingFunction( GetRescaleResultDataSet(), predicted_results);
        current_model->SetRescaleOutput( m_RescaleOutputLastRound);
      }

      // combine it with the objective function evaluation
      this->GetTracker().Track
      (
        util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
        (
          new storage::Pair< util::ShPtr< Interface>, float>( current_model, objective_result)
        )
      );
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorNeuralNetworkSelective::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "trains a neural network (see http://en.wikipedia.org/wiki/Artificial_neural_network)"
      );

      parameters.AddInitializer
      (
        "transfer function",
        "function that translates input from neurons in the prior layer into the output of each hidden layer",
        io::Serialization::GetAgent( &m_TransferFunction),
        "Sigmoid"
      );
      parameters.AddInitializer
      (
        "weight update",
        "algorithm used to update the weights",
        io::Serialization::GetAgent( &m_WeightUpdateType),
        "Resilient"
      );
      parameters.AddOptionalInitializer
      (
        "bias update",
        "algorithm used to update the biases; if omitted, the weight update type is used",
        io::Serialization::GetAgent( &m_BiasUpdateType)
      );
      parameters.AddOptionalInitializer
      (
        "iteration weight update",
        "algorithm used to update the weights after every run through the data",
        io::Serialization::GetAgent( &m_IterationWeightUpdateType)
      );
      parameters.AddInitializer
      (
        "steps per update",
        "# of features seen between each update of weights (set to 0 to use the size of the training data set)",
        io::Serialization::GetAgent( &m_UpdateEveryNthFeature),
        "0"
      );
      parameters.AddOptionalInitializer
      (
        "hidden architecture",
        "# of neurons in each hidden layer, e.g. hidden architecture(100) declares a single hidden layer w/ 100 neurons. "
        "This variable can be omitted to train a single layer perception (logistic regression)",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_HiddenArchitecture,
          io::Serialization::GetAgentWithMin( size_t( 1)),
          0
        )
      );
      parameters.AddInitializer
      (
        "rmsd report frequency",
        "# of iterations between reports of rmsd on the last training batch",
        io::Serialization::GetAgentWithMin( &m_IterationsPerRMSDMessage, size_t( 1)),
        "1"
      );
      parameters.AddOptionalInitializer
      (
        "initial network",
        "method by which to create the initial network, if not randomly",
        io::Serialization::GetAgent( &m_Pretrainer)
      );
      parameters.AddInitializer
      (
        "shuffle",
        "primarily for non-batch update; if true, shuffle the order or data points between each run through the data",
        io::Serialization::GetAgent( &m_Shuffle),
        "False"
      );
      parameters.AddInitializer
      (
        "data selector",
        "method to use to select which data to backpropagate",
        io::Serialization::GetAgent( &m_DataSelector),
        "All"
      );
      parameters.AddInitializer
      (
        "connection density",
        "Fraction of connections made between input and hidden layers",
        io::Serialization::GetAgentWithRange( &m_ConnectionDensity, 0.0, 1.0),
        "1"
      );
      parameters.AddInitializer
      (
        "balance",
        "Whether to automatically balance each class (as defined by the objective function's cutoff, if applicable)",
        io::Serialization::GetAgent( &m_Balance),
        "False"
      );
      parameters.AddInitializer
      (
        "balance max repeats",
        "Applies only if balance=True; absolute maximum number of times that a feature can be repeated in order to reach"
        "the targeted ratio of positives to negatives",
        io::Serialization::GetAgent( &m_BalanceMaxRepeatedFeatures),
        "1000000"
      );
      parameters.AddInitializer
      (
        "balance target ratio",
        "Applies only if balance=True; target ratio between most-common and underrepresented class in the dataset "
        "achieved by data replication; to simulate normal balancing; this should be 1, but smaller values may yield "
        "more general models",
        io::Serialization::GetAgentWithRange( &m_BalanceMaxOversampling, float( 0.0), float( 1.0)),
        "1"
      );
      parameters.AddOptionalInitializer
      (
        "dropout",
        "fraction of neurons to have \"dropout\" (set to 0) for an entire weight-update session, per layer",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_Dropout,
          io::Serialization::GetAgentWithRange( float( 0.0), float( 1.0)),
          0
        )
      );
      parameters.AddOptionalInitializer
      (
        "dropout partitions",
        "# of partitions between neurons in each layer.  The number of neurons dropped within each partition will always"
        "be the same, so this allows difference regions of hidden neurons to learn different functionalities, which is"
        "particularly useful for multi-output ANNs",
        io::Serialization::GetAgentWithSizeLimits( &m_DropoutPartitions, size_t( 0))
      );
      parameters.AddInitializer
      (
        "scaling",
        "Type of input scaling. Normally AveStd works best, but MinMax and None may also be used in some circumstances",
        io::Serialization::GetAgent( &m_RescaleType),
        "AveStd"
      );
      parameters.AddInitializer
      (
        "input noise",
        "Amount of noise to apply to input features, in units of z-score",
        io::Serialization::GetAgentWithRange( &m_NoiseZScore, 0.0, 3.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "rescale output dynamic range",
        "if true, output will be rescaled to the dynamic range of the transfer function. Useful for "
        "regression, where the output of the function being predicted could feasibly extend beyond the range seen in the "
        "the training set by a marginal amount. It may also result in smaller weights. "
        "Classification problems often benefit from using the False setting",
        io::Serialization::GetAgent( &m_RescaleOutputDynamicRange),
        "True"
      );
      parameters.AddInitializer
      (
        "input dropout type",
        "Method of applying dropout to the input layer. Descriptions for each method:\n"
        "Zero                         <-- Set the dropped neurons to zero; fast, but biased\n"
        "Noise,                       <-- Set the neurons to a random gaussian value sampled according to the input mean/std\n"
        "CopySingleRandomFeature      <-- Copy dropped values from a single randomly-selected feature\n"
        "CopyEachRandomFeature        <-- Copy each dropped value from a randomly-selected feature (training example)\n"
        "CopySingleRandomPeerFeature  <-- Copy values from a single, randomly-selected peer feature\n"
        "CopyEachRandomPeerFeature    <-- Copy each value from a randomly-selected peer feature\n"
        "A peer is a training example that has the same class for every result column.  Use of input dropout types that "
        "have Peer in the name requires that a classification-based objective function is used (typically any objective "
        "function with a cutoff). The output class is a binary value of whether the column is above or below the cutoff"
        ". The balanced variants (CopySingleBalancedFeature, CopyEachBalancedFeature) select a random feature of a "
        " random class.",
        io::Serialization::GetAgent( &m_InputDropoutType),
        "Zero"
      );

      if( !m_IsPretrainer)
      {
        parameters.AddInitializer
        (
          "align cutoff",
          "Adjust output scaling automatically according to the objective functions.  This will only have an effect for "
          "objective functions such as enrichment and fpp vs ppv that have an internal FPR or similar cutoff.",
          io::Serialization::GetAgent( &m_AlignCutoff),
          "False"
        );
        parameters.AddInitializer
        (
          "objective function",
          "function that evaluates the model after each batch step",
          io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
          "RMSD"
        );
      }
      else
      {
        // pretraining-specific parameters
        parameters.AddInitializer
        (
          "iterations",
          "Number of complete iterations through the data to pretrain the network",
          io::Serialization::GetAgent( &m_NumIterations)
        );
      }

      return parameters;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorNeuralNetworkSelective::Read( std::istream &ISTREAM)
    {
      descriptor::Dataset::s_Instance.IsDefined();

      // read members
      io::Serialize::Read( m_UpdateEveryNthFeature, ISTREAM);
      io::Serialize::Read( m_HiddenArchitecture, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Weight, ISTREAM);
      io::Serialize::Read( m_TransferFunction, ISTREAM);
      io::Serialize::Read( m_WeightUpdateType, ISTREAM);
      io::Serialize::Read( m_IterationWeightUpdateType, ISTREAM);
      io::Serialize::Read( m_BiasUpdaters, ISTREAM);
      io::Serialize::Read( m_WeightUpdaters, ISTREAM);
      io::Serialize::Read( m_IterationWeightUpdaters, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorNeuralNetworkSelective::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_UpdateEveryNthFeature, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HiddenArchitecture, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Weight, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TransferFunction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WeightUpdateType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IterationWeightUpdateType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BiasUpdaters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WeightUpdaters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IterationWeightUpdaters, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! run forward through ANN and compute hidden terms m_Hidden (test)
    void ApproximatorNeuralNetworkSelective::CalcHiddenTerms
    (
      const FeatureReference< float> &FEATURE,
      const size_t &THREAD_ID,
      const size_t &FEATURE_ID
    )
    {
      // get the thread-local errors, hidden, and hidden input arrays
      storage::Vector< linal::Vector< float> > &hidden( m_Hidden( THREAD_ID));
      storage::Vector< linal::Vector< float> > &hidden_input( m_HiddenInput( THREAD_ID));

      // input layer
      // perform hidden_input( 0) = m_Weight( 0) * FEATURE + m_Bias( 0); without creating new vectors
      hidden_input( 0) = m_Bias( 0);
      if( m_NumberToDrop( 0) || m_NoiseZScore)
      {
        linal::Vector< float> &feature_with_dropout( m_FeatureWithDropout( THREAD_ID));
        std::copy( FEATURE.Begin(), FEATURE.End(), feature_with_dropout.Begin());
        linal::MatrixConstReference< float> allfeatures( m_TrainingData->GetFeaturesReference());
        if( m_NoiseZScore)
        {
          random::UniformDistribution &rng( m_ThreadRandomNumberGenerators( THREAD_ID));
          for( float *itr( feature_with_dropout.Begin()), *itr_end( feature_with_dropout.End()); itr != itr_end; ++itr)
          {
            *itr += rng.RandomGaussian( 0.0, m_NoiseZScore);
          }
        }
        if( m_NumberToDrop( 0))
        {
          const storage::Vector< size_t> &drop_ids( m_ChosenDropped( 0));
          if( m_InputDropoutType == e_Noise)
          {
            random::UniformDistribution &rng( m_ThreadRandomNumberGenerators( THREAD_ID));
            for
            (
              storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
              itr_drop != itr_drop_end;
              ++itr_drop
            )
            {
              feature_with_dropout( *itr_drop) = rng.RandomGaussian( 0.0, 0.5);
            }
          }
          else if( m_InputDropoutType == e_Zero)
          {
            for
            (
              storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
              itr_drop != itr_drop_end;
              ++itr_drop
            )
            {
              feature_with_dropout( *itr_drop) = 0.0;
            }
          }
          else if( m_InputDropoutType == e_CopySingleRandomFeature || m_InputDropoutType == e_CopyEachRandomFeature)
          {
            random::UniformDistribution &rng( m_ThreadRandomNumberGenerators( THREAD_ID));
            const size_t max_feature_to_pick( m_TrainingData->GetSize() - 1);
            if( m_InputDropoutType == e_CopySingleRandomFeature)
            {
              const size_t random_feature( rng.Random( max_feature_to_pick));
              linal::VectorConstReference< float> selected_feature( allfeatures.GetRow( random_feature));
              for
              (
                storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
                itr_drop != itr_drop_end;
                ++itr_drop
              )
              {
                feature_with_dropout( *itr_drop) = selected_feature( *itr_drop);
              }
            }
            else
            {
              for
              (
                storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
                itr_drop != itr_drop_end;
                ++itr_drop
              )
              {
                feature_with_dropout( *itr_drop) = allfeatures( rng.Random( max_feature_to_pick), *itr_drop);
              }
            }
          }
          else if
          (
            m_InputDropoutType == e_CopySingleRandomPeerFeature
            || m_InputDropoutType == e_CopyEachRandomPeerFeature
          )
          {
            // peer-based dropout types
            random::UniformDistribution &rng( m_ThreadRandomNumberGenerators( THREAD_ID));
            const storage::Vector< size_t> &peers( m_PeerFeatures( m_ResultClass( FEATURE_ID)));
            const size_t max_feature_to_pick( peers.GetSize() - 1);
            if( m_InputDropoutType == e_CopySingleRandomPeerFeature)
            {
              const size_t random_feature( peers( rng.Random( max_feature_to_pick)));
              linal::VectorConstReference< float> selected_feature( allfeatures.GetRow( random_feature));
              for
              (
                storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
                itr_drop != itr_drop_end;
                ++itr_drop
              )
              {
                feature_with_dropout( *itr_drop) += selected_feature( *itr_drop);
              }
            }
            else // if( m_InputDropoutType == e_CopyEachRandomPeerFeature)
            {
              for
              (
                storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
                itr_drop != itr_drop_end;
                ++itr_drop
              )
              {
                feature_with_dropout( *itr_drop) = allfeatures( peers( rng.Random( max_feature_to_pick)), *itr_drop);
              }
            }
          }
          else if( m_InputDropoutType == e_CopySingleBalancedFeature)
          {
            // single feature, selected randomly from randomly selected class dropout types
            random::UniformDistribution &rng( m_ThreadRandomNumberGenerators( THREAD_ID));
            // choose which class to select the feature from
            const storage::Vector< size_t> &chosen_class( m_PeerFeatures( rng.Random( m_PeerFeatures.GetSize() - 1)));
            const size_t random_feature_id( chosen_class( rng.Random( chosen_class.GetSize() - 1)));
            linal::VectorConstReference< float> selected_feature( allfeatures.GetRow( random_feature_id));
            for
            (
              storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
              itr_drop != itr_drop_end;
              ++itr_drop
            )
            {
              feature_with_dropout( *itr_drop) = selected_feature( *itr_drop);
            }
          }
          else if( m_InputDropoutType == e_CopyEachBalancedFeature)
          {
            // single feature, selected randomly from randomly selected class dropout types
            random::UniformDistribution &rng( m_ThreadRandomNumberGenerators( THREAD_ID));
            // choose which class to select the feature from
            const size_t max_class( m_PeerFeatures.GetSize() - 1);
            for
            (
              storage::Vector< size_t>::const_iterator itr_drop( drop_ids.Begin()), itr_drop_end( drop_ids.End());
              itr_drop != itr_drop_end;
              ++itr_drop
            )
            {
              const storage::Vector< size_t> &chosen_class( m_PeerFeatures( rng.Random( max_class)));
              const size_t random_feature_id( chosen_class( rng.Random( chosen_class.GetSize() - 1)));
              feature_with_dropout( *itr_drop) = allfeatures( random_feature_id, *itr_drop);
            }
          }
        }
        linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( 0), m_Weight( 0), feature_with_dropout);
      }
      else
      {
        linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( 0), m_Weight( 0), FEATURE);
      }

      // run hidden_input through the transfer functions to get the 1st hidden layer output
      m_TransferFunction->F( hidden( 0), hidden_input( 0));

      // remaining layers
      for( size_t k( 1), number_hidden_layers( hidden.GetSize()); k < number_hidden_layers; ++k)
      {
        // perform hidden_input( k) = m_Weight( k) * hidden(k-1) + m_Bias( k); without creating new vectors
        hidden_input( k) = m_Bias( k);
        if( m_NumberToDrop( k))
        {
          linal::Vector< float> &hidden_dropped( hidden( k - 1));
          const storage::Vector< size_t> &drop_ids( m_ChosenDropped( k));
          for
          (
            storage::Vector< size_t>::const_iterator itr_dropped( drop_ids.Begin()), itr_dropped_end( drop_ids.End());
            itr_dropped != itr_dropped_end;
            ++itr_dropped
          )
          {
            hidden_dropped( *itr_dropped) = 0.0;
          }
        }
        linal::GetDefaultOperations< float>().VectorPlusEqualsMatrixTimesVector( hidden_input( k), m_Weight( k), hidden( k - 1));

        // run hidden_input through the transfer functions to get the k+1th hidden layer output
        m_TransferFunction->F( hidden( k), hidden_input( k));

      }
    } // CalcHiddenTerms

    //! run backward through ANN and compute err terms m_Errors (train)
    //! @return true if the error should be backpropagated
    bool ApproximatorNeuralNetworkSelective::CalcErrorTerms
    (
      const FeatureReference< float> &RESULT,
      const size_t &THREAD_ID,
      const size_t &FEATURE_ID
    )
    {
      // get the thread-local errors, hidden, and hidden input arrays
      storage::Vector< linal::Vector< float> > &errors( m_Errors( THREAD_ID));
      const storage::Vector< linal::Vector< float> > &hidden( m_Hidden( THREAD_ID));
      const storage::Vector< linal::Vector< float> > &hidden_input( m_HiddenInput( THREAD_ID));

      // compute difference
      const linal::Vector< float> &predicted( hidden.LastElement());
      linal::Vector< float> &storage( errors.LastElement());
      storage = RESULT;
      storage -= predicted;
      const float error( storage.Norm());
      m_RMSDError( THREAD_ID) += error;

      // decide whether or not to backprop this feature
      // may also decide not to backprop some of the oclumns
      if( !m_DataSelector->ShouldBackpropagate( predicted, storage, FEATURE_ID, THREAD_ID))
      {
        return false;
      }

      // output layer

      // the following code performs this operation
      // errors.LastElement() = ( RESULT - hidden.LastElement())
      //                          * m_TransferFunction->dF( hidden_input.LastElement(), hidden.LastElement())
      // without creating any (expensive) temporary vectors
      m_TransferFunction->MultiplyBydF
      (
        errors.LastElement(),
        hidden_input.LastElement(),
        hidden.LastElement()
      );

      // all other layers
      for( size_t i( hidden.GetSize() - 1); i > 0; --i)
      {
        // errors( i - 1) = errors( i) * m_Weight( i)
        linal::GetDefaultOperations< float>().VectorEqualsVectorTimesMatrix( errors( i - 1), errors( i), m_Weight( i));
        m_TransferFunction->MultiplyBydF
        (
          errors( i - 1),
          hidden_input( i - 1),
          hidden( i - 1)
        );
      }
      return true;
    } // CalcErrorTerms

    //! compute changes m_SlopeBias/Weight to be applied on ANN (train)
    void ApproximatorNeuralNetworkSelective::CalcChangeTerms( const FeatureReference< float> &FEATURE, const size_t &THREAD_ID)
    {
      // get the thread-local arrays
      const storage::Vector< linal::Vector< float> > &errors( m_Errors( THREAD_ID));
      const storage::Vector< linal::Vector< float> > &hidden( m_Hidden( THREAD_ID));
      storage::Vector< linal::Vector< float> > &slopes_bias( m_SlopesBias( THREAD_ID));
      storage::Vector< linal::Matrix< float> > &slopes_weight( m_SlopesWeight( THREAD_ID));

      // first layer
      // m_SlopesWeight( THREAD_ID)( 0) += math::OuterProduct( errors( 0), FEATURE)
      linal::AddOuterProductToMatrix( slopes_weight( 0), errors( 0), FEATURE);

      // all other layers
      for( size_t i( 1), number_hidden_layers( hidden.GetSize()); i < number_hidden_layers; ++i)
      {
        // m_SlopesWeight( THREAD_ID)( i) += math::OuterProduct( errors( i), FEATURE)
        linal::AddOuterProductToMatrix( slopes_weight( i), errors( i), hidden( i - 1));
      }

      // add the errors to the slopes for the biases directly
      for( size_t i( 0), number_hidden_layers( hidden.GetSize()); i < number_hidden_layers; ++i)
      {
        slopes_bias( i) += errors( i);
      }
    } // CalcChangeTerms

    //! train ANN with a feature
    void ApproximatorNeuralNetworkSelective::TrainThread( const size_t &THREAD_ID)
    {
      // zero out the slopes vectors / matrices
      for
      (
        size_t layer_number( 0), number_hidden_layers( m_SlopesBias( THREAD_ID).GetSize());
        layer_number < number_hidden_layers;
        ++layer_number
      )
      {
        m_SlopesBias( THREAD_ID)( layer_number) = float( 0);
        m_SlopesWeight( THREAD_ID)( layer_number) = float( 0);
      }

      // only run the data for this range if it actually existed
      // it may not if m_DataSetRanges.GetSize() is not a multiple of the number of threads
      if( THREAD_ID + m_DataSetRangePosition < m_DataSetRanges.GetSize())
      {
        // the data position is influenced by the data set range position and the thread id
        const math::Range< size_t> &data( m_DataSetRanges( THREAD_ID + m_DataSetRangePosition));
        for
        (
          size_t feature_id( data.GetMin()), feature_end( data.GetMax());
          feature_id < feature_end;
          ++feature_id
        )
        {
          // training step
          const size_t eff_id( m_Order( feature_id));
          const FeatureReference< float> &feature( m_TrainingData->GetFeaturesPtr()->operator()( eff_id));
          CalcHiddenTerms( feature, THREAD_ID, eff_id);
          if( CalcErrorTerms( m_TrainingData->GetResultsPtr()->operator()( eff_id), THREAD_ID, eff_id))
          {
            CalcChangeTerms( feature, THREAD_ID);
            ++m_NumberFeaturesUpdated( THREAD_ID);
          }
        }
        // TODO: Test what happens if we always backpropagate
        m_DataSelector->FinalizeConformation();
        for
        (
          size_t feature_id( data.GetMin()), feature_end( data.GetMax());
          feature_id < feature_end;
          ++feature_id
        )
        {
          // training step
          const size_t eff_id( m_Order( feature_id));
          const FeatureReference< float> &feature( m_TrainingData->GetFeaturesPtr()->operator()( eff_id));
          CalcHiddenTerms( feature, THREAD_ID, eff_id);
          if( CalcErrorTerms( m_TrainingData->GetResultsPtr()->operator()( eff_id), THREAD_ID, eff_id))
          {
            CalcChangeTerms( feature, THREAD_ID);
            ++m_NumberFeaturesUpdated( THREAD_ID);
          }
        }
        m_DataSelector->FinalizeConformation();
      }
    } // TrainThread

    void ApproximatorNeuralNetworkSelective::UpdateWeights()
    {
      // handle dropout for all hidden neurons (except the output layer)
      for
      (
        size_t layer_number( 1), output_layer_number( m_NumberToDrop.GetSize());
        layer_number < output_layer_number;
        ++layer_number
      )
      {
        if( !m_NumberToDrop( layer_number))
        {
          // nothing to drop, continue
          continue;
        }
        // shuffle the first n-indices in the m_NeuronIndices array
        storage::Vector< size_t>::const_iterator
          itr_chosen( m_ChosenDropped( layer_number).Begin()), itr_chosen_end( m_ChosenDropped( layer_number).End());
        linal::Matrix< float> &weight_matrix( m_SlopesWeight( 0)( layer_number - 1));
        linal::Vector< float> &bias_vector( m_SlopesBias( 0)( layer_number - 1));
        for( ; itr_chosen != itr_chosen_end; ++itr_chosen)
        {
          weight_matrix.GetRow( *itr_chosen) = util::GetUndefined< float>();
          bias_vector( *itr_chosen) = util::GetUndefined< float>();
        }
      }
      for
      (
        size_t layer_number( 0), output_layer_number( m_NumberToDrop.GetSize());
        layer_number < output_layer_number;
        ++layer_number
      )
      {
        if( !m_NumberToDrop( layer_number))
        {
          // nothing to drop, continue
          continue;
        }
        // shuffle the first n-indices in the m_NeuronIndices array
        storage::Vector< size_t>::const_iterator
          itr_chosen( m_ChosenDropped( layer_number).Begin()), itr_chosen_end( m_ChosenDropped( layer_number).End());
        linal::Matrix< float> &weight_matrix( m_SlopesWeight( 0)( layer_number));
        for( ; itr_chosen != itr_chosen_end; ++itr_chosen)
        {
          weight_matrix.SetCol( *itr_chosen, util::GetUndefined< float>());
        }
      }

      for
      (
        size_t layer_number( 0), output_layer_number( m_Weight.GetSize());
        layer_number < output_layer_number;
        layer_number++
      )
      {
        // update the bias using the desired propagation method
        ( *m_WeightUpdaters( layer_number))
        (
          m_Weight( layer_number).Begin(),
          m_SlopesWeight( 0)( layer_number).Begin()
        );

        // update the bias using the desired propagation method
        ( *m_BiasUpdaters( layer_number))
        (
          m_Bias( layer_number).Begin(),
          m_SlopesBias( 0)( layer_number).Begin()
        );
      }
      UpdateDroppedNeurons();
    } // UpdateWeights

    void ApproximatorNeuralNetworkSelective::UpdateDroppedNeurons()
    {
      for
      (
        size_t layer_number( 0), output_layer_number( m_Weight.GetSize());
        layer_number < output_layer_number;
        layer_number++
      )
      {
        if( !m_NumberToDrop( layer_number))
        {
          // nothing to drop, continue
          continue;
        }
        storage::Vector< size_t> &chosen_array( m_ChosenDropped( layer_number));
        const size_t nr_blocks( layer_number ? m_DropoutPartitions( layer_number - 1) + 1 : 1);
        const size_t n_to_drop_per_block( m_NumberToDrop( layer_number));
        for( size_t block_number( 0), neuron( 0); block_number < nr_blocks; ++block_number)
        {
          // shuffle the first n-indices in the m_NeuronIndices array
          storage::Vector< size_t> &indices_array( m_NeuronIndices( layer_number)( block_number));
          const size_t layer_size( indices_array.GetSize());
          for( size_t block_neuron( 0); block_neuron < n_to_drop_per_block; ++block_neuron, ++neuron)
          {
            std::swap( indices_array( block_neuron), indices_array( random::GetGlobalRandom().Random( block_neuron, layer_size - 1)));
            chosen_array( neuron) = indices_array( block_neuron);
          }
        }
        if( m_UpdateEveryNthFeature == 0 || m_UpdateEveryNthFeature > chosen_array.GetSize())
        {
          chosen_array.Sort( std::less< size_t>());
        }
      }
    } // UpdateWeights

    //! @brief sets up the neural network with a particular architecture, after training data was set
    void ApproximatorNeuralNetworkSelective::SetupArchitecture()
    {
      BCL_Assert
      (
        m_TrainingData.IsDefined(),
        "SetupArchitecture requires valid training data!"
      );

      // create the complete architecture, including input and output neurons
      // start by making a vector with just the input neurons
      storage::Vector< size_t> architecture( 1, m_TrainingData->GetFeatureSize());

      // append the hidden neurons
      if( !m_HiddenArchitecture.IsEmpty())
      {
        architecture.Append( m_HiddenArchitecture);
      }

      // add the output neurons
      architecture.PushBack( m_TrainingData->GetResultSize());

      m_Bias.Reset();
      m_Weight.Reset();
      m_Bias.AllocateMemory( architecture.GetSize() - 1);
      m_Weight.AllocateMemory( architecture.GetSize() - 1);

      if( m_TrainingData->GetSize() == 0)
      {
        // no training data, can't set the architecture up yet
        return;
      }

      // ensure that there are no empty layers
      BCL_Assert
      (
        math::Statistics::MinimumValue( architecture.Begin(), architecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );

      // initialize data for the given architecture
      for( size_t i( 1); i < architecture.GetSize(); ++i)
      {
        m_Bias.PushBack( linal::Vector< float>( architecture( i), 0.0));
        m_Weight.PushBack( linal::Matrix< float>( architecture( i), architecture( i - 1), 0.0));
      }

      // randomize the initial values of bias and weight
      storage::Vector< size_t>::const_iterator itr_arch( architecture.Begin());
      for
      (
        storage::Vector< linal::Vector< float> >::iterator itr( m_Bias.Begin()), itr_end( m_Bias.End());
        itr != itr_end;
        ++itr, ++itr_arch
      )
      {
        for( float *itr_bias( itr->Begin()), *itr_bias_end( itr->End()); itr_bias < itr_bias_end; ++itr_bias)
        {
          *itr_bias = random::GetGlobalRandom().RandomGaussian( 0.0, 0.1);
        }
      }

      itr_arch = architecture.Begin();
      for
      (
        storage::Vector< linal::Matrix< float> >::iterator itr( m_Weight.Begin()), itr_end( m_Weight.End());
        itr != itr_end;
        ++itr, ++itr_arch
      )
      {
        float std( 1.0 / math::Sqrt( float( *itr_arch)));
        BCL_MessageStd( "std: " + util::Format()( std));
        if( m_ConnectionDensity >= float( 1.0) || itr == m_Weight.Last())
        {
          for( float *itr_weight( itr->Begin()), *itr_weight_end( itr->End()); itr_weight < itr_weight_end; ++itr_weight)
          {
            *itr_weight = random::GetGlobalRandom().RandomGaussian( 0.0, std);
          }
        }
        else
        {
          for( float *itr_weight( itr->Begin()), *itr_weight_end( itr->End()); itr_weight < itr_weight_end; ++itr_weight)
          {
            if( random::GetGlobalRandom().Double() < m_ConnectionDensity)
            {
              *itr_weight = random::GetGlobalRandom().RandomGaussian( 0.0, std);
            }
          }
        }
      }
    }

    //! @brief sets up the weight updaters and change weights
    //!        this has to be done unless the iterate is read from a file
    void ApproximatorNeuralNetworkSelective::InitializeWeightUpdaters()
    {
      // create the architecture, excluding input neurons (for which there are no weight updaters)
      storage::Vector< size_t> architecture( m_HiddenArchitecture);

      // add the output neurons
      architecture.PushBack( m_TrainingData->GetResultSize());

      m_BiasUpdaters.Reset();
      m_WeightUpdaters.Reset();
      m_IterationWeightUpdaters.Reset();
      m_BiasUpdaters.AllocateMemory( architecture.GetSize() - 1);
      m_WeightUpdaters.AllocateMemory( architecture.GetSize() - 1);
      m_IterationWeightUpdaters.AllocateMemory( architecture.GetSize() - 1);

      // ensure that there are no empty layers
      BCL_Assert
      (
        math::Statistics::MinimumValue( architecture.Begin(), architecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );

      if( !m_BiasUpdateType.IsDefined())
      {
        m_BiasUpdateType = m_WeightUpdateType;
      }

      // initialize data for the given architecture
      for( size_t i( 0); i < architecture.GetSize(); ++i)
      {
        // initialize weight updaters for the given architecture
        m_BiasUpdaters.PushBack( util::CloneToShPtr( *m_BiasUpdateType));
        m_WeightUpdaters.PushBack( util::CloneToShPtr( *m_WeightUpdateType));
        m_BiasUpdaters.LastElement()->Initialize( m_Bias( i).GetSize());
        m_WeightUpdaters.LastElement()->Initialize( m_Weight( i).GetNumberOfElements());
        if( m_IterationWeightUpdateType.IsDefined())
        {
          m_IterationWeightUpdaters.PushBack( util::CloneToShPtr( *m_IterationWeightUpdateType));
        }
      }
    }

    //! @brief sets the data set ranges up for the # of threads
    void ApproximatorNeuralNetworkSelective::SetupDataSetRanges()
    {
      BCL_Assert
      (
        m_TrainingData.IsDefined(),
        "SetupDataSetRanges requires valid training data!"
      );

      size_t features_per_weight_update( m_UpdateEveryNthFeature);
      // 0 indicates batch mode, so if the m_UpdateEveryNthFeature was 0, use the data set size
      if( m_UpdateEveryNthFeature == size_t( 0))
      {
        features_per_weight_update = m_Order.GetSize();
      }

      // redetermine the # of threads, since this influences how the data is split up
      m_NumberThreads = std::min( features_per_weight_update, sched::GetNumberCPUs());

      BCL_MessageStd( "Set up data set ranges with # threads: " + util::Format()( m_NumberThreads));

      // set up the training ranges for each thread
      m_DataSetRanges.Reset();

      // calculate the # of weigh updates per run through the data set
      const size_t number_epochs_per_data_set
      (
        ( m_Order.GetSize() - 1) / features_per_weight_update + 1
      );
      features_per_weight_update = m_Order.GetSize() / number_epochs_per_data_set;
      m_DataSetRanges.AllocateMemory( number_epochs_per_data_set * m_NumberThreads);

      // if m_TrainingData->GetSize() is not a multiple of features_per_weight_update, distribute the extra features
      // over the initial number_epochs_with_extra_feature epochs
      const size_t number_epochs_with_extra_feature( m_Order.GetSize() % features_per_weight_update);

      size_t current_index( 0);

      for( size_t epoch_number( 0); epoch_number < number_epochs_per_data_set; ++epoch_number)
      {
        // determine the # of features that will be examined in this epoch
        const size_t features_in_epoch( features_per_weight_update + size_t( epoch_number < number_epochs_with_extra_feature));

        // # of features per thread per epoch
        const size_t min_features_per_thread( features_in_epoch / m_NumberThreads);

        // if features_per_weight_update is not evenly divisible by m_NumberThreads, then the extra features are distributed
        // over the initial number_threads_with_extra_feature threads
        const size_t number_threads_with_extra_feature( features_in_epoch % m_NumberThreads);

        // make data set ranges for each thread for this epoch
        for( size_t current( 0); current < m_NumberThreads; ++current)
        {
          const size_t current_width( min_features_per_thread + size_t( current < number_threads_with_extra_feature));

          m_DataSetRanges.PushBack
          (
            math::Range< size_t>
            (
              current_index,
              current_index + current_width
            )
          );
          current_index += current_width;
        }
      }

      m_RMSDError = linal::Vector< float>( m_NumberThreads, 0.0);

      // create the complete architecture, including input and output neurons
      // start by making a vector with just the input neurons
      storage::Vector< size_t> architecture( 1, m_TrainingData->GetFeatureSize());

      // append the hidden neurons
      if( !m_HiddenArchitecture.IsEmpty())
      {
        architecture.Append( m_HiddenArchitecture);
      }

      // add the output neurons
      const size_t result_size( m_TrainingData->GetResultSize());
      architecture.PushBack( result_size);

      // set up separate temporary vectors and matrices for each thread
      storage::Vector< linal::Vector< float> > temp_vector_prototype( architecture.GetSize() - 1);
      storage::Vector< linal::Matrix< float> > temp_matrix_prototype( architecture.GetSize() - 1);

      // ensure that there are no empty layers
      BCL_Assert
      (
        math::Statistics::MinimumValue( architecture.Begin(), architecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );

      // initialize data for the given architecture
      for( size_t i( 1); i < architecture.GetSize(); ++i)
      {
        temp_vector_prototype( i - 1) = linal::Vector< float>( architecture( i), 0.0);
        temp_matrix_prototype( i - 1) = linal::Matrix< float>( architecture( i), architecture( i - 1), 0.0);
      }

      // set the temporary vectors equal to each other
      if( ( !m_Dropout.IsEmpty() && size_t( m_Dropout( 0) * architecture( 0))) || m_NoiseZScore)
      {
        m_FeatureWithDropout
          = storage::Vector< linal::Vector< float> >( m_NumberThreads, linal::Vector< float>( architecture( 0), 0.0));
      }
      m_SlopesBias = m_Errors = m_HiddenInput = m_Hidden =
        storage::Vector< storage::Vector< linal::Vector< float> > >( m_NumberThreads, temp_vector_prototype);
      m_SlopesWeight =
        storage::Vector< storage::Vector< linal::Matrix< float> > >( m_NumberThreads, temp_matrix_prototype);
      m_NumberFeaturesUpdated = linal::Vector< size_t>( m_NumberThreads, size_t( 0));

      BCL_MessageStd( "Set up data set ranges with # ranges: " + util::Format()( m_DataSetRanges.GetSize()));
      BCL_MessageStd( "Set up data set ranges for updating every nth feature: " + util::Format()( features_per_weight_update));
      m_ThreadRandomNumberGenerators
        = storage::Vector< random::UniformDistribution>
          (
            m_NumberThreads,
            random::UniformDistribution( random::GetGlobalRandom().GetSeed())
          );
    } // SetupDataSetRanges

    //! @brief balance features based on the objective function
    //! @param MAX_REPEATS the maximum # of repeats allowed for any given feature; used for multi-column outputs, when
    //!        otherwise an output column with extremely few positives would see severe over-representation of certain
    //!        features
    //! @param MAX_BALANCE_RATIO 0-1; indicates the maximum ratio that the under-represented features will be balanced
    //! If underepresented features are already balanced to at least MAX_BALANCE_RATIO before calling this function,
    //! no features will be removed.  If the value of MAX_BALANCE_RATIO is greater than MAX_REPEATS would allow for the
    //! column; then each feature will just be replicated MAX_REPEATS times
    void ApproximatorNeuralNetworkSelective::BalanceFeatures
    (
      const size_t &MAX_REPEATS,
      const float &MAX_BALANCE_RATIO
    )
    {
      m_Order.Resize( m_TrainingData->GetSize());
      for( size_t i( 0), n_features( m_TrainingData->GetSize()); i < n_features; ++i)
      {
        m_Order( i) = i;
      }
      // get scaling and cutoff information
      const float cutoff( m_ObjectiveFunction->GetThreshold());
      if( !util::IsDefined( cutoff))
      {
        BCL_MessageCrt( "Cannot balance a dataset when a non-classification type objective is used!");
        return;
      }

      // determine the class of each training example result
      DetermineResultClasses();

      // create reflecting iterators for each class and determine the size of the most populated class
      const size_t n_classes( m_PeerFeatures.GetSize());
      BCL_Assert( n_classes > 1, "Only a single class of features was found; balancing is impossible!");
      size_t max_class_size( m_PeerFeatures.FirstElement().GetSize());
      storage::Vector< size_t> class_iterators( n_classes, size_t( 0));
      for( size_t class_id( 0); class_id < n_classes; ++class_id)
      {
        max_class_size = std::max( max_class_size, m_PeerFeatures( class_id).GetSize());
      }
      size_t second_most_common_class_size( 0);
      for( size_t class_id( 0); class_id < n_classes; ++class_id)
      {
        if( m_PeerFeatures( class_id).GetSize() != max_class_size)
        {
          second_most_common_class_size = std::max( second_most_common_class_size, m_PeerFeatures( class_id).GetSize());
        }
      }

      size_t max_target_size( std::min( size_t( max_class_size * m_BalanceMaxOversampling + 1), max_class_size));
      storage::Vector< size_t> max_this_class( n_classes, max_class_size);
      float oversampling_factor
      (
        std::min( float( m_BalanceMaxRepeatedFeatures), float( max_target_size) / float( second_most_common_class_size))
      );
      for( size_t class_id( 0); class_id < n_classes; ++class_id)
      {
        if( m_PeerFeatures( class_id).GetSize() < max_class_size)
        {
          max_this_class( class_id) =
            std::min( m_PeerFeatures( class_id).GetSize() * oversampling_factor, float( max_target_size));
        }
      }
      m_Order.Reset();
      m_Order.AllocateMemory( n_classes * max_class_size);
      for( size_t class_feature( 0); class_feature < max_class_size; ++class_feature)
      {
        for( size_t class_id( 0); class_id < n_classes; ++class_id)
        {
          if( class_feature < max_this_class( class_id))
          {
            const storage::Vector< size_t> &peers( m_PeerFeatures( class_id));
            m_Order.PushBack( peers( class_iterators( class_id) % peers.GetSize()));
            ++class_iterators( class_id);
          }
        }
      }
      m_Order.Shuffle();
    }

    //! @brief helper function to initialize dropout-related vectors
    void ApproximatorNeuralNetworkSelective::InitializeDropout()
    {
      // every hidden layer and the input layer can incorporate dropout
      const size_t number_potential_dropout_layers( m_HiddenArchitecture.GetSize() + 1);
      // handle dropout
      if( m_Dropout.IsEmpty())
      {
        m_Dropout = linal::Vector< float>( number_potential_dropout_layers, float( 0));
      }
      else if( m_Dropout.GetSize() > number_potential_dropout_layers)
      {
        BCL_Exit( "Too many dropout layers specified!", -1);
      }
      else if( m_Dropout.GetSize() < number_potential_dropout_layers)
      {
        BCL_MessageCrt( "Assuming layers > " + util::Format()( m_Dropout.GetSize()) + " have no dropout!");
        linal::Vector< float> tmp( number_potential_dropout_layers, float( 0));
        std::copy( m_Dropout.Begin(), m_Dropout.End(), tmp.Begin());
        m_Dropout = tmp;
      }

      // the only thing special about dropout-immune neurons is that they cannot occur in the visible layer, since
      // that really does not make any sense
      if( m_DropoutPartitions.IsEmpty())
      {
        m_DropoutPartitions = linal::Vector< size_t>( m_HiddenArchitecture.GetSize(), size_t( 0));
      }
      else if( m_DropoutPartitions.GetSize() > m_HiddenArchitecture.GetSize())
      {
        BCL_Exit( "Too many dropout-immune layers specified!", -1);
      }
      else if( m_DropoutPartitions.GetSize() < m_HiddenArchitecture.GetSize())
      {
        linal::Vector< size_t> tmp( m_HiddenArchitecture.GetSize() - 1, size_t( 0));
        std::copy( m_DropoutPartitions.Begin(), m_DropoutPartitions.End(), tmp.Begin());
        m_DropoutPartitions = tmp;
      }

      // create the complete architecture, including input and output neurons
      // start by making a vector with just the input neurons
      storage::Vector< size_t> architecture( 1, m_TrainingData->GetFeatureSize());

      // append the hidden neurons
      if( !m_HiddenArchitecture.IsEmpty())
      {
        architecture.Append( m_HiddenArchitecture);
      }

      m_NumberToDrop = linal::Vector< size_t>( number_potential_dropout_layers, size_t( 0));
      m_NumberToDrop( 0) = architecture( 0) * m_Dropout( 0);
      for( size_t layer( 1), n_layers( architecture.GetSize()); layer < n_layers; ++layer)
      {
        const size_t n_dropout_blocks_this_layer( m_DropoutPartitions( layer - 1) + 1);
        m_NumberToDrop( layer) = architecture( layer) * m_Dropout( layer) / float( n_dropout_blocks_this_layer);
        if( m_Dropout( layer) > 0.0 && m_NumberToDrop( layer) == 0)
        {
          BCL_MessageCrt
          (
            "Warning; a dropout value was specified it is too small for at least one neuron to be dropped from every "
            "partition of layer " + util::Format()( layer)
          );
        }
      }
      m_NeuronIndices = storage::Vector< storage::Vector< storage::Vector< size_t> > >( number_potential_dropout_layers);
      m_ChosenDropped = storage::Vector< storage::Vector< size_t> >( number_potential_dropout_layers);
      for( size_t layer( 0), n_layers( architecture.GetSize()); layer < n_layers; ++layer)
      {
        const size_t dropout_prone( architecture( layer));
        const size_t n_dropout_blocks( layer ? m_DropoutPartitions( layer - 1) + 1 : 1);
        m_NeuronIndices( layer).Resize( n_dropout_blocks, storage::Vector< size_t>( dropout_prone));

        m_ChosenDropped( layer).Resize( m_NumberToDrop( layer) * n_dropout_blocks);
        const size_t base_neurons_per_block( dropout_prone / n_dropout_blocks);
        const size_t nr_blocks_extra_neuron( dropout_prone % n_dropout_blocks);
        for( size_t block_nr( 0), neuron( 0); block_nr < n_dropout_blocks; ++block_nr)
        {
          const size_t nr_neurons_this_block( base_neurons_per_block + ( block_nr < nr_blocks_extra_neuron));
          for( size_t block_neuron( 0); block_neuron < nr_neurons_this_block; ++block_neuron, ++neuron)
          {
            m_NeuronIndices( layer)( block_nr)( block_neuron) = neuron;
          }
        }
      }
    }

    //! @brief determine result classes; only called if input dropout type requires peers
    void ApproximatorNeuralNetworkSelective::DetermineResultClasses()
    {
      if( !m_ResultClass.IsEmpty())
      {
        return;
      }
      const float cutoff( m_ObjectiveFunction->GetThreshold());
      const linal::MatrixConstReference< float> results( m_TrainingData->GetResultsReference());

      const size_t results_size( m_TrainingData->GetResultSize());
      // Rescaled cutoffs, one for each result
      linal::Vector< float> scaled_cutoffs( results_size, cutoff);

      // compute the rescaled cutoffs
      const util::SiPtr< const RescaleFeatureDataSet> results_scaling( m_TrainingData->GetResultsPtr()->GetScaling());
      for( size_t result( 0); result < results_size; ++result)
      {
        scaled_cutoffs( result) = results_scaling->RescaleValue( result, cutoff);
      }

      const size_t bitsize( sizeof( size_t) * 8);
      const size_t n_sizets( ( results_size - 1) / bitsize + 1);
      storage::Vector< size_t> result_class( n_sizets, size_t( 0));
      m_ResultClass = storage::Vector< size_t>( m_TrainingData->GetSize(), size_t( 0));
      // for each result, compute the result class
      size_t n_classes_so_far( 0);
      storage::Map< storage::Vector< size_t>, size_t> class_to_index;
      storage::Vector< storage::Vector< size_t> > classes;
      for( size_t feature( 0), n_features( m_TrainingData->GetSize()); feature < n_features; ++feature)
      {
        result_class.SetAllElements( 0);
        linal::VectorConstReference< float> result_row( results.GetRow( feature));
        for( size_t result_index( 0); result_index < results_size; ++result_index)
        {
          if( result_row( result_index) >= scaled_cutoffs( result_index))
          {
            result_class( result_index / bitsize) |= ( 1 << ( result_index % bitsize));
          }
        }
        storage::Map< storage::Vector< size_t>, size_t>::const_iterator itr_map( class_to_index.Find( result_class));
        if( itr_map == class_to_index.End())
        {
          class_to_index[ result_class] = m_ResultClass( feature) = n_classes_so_far++;
          classes.PushBack( result_class);
        }
        else
        {
          m_ResultClass( feature) = itr_map->second;
        }
      }
      // accumulate all classes
      m_PeerFeatures = storage::Vector< storage::Vector< size_t> >( n_classes_so_far);
      for( size_t feature( 0), n_features( m_TrainingData->GetSize()); feature < n_features; ++feature)
      {
        m_PeerFeatures( m_ResultClass( feature)).PushBack( feature);
      }
      std::stringstream info_stream;
      info_stream << "Found " << n_classes_so_far << " classes of training points. Counts per class: ";
      for( size_t class_id( 0); class_id < n_classes_so_far; ++class_id)
      {
        info_stream << m_PeerFeatures( class_id).GetSize() << ' ';
      }
      info_stream << ".  Class Binary IDs: ";
      for( size_t class_id( 0); class_id < n_classes_so_far; ++class_id)
      {
        for( size_t result_index( 0); result_index < results_size; ++result_index)
        {
          info_stream << int( bool( classes( class_id)( result_index / bitsize) & ( 1 << ( result_index % bitsize))));
        }
        info_stream << ' ';
      }
      BCL_MessageStd( info_stream.str());
    }

  } // namespace model
} // namespace bcl
