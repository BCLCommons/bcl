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
#include "model/bcl_model_approximator_restricted_boltzmann_machine.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_statistics.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ApproximatorRestrictedBoltzmannMachine::s_IterateInstance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorRestrictedBoltzmannMachine( false))
    );
    const util::SiPtr< const util::ObjectInterface> ApproximatorRestrictedBoltzmannMachine::s_PretrainInstance
    (
      util::Enumerated< PretrainNeuralNetworkInterface>::AddInstance( new ApproximatorRestrictedBoltzmannMachine( true))
    );

    //! @brief default constructor
    //! @param PRETRAIN whether this object is being constructed as a pretrainer or an iterate
    ApproximatorRestrictedBoltzmannMachine::ApproximatorRestrictedBoltzmannMachine( const bool &PRETRAIN) :
      m_UpdateEveryNthFeature( 0),
      m_NumberThreads(),
      m_IterationsPerRMSDMessage( 1),
      m_Type( RestrictedBoltzmannMachineLayer::e_StochasticSigmoid),
      m_Shuffle( false),
      m_NumberStochasticSteps( 3),
      m_DataSetRangePosition( 0),
      m_IsPretrainer( PRETRAIN),
      m_NumIterations( 0)
    {
    }

    //! @brief constructor from training data, transfer, rescale, and objective functions
    //! @param TRAINING_DATA data to train the NeuralNetwork on
    //! @param UPDATE_EVERY_NTH_FEATURE how often the weights get updated
    //! @param ARCHITECTURE the # of neurons in each hidden layer of the network
    //! @param TRANSFER_FUNCTION ShPtr to the transfer function between input and output of each neuron
    //! @param OBJECTIVE_FUNCTION ShPtr to objective function
    //! @param WEIGHT_UPDATE_FUNCTION method by which to update the weights
    //! @param BIAS_UPDATE_FUNCTION method by which to update the biases
    //! @param ITERATIONS_PER_RMSD_REPORT # iterations per report of the rmsd
    ApproximatorRestrictedBoltzmannMachine::ApproximatorRestrictedBoltzmannMachine
    (
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t UPDATE_EVERY_NTH_FEATURE,
      const storage::Vector< size_t> &ARCHITECTURE,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
      const util::Implementation< NeuralNetworkUpdateWeightsInterface> &WEIGHT_UPDATE_FUNCTION,
      const util::Implementation< NeuralNetworkUpdateWeightsInterface> &BIAS_UPDATE_FUNCTION,
      const size_t &ITERATIONS_PER_RMSD_REPORT,
      const RestrictedBoltzmannMachineLayer::Type &TYPE
    ) :
      ApproximatorBase( OBJECTIVE_FUNCTION),
      m_UpdateEveryNthFeature( UPDATE_EVERY_NTH_FEATURE),
      m_IterationsPerRMSDMessage( ITERATIONS_PER_RMSD_REPORT),
      m_Type( TYPE),
      m_Shuffle( UPDATE_EVERY_NTH_FEATURE), // Always shuffle except in batch-update mode
      m_NumberStochasticSteps( 3),
      m_DataSetRangePosition( 0),
      m_HiddenArchitecture( ARCHITECTURE),
      m_WeightUpdateType( WEIGHT_UPDATE_FUNCTION),
      m_BiasUpdateType( BIAS_UPDATE_FUNCTION),
      m_WeightDecay( 0.0001),
      m_IsPretrainer( false),
      m_NumIterations( 0)
    {
      // set and rescale training data set
      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorRestrictedBoltzmannMachine copied from this instance
    ApproximatorRestrictedBoltzmannMachine *ApproximatorRestrictedBoltzmannMachine::Clone() const
    {
      return new ApproximatorRestrictedBoltzmannMachine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorRestrictedBoltzmannMachine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorRestrictedBoltzmannMachine::GetAlias() const
    {
      static const std::string s_Name( "RestrictedBoltzmannMachine");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorRestrictedBoltzmannMachine::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      m_TrainingData = DATA;
      DATA->GetFeatures().DeScale();
      DATA->GetFeatures().Rescale( RestrictedBoltzmannMachineLayer::s_DefaultInputRange, RescaleFeatureDataSet::e_MinMax);
      DATA->GetResults().DeScale();
      DATA->GetResults().Rescale( TransferSigmoid().GetDynamicOutputRange());

      BCL_MessageStd
      (
        "Setting up training data with " + util::Format()( DATA->GetSize()) + " points"
      );

      m_Order.Resize( DATA->GetSize());
      for( size_t i( 0), n_features( DATA->GetSize()); i < n_features; ++i)
      {
        m_Order( i) = i;
      }

      // set the data ranges up for threading (if applicable)
      SetupDataSetRanges();

      // if the current instance is NOT used through IterateInterfaceFromFile
      // then initialize the Iterate properly
      if( !IsTrainingContinued())
      {
        // weight updaters need to be initialized regardless
        InitializeWeightUpdaters();
      }
    } // SetTrainingData

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorRestrictedBoltzmannMachine::GetCurrentModel() const
    {
      // make a new model out of the current data members
      const size_t number_layers( m_Layers.GetSize()), number_vis_layers_no_output( number_layers - 1);
      storage::Vector< linal::Vector< float> > hidden_biases( number_layers + 1);
      storage::Vector< linal::Matrix< float> > weights( number_layers + 1);
      for( size_t layer( 0); layer < number_vis_layers_no_output; ++layer)
      {
        hidden_biases( layer) = m_Layers( layer).GetHiddenBias();
        weights( layer) = m_Layers( layer).GetWeight().Transposed();
      }
      hidden_biases( number_vis_layers_no_output) = m_Layers( number_vis_layers_no_output).GetHiddenBias();

      // get the number of real hidden neurons in the next-to-last layer
      const size_t number_real_visible_neurons
      (
        m_HiddenArchitecture.GetSize() > size_t( 1)
        ? m_HiddenArchitecture( m_HiddenArchitecture.GetSize() - 2)
        : m_TrainingData->GetFeatureSize()
      );

      // number of final hidden neurons
      const size_t number_final_hidden_neurons( m_HiddenArchitecture.LastElement());

      // determine the number of output neurons
      const size_t number_output_neurons( m_TrainingData->GetResultSize());

      // get the biases from the final visible layer for the output neurons
      hidden_biases( number_layers)
        = linal::VectorConstReference< float>
          (
            number_output_neurons,
            m_Layers( number_vis_layers_no_output).GetVisibleBias().Begin() + number_real_visible_neurons
          );

      // create a reference to the weights that are hidden vs.
      weights( number_vis_layers_no_output)
        = linal::MatrixConstReference< float>
          (
            number_real_visible_neurons,
            number_final_hidden_neurons,
            m_Layers( number_vis_layers_no_output).GetWeight().Begin()
          ).Transposed();

      weights( number_layers)
        = linal::MatrixConstReference< float>
          (
            number_output_neurons,
            number_final_hidden_neurons,
            m_Layers( number_vis_layers_no_output).GetWeight()[ number_real_visible_neurons]
          );

      util::ShPtr< Interface> model
      (
        new NeuralNetwork
        (
          GetRescaleFeatureDataSet(),
          GetRescaleResultDataSet(),
          hidden_biases,
          weights,
          util::Implementation< TransferFunctionInterface>( TransferSigmoid())
        )
      );
      return model;
    }

    //! @brief the main operation, pretrains a neural network
    //! @param DATA the data for use in pretraining
    //! @param OBJECTIVE ShPtr to the objective function for the network
    util::ShPtr< NeuralNetwork> ApproximatorRestrictedBoltzmannMachine::PretrainNetwork
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

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> >
      ApproximatorRestrictedBoltzmannMachine::GetCurrentApproximation() const
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

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorRestrictedBoltzmannMachine::Next()
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
                ApproximatorRestrictedBoltzmannMachine
              >
              (
                group_id,
                *this,
                &ApproximatorRestrictedBoltzmannMachine::TrainThread,
                thread_ids( thread_number),
                sched::JobInterface::e_READY,
                NULL
              )
            );
      }

      m_DataSetRangePosition = 0;
      while( m_DataSetRangePosition < m_DataSetRanges.GetSize())
      {
        for( size_t thread_number( 0); thread_number < m_NumberThreads; thread_number++)
        {
          sched::GetScheduler().RunJob( all_jobs( thread_number));
        }

        sched::GetScheduler().Join( all_jobs( 0));
        for( size_t thread_number( 1); thread_number < m_NumberThreads; thread_number++)
        {
          sched::GetScheduler().Join( all_jobs( thread_number));
          // accumulate the slopes from the newly joined jobs with the slopes from the earlier jobs
          for
          (
            size_t hidden_layer_number( 0), number_hidden_layers( m_Trainers( thread_number).GetSize());
            hidden_layer_number < number_hidden_layers;
            ++hidden_layer_number
          )
          {
            m_Trainers( 0)( hidden_layer_number).AccumulateChangesFrom( m_Trainers( thread_number)( hidden_layer_number));
          }
        }

        UpdateWeights();

        // move on the the next epoch
        m_DataSetRangePosition += m_NumberThreads;
      }

      if( ( this->GetTracker().GetIteration() % m_IterationsPerRMSDMessage) == size_t( 0))
      {
        for( size_t thread_number( 1); thread_number < m_NumberThreads; thread_number++)
        {
          // accumulate errors
          m_RMSDError( 0) += m_RMSDError( thread_number);
          m_ReconstructionError( 0) += m_ReconstructionError( thread_number);
        }
        m_RMSDError( 0) /= m_TrainingData->GetSize();
        m_RMSDError( 0) = std::max( m_RMSDError( 0), float( 0.0));
        m_ReconstructionError( 0) /= m_TrainingData->GetSize();
        BCL_MessageStd
        (
          "Training relative RMSD for iteration " + util::Format()( this->GetTracker().GetIteration())
          + ": " + util::Format()( m_RMSDError( 0))
          + " reconstruction error in 1st layer: " + util::Format()( m_ReconstructionError( 0))
        );
      }

      m_DataSetRangePosition = 0;

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

      // get current model
      util::ShPtr< Interface> current_model( GetCurrentModel());

      const float objective_result( m_ObjectiveFunction->operator()( current_model));

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
    io::Serializer ApproximatorRestrictedBoltzmannMachine::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "trains a neural network (see http://en.wikipedia.org/wiki/Artificial_neural_network)"
      );

      parameters.AddInitializer
      (
        "type",
        "Type of RBM to train; StochasticSigmoid is a noisy sigmoid approx for continous data, see "
        "http://www.ee.nthu.edu.tw/~hchen/pubs/iee2003.pdf . Binomial is from http://www.cs.toronto.edu/~hinton/science.pdf",
        io::Serialization::GetAgent( &m_Type),
        "StochasticSigmoid"
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
      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each batch step",
        io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      parameters.AddInitializer
      (
        "steps per update",
        "# of features seen between each update of weights (set to 0 to use the size of the training data set)",
        io::Serialization::GetAgent( &m_UpdateEveryNthFeature),
        "0"
      );
      parameters.AddInitializer
      (
        "hidden architecture",
        "# of neurons in each hidden layer, e.g. (100) declares that there is 1 hidden layer w/ 100 neurons",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_HiddenArchitecture,
          io::Serialization::GetAgentWithMin( size_t( 1))
        )
      );
      parameters.AddInitializer
      (
        "rmsd report frequency",
        "# of iterations between reports of rmsd on the last training batch",
        io::Serialization::GetAgentWithMin( &m_IterationsPerRMSDMessage, size_t( 1)),
        "1"
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
        "stochastic steps",
        "0 gives fast convergence but poor generalizability, 3 gives the most generality but can result in slow "
        "convergence, especially for binomial networks.",
        io::Serialization::GetAgentWithRange( &m_NumberStochasticSteps, 0, 3),
        "3"
      );
      parameters.AddInitializer
      (
        "weight decay",
        "Weight decay parameter; large values (say 0.001) may result in better networks but slower training times"
        ". This number should be scaled with # steps per update (batch size); so if multiplying batch size by X, multiply "
        "this number by X too for consistency",
        io::Serialization::GetAgentWithRange( &m_WeightDecay, 0.0, 0.1),
        "0.00001"
      );
      if( m_IsPretrainer)
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
    std::istream &ApproximatorRestrictedBoltzmannMachine::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_IsPretrainer, ISTREAM);
      util::ObjectDataLabel label;
      io::Serialize::Read( label, ISTREAM);
      BCL_Assert( ApproximatorBase::TryRead( label, util::GetLogger()), "Could not read iterate");

      // read members
      io::Serialize::Read( m_Layers,               ISTREAM);
      io::Serialize::Read( m_Trainers,             ISTREAM);
      io::Serialize::Read( m_HiddenArchitecture,   ISTREAM);
      io::Serialize::Read( m_HiddenBiasUpdaters,   ISTREAM);
      io::Serialize::Read( m_VisibleBiasUpdaters,  ISTREAM);
      io::Serialize::Read( m_HiddenNoiseUpdaters,  ISTREAM);
      io::Serialize::Read( m_VisibleNoiseUpdaters, ISTREAM);
      io::Serialize::Read( m_WeightUpdaters,       ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorRestrictedBoltzmannMachine::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members.  Note that this is probably insufficient;
      io::Serialize::Write( m_IsPretrainer, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( ApproximatorBase::GetLabel(), OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Layers, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Trainers, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HiddenArchitecture, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HiddenBiasUpdaters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_VisibleBiasUpdaters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HiddenNoiseUpdaters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_VisibleNoiseUpdaters, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WeightUpdaters, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

    //! train ANN with a feature
    void ApproximatorRestrictedBoltzmannMachine::TrainThread( const size_t &THREAD_ID)
    {
      // get a reference to the array of trainers
      storage::Vector< TrainRestrictedBoltzmannMachineLayer> &trainers( m_Trainers( THREAD_ID));

      // determine the # of layers
      const size_t number_layers( m_Layers.GetSize());

      // reset the trainers to prepare for the next round
      for( size_t layer_number( 0); layer_number < number_layers; ++layer_number)
      {
        trainers( layer_number).Reset();
      }

      // only run the data for this range if it actually existed
      // it may not if m_DataSetRanges.GetSize() is not a multiple of the number of threads
      if( THREAD_ID + m_DataSetRangePosition >= m_DataSetRanges.GetSize())
      {
        return;
      }

      // determine the # of layers prior to the output layer
      const size_t number_layers_before_output( number_layers - 1);

      // get the feature matrix
      const linal::MatrixConstReference< float> features( m_TrainingData->GetFeaturesPtr()->GetMatrix());
      const linal::MatrixConstReference< float> results( m_TrainingData->GetResultsPtr()->GetMatrix());

      // create temporary vectors
      linal::Vector< float> first_hidden_tmp( m_HiddenArchitecture( 0));
      linal::Vector< float> first_visible_tmp( m_TrainingData->GetFeatureSize());

      // the data position is influenced by the data set range position and the thread id
      const math::Range< size_t> &data( m_DataSetRanges( THREAD_ID + m_DataSetRangePosition));
      float &rmsd_location( m_RMSDError( THREAD_ID));
      float &reconstruction_location( m_ReconstructionError( THREAD_ID));
      for
      (
        size_t feature_id( data.GetMin()), feature_end( data.GetMax());
        feature_id < feature_end;
        ++feature_id
      )
      {
        // training step
        const size_t eff_id( m_Order( feature_id));
        linal::VectorConstReference< float> feature( m_TrainingData->GetFeaturesPtr()->operator()( eff_id));

        // get the 1st layer error
        reconstruction_location += m_Layers( 0).ComputeReconstructionError( feature, first_hidden_tmp, first_visible_tmp);

        linal::VectorConstReference< float> result( results.GetRow( eff_id));
        const linal::VectorConstReference< float> original_result( result);

        // propagate the feature through the various layers
        for( size_t layer_number( 0); layer_number < number_layers_before_output; ++layer_number)
        {
          feature = trainers( layer_number).Train( feature);
        }
        // propagate through to the result
        result = trainers( number_layers_before_output).Train( feature, result);

        // track the RMSD between the reconstructed points
        float err_tracker( 0);
        // add the reconstruction error
        for
        (
          const float *itr_reconstructed( result.Begin()), *itr_original( original_result.Begin()), *itr_reconstructed_end( result.End());
          itr_reconstructed != itr_reconstructed_end;
          ++itr_reconstructed, ++itr_original
        )
        {
          err_tracker += math::Sqr( *itr_reconstructed - *itr_original);
        }
        rmsd_location += math::Sqrt( err_tracker / float( result.GetSize()));
      }
    } // TrainThread

    void ApproximatorRestrictedBoltzmannMachine::UpdateWeights()
    {
      for( size_t layer_number( 0), number_layers( m_Layers.GetSize()); layer_number < number_layers; ++layer_number)
      {
        // Update all features
        m_Trainers( 0)( layer_number).UpdateLayer
        (
          m_WeightDecay,
          *m_WeightUpdaters( layer_number),
          *m_VisibleBiasUpdaters( layer_number),
          *m_HiddenBiasUpdaters( layer_number),
          *m_VisibleNoiseUpdaters( layer_number),
          *m_HiddenNoiseUpdaters( layer_number)
        );
      }
    } // UpdateWeights

    //! @brief sets up the weight updaters and change weights
    //!        this has to be done unless the iterate is read from a file
    void ApproximatorRestrictedBoltzmannMachine::InitializeWeightUpdaters()
    {
      // create the architecture, excluding input neurons (for which there are no weight updaters)
      storage::Vector< size_t> visible_architecture( size_t( 1), m_TrainingData->GetFeatureSize());
      visible_architecture.Append( m_HiddenArchitecture);
      visible_architecture.PopBack();
      visible_architecture.LastElement() += m_TrainingData->GetResultSize();

      m_HiddenBiasUpdaters.Reset();
      m_HiddenNoiseUpdaters.Reset();
      m_WeightUpdaters.Reset();
      m_VisibleBiasUpdaters.Reset();
      m_VisibleNoiseUpdaters.Reset();
      const size_t number_layers( m_HiddenArchitecture.GetSize());
      m_HiddenBiasUpdaters.AllocateMemory( number_layers);
      m_HiddenNoiseUpdaters.AllocateMemory( number_layers);
      m_VisibleBiasUpdaters.AllocateMemory( number_layers);
      m_VisibleNoiseUpdaters.AllocateMemory( number_layers);
      m_WeightUpdaters.AllocateMemory( number_layers);

      // ensure that there are no empty layers
      BCL_Assert
      (
        math::Statistics::MinimumValue( visible_architecture.Begin(), visible_architecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );
      BCL_Assert
      (
        math::Statistics::MinimumValue( m_HiddenArchitecture.Begin(), m_HiddenArchitecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );

      if( !m_BiasUpdateType.IsDefined())
      {
        m_BiasUpdateType = m_WeightUpdateType;
      }

      // initialize data for the given architecture
      for( size_t i( 0); i < number_layers; ++i)
      {
        // initialize weight updaters for the given architecture
        m_VisibleBiasUpdaters.PushBack( util::CloneToShPtr( *m_BiasUpdateType));
        m_HiddenBiasUpdaters.PushBack( util::CloneToShPtr( *m_BiasUpdateType));
        m_WeightUpdaters.PushBack( util::CloneToShPtr( *m_WeightUpdateType));
        m_VisibleBiasUpdaters.LastElement()->Initialize( visible_architecture( i));
        m_HiddenBiasUpdaters.LastElement()->Initialize( m_HiddenArchitecture( i));
        m_WeightUpdaters.LastElement()->Initialize( visible_architecture( i) * m_HiddenArchitecture( i));
        if( m_Type == RestrictedBoltzmannMachineLayer::e_StochasticSigmoid)
        {
          m_VisibleNoiseUpdaters.PushBack( util::CloneToShPtr( *m_BiasUpdateType));
          m_HiddenNoiseUpdaters.PushBack( util::CloneToShPtr( *m_BiasUpdateType));
          m_VisibleNoiseUpdaters.LastElement()->Initialize( visible_architecture( i));
          m_HiddenNoiseUpdaters.LastElement()->Initialize( m_HiddenArchitecture( i));
        }
      }
      // for binomial types, it is still necessary to have updaters for the noise, even though they are never used
      if( m_Type != RestrictedBoltzmannMachineLayer::e_StochasticSigmoid)
      {
        m_VisibleNoiseUpdaters = m_VisibleBiasUpdaters;
        m_HiddenNoiseUpdaters = m_HiddenBiasUpdaters;
      }
    }

    //! @brief sets the data set ranges up for the # of threads
    void ApproximatorRestrictedBoltzmannMachine::SetupDataSetRanges()
    {
      BCL_Assert
      (
        m_TrainingData.IsDefined(),
        "SetupDataSetRanges requires valid training data!"
      );

      size_t features_per_weight_update( std::min( m_TrainingData->GetSize(), m_UpdateEveryNthFeature));
      // 0 indicates batch mode, so if the m_UpdateEveryNthFeature was 0, use the data set size
      if( m_UpdateEveryNthFeature == size_t( 0))
      {
        features_per_weight_update = m_TrainingData->GetSize();
      }

      // redetermine the # of threads, since this influences how the data is split up
      m_NumberThreads = std::min( features_per_weight_update, sched::GetNumberCPUs());

      BCL_MessageStd( "Set up data set ranges with # threads: " + util::Format()( m_NumberThreads));

      // set up the training ranges for each thread
      m_DataSetRanges.Reset();

      // calculate the # of weigh updates per run through the data set
      const size_t number_epochs_per_data_set
      (
        ( m_TrainingData->GetSize() - 1) / features_per_weight_update + 1
      );
      features_per_weight_update = m_TrainingData->GetSize() / number_epochs_per_data_set;
      m_DataSetRanges.AllocateMemory( number_epochs_per_data_set * m_NumberThreads);

      // if m_TrainingData->GetSize() is not a multiple of features_per_weight_update, distribute the extra features
      // over the initial number_epochs_with_extra_feature epochs
      const size_t number_epochs_with_extra_feature( m_TrainingData->GetSize() % features_per_weight_update);

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
      m_ReconstructionError = linal::Vector< float>( m_NumberThreads, 0.0);

      // create the complete architecture, including input and output neurons
      // start by making a vector with just the input neurons
      storage::Vector< size_t> architecture( 1, m_TrainingData->GetFeatureSize());

      // append the hidden neurons
      if( !m_HiddenArchitecture.IsEmpty())
      {
        architecture.Append( m_HiddenArchitecture);
      }

      // create the architecture, excluding input neurons (for which there are no weight updaters)
      storage::Vector< size_t> visible_architecture( size_t( 1), m_TrainingData->GetFeatureSize());
      visible_architecture.Append( m_HiddenArchitecture);
      visible_architecture.PopBack();
      visible_architecture.LastElement() += m_TrainingData->GetResultSize();

      // ensure that there are no empty layers
      BCL_Assert
      (
        math::Statistics::MinimumValue( visible_architecture.Begin(), visible_architecture.End()) > 0
        && math::Statistics::MinimumValue( m_HiddenArchitecture.Begin(), m_HiddenArchitecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );

      // initialize data for the given architecture
      m_Layers.Reset();
      const size_t number_layers( visible_architecture.GetSize());
      m_Layers.Resize( number_layers);
      m_Trainers.Resize( m_NumberThreads);
      for( size_t thread( 0); thread < m_NumberThreads; ++thread)
      {
        m_Trainers( thread).Reset();
        m_Trainers( thread).Resize( number_layers);
      }
      for( size_t i( 0); i < number_layers; ++i)
      {
        m_Layers( i) = RestrictedBoltzmannMachineLayer( visible_architecture( i), m_HiddenArchitecture( i), m_Type);
        for( size_t thread( 0); thread < m_NumberThreads; ++thread)
        {
          m_Trainers( thread)( i) = TrainRestrictedBoltzmannMachineLayer( m_Layers( i), m_NumberStochasticSteps);
        }
      }

      BCL_MessageStd( "Set up data set ranges with # ranges: " + util::Format()( m_DataSetRanges.GetSize()));
      BCL_MessageStd( "Set up data set ranges for updating every nth feature: " + util::Format()( features_per_weight_update));
    } // SetupDataSetRanges

  } // namespace model
} // namespace bcl
