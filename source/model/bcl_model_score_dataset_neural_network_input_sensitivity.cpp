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

// include forward header of this class
#include "model/bcl_model_score_dataset_neural_network_input_sensitivity.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"
#include "model/bcl_model_feature_data_reference.h"
#include "model/bcl_model_neural_network.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetNeuralNetworkInputSensitivity::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetNeuralNetworkInputSensitivity())
    );

    // Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
    // gigantic models repetitively
    storage::Map< std::string, storage::Vector< util::OwnPtr< NeuralNetwork> > > ScoreDatasetNeuralNetworkInputSensitivity::s_Models =
      storage::Map< std::string, storage::Vector< util::OwnPtr< NeuralNetwork> > >();

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetNeuralNetworkInputSensitivity
    ScoreDatasetNeuralNetworkInputSensitivity *ScoreDatasetNeuralNetworkInputSensitivity::Clone() const
    {
      return new ScoreDatasetNeuralNetworkInputSensitivity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetNeuralNetworkInputSensitivity::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetNeuralNetworkInputSensitivity>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetNeuralNetworkInputSensitivity::GetAlias() const
    {
      static const std::string s_Name( "InputSensitivityNeuralNetwork");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief score features and results
    //! @param FEATURES matrix of features
    //! @param RESULTS matrix of results
    //! @return scores of the dataset
    linal::Vector< float> ScoreDatasetNeuralNetworkInputSensitivity::Score( const descriptor::Dataset &DATASET) const
    {
      const size_t feature_size( DATASET.GetFeatureSize());
      const size_t dataset_size( DATASET.GetSize());

      // handle empty dataset
      if( dataset_size == size_t( 0))
      {
        return linal::Vector< float>( feature_size, float( 0.0));
      }

      // get a reference to the features and results
      linal::MatrixConstReference< float> features_ref( DATASET.GetFeaturesReference());
      linal::MatrixConstReference< float> results_ref( DATASET.GetResultsReference());

      m_DatasetPtr    = util::ToSiPtr( DATASET);
      m_FeatureNumber = 0;

      // make predictions with all models
      const size_t n_models( m_Models->GetSize());
      m_PredictionClassifications = storage::Vector< linal::Matrix< char> >( n_models);
      FeatureDataReference< float> experimental_results_fds( DATASET.GetResultsReference());
      if( DATASET.GetIdsPtr().IsDefined())
      {
        m_Objective->SetData( experimental_results_fds, *DATASET.GetIdsPtr());
      }
      else
      {
        m_Objective->SetData( experimental_results_fds);
      }
      for( size_t model_number( 0); model_number < n_models; ++model_number)
      {
        // get the predictions for this model
        FeatureDataSet< float> predictions
        (
          m_Models->operator ()( model_number)->operator()
          (
            FeatureDataReference< float>( DATASET.GetFeaturesReference())
          )
        );

        // have the objective function compute the classifications
        m_PredictionClassifications( model_number)
          = m_Objective->GetFeaturePredictionClassifications( experimental_results_fds, predictions);
      }
      m_Scorer.InitializeBalancing( m_PredictionClassifications, m_DatasetPtr->GetFeatureSize());

      // determine the # of threads to use
      size_t n_threads( std::max( size_t( 1), sched::GetScheduler().GetNumberUnusedCPUS()));

      m_AveDescriptorScores =
        storage::Vector< math::RunningAverage< linal::Vector< float> > >
        (
          n_threads,
          math::RunningAverage< linal::Vector< float> >( linal::Vector< float>( feature_size, 0.0))
        );

      // initially, just predict all results
      m_FeaturesToComputeSize = 0;

      util::ShPtrVector< sched::JobInterface> jobs( n_threads);
      storage::Vector< size_t> thread_ids( n_threads, size_t( 0));

      // initialize data for each thread
      for( size_t thread_number( 0); thread_number < n_threads; ++thread_number)
      {
        thread_ids( thread_number) = thread_number;
        jobs( thread_number) =
          util::ShPtr< sched::JobInterface>
          (
            new sched::UnaryFunctionJobWithData< const size_t, void, ScoreDatasetNeuralNetworkInputSensitivity>
            (
              0,
              *this,
              &ScoreDatasetNeuralNetworkInputSensitivity::RunThread,
              thread_ids( thread_number),
              sched::JobInterface::e_READY,
              NULL
            )
          );
        sched::GetScheduler().SubmitJob( jobs( thread_number));
      }
      sched::GetScheduler().Join( jobs( 0));
      for( size_t thread_number( 1); thread_number < n_threads; ++thread_number)
      {
        sched::GetScheduler().Join( jobs( thread_number));
        if( m_AveDescriptorScores( thread_number).GetWeight())
        {
          if( !m_AveDescriptorScores( 0).GetWeight())
          {
            m_AveDescriptorScores( 0) = m_AveDescriptorScores( thread_number);
          }
          else
          {
            m_AveDescriptorScores( 0).AddWeightedObservation
            (
              m_AveDescriptorScores( thread_number).GetAverage(),
              m_AveDescriptorScores( thread_number).GetWeight()
            );
          }
        }
      }
      m_Scorer.AddUtilityScore( m_AveDescriptorScores( 0));

      // return input-sensitivity for every feature column in dataset
      return m_AveDescriptorScores( 0).GetAverage();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool ScoreDatasetNeuralNetworkInputSensitivity::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // get the models if they are not already in memory
      storage::Vector< util::OwnPtr< NeuralNetwork> > &models( s_Models[ GetString()]);

      if( m_Key.empty())
      {
        // get all the keys for this storage
        storage::Vector< std::string> keys( m_ModelStorage->GetAllKeys());

        BCL_MessageStd( "Found " + util::Format()( keys.GetSize()) + " models for input sensitivity");

        if( models.IsEmpty() && !keys.IsEmpty())
        {
          // no models already loaded
          RetrieveInterface::t_Container interfaces( m_ModelStorage->RetrieveEnsemble( keys));
          models = storage::Vector< util::OwnPtr< NeuralNetwork> >( interfaces.Begin(), interfaces.End());
        }
        if( keys.GetSize() != models.GetSize())
        {
          ERR_STREAM << "# of models in " << m_ModelStorage->GetString() << " changed; aborting";
          return false;
        }
        else if( keys.IsEmpty())
        {
          ERR_STREAM << "No models found in storage! " << m_ModelStorage->GetString();
          return false;
        }
      }
      else if( models.IsEmpty())
      {
        // no models already loaded
        models.PushBack( m_ModelStorage->Retrieve( m_Key));
        if( !models.LastElement().IsDefined())
        {
          ERR_STREAM << m_Key << " is not a key for " << m_ModelStorage->GetString();
          return false;
        }
      }
      if( !m_Objective.TryRead( m_ModelStorage->RetrieveEnsembleCVInfo().FirstElement().GetObjective(), ERR_STREAM))
      {
        ERR_STREAM << "Could not read objective function, aborting";
        return false;
      }

      m_Models = util::ToSiPtr( models);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetNeuralNetworkInputSensitivity::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Calculates the sensitivity of the model to a change in inputs");

      parameters.AddInitializer
      (
        "storage",
        "type of storage for models",
        io::Serialization::GetAgent( &m_ModelStorage)
      );

      parameters.AddInitializer
      (
        "key",
        "optional key of the desired model in the storage",
        io::Serialization::GetAgent( &m_Key),
        ""
      );

      parameters.AddInitializer
      (
        "weights",
        "Change the weighting of various measures the derivative scores",
        io::Serialization::GetAgent( &m_Scorer)
      );

      return parameters;
    }

    //! @brief determine the performance of the model ensemble on a particular feature
    //! @param PREDICTION_CLASSIFICATIONS vector of model prediction classifications
    //! @param FEATURE_NR feature number of interest
    //! @return a vector containing strings with P|p|N|n|\0 for whether each result is a TP, FP, TN, FN< or NA for each model
    storage::Vector< std::string> ScoreDatasetNeuralNetworkInputSensitivity::PartitionModels
    (
      const storage::Vector< linal::Matrix< char> > &PREDICTION_CLASSIFICATIONS,
      const size_t &FEATURE_NR
    )
    {
      const size_t n_models( PREDICTION_CLASSIFICATIONS.GetSize());
      const size_t result_size( PREDICTION_CLASSIFICATIONS( 0).GetNumberCols());
      storage::Vector< std::string> partition( n_models, std::string( result_size, '\0'));
      // transfer model correctness into the desired format
      for( size_t model_number( 0); model_number < n_models; ++model_number)
      {
        for( size_t result_number( 0); result_number < result_size; ++result_number)
        {
          partition( model_number)[ result_number] = PREDICTION_CLASSIFICATIONS( model_number)( FEATURE_NR, result_number);
        }
      }
      return partition;
    }

    //! @brief determine the side of the actual result
    //! @param ACTUAL actual / experimental value
    //! @param CUTOFF the cutoff value
    //! @return a vector containing size_ts 0/1/undefined depending on whether the model was bad, good, or unknown
    storage::Vector< size_t> ScoreDatasetNeuralNetworkInputSensitivity::GetCutoffSides
    (
      const linal::VectorConstInterface< float> &ACTUAL,
      const float &CUTOFF
    )
    {
      const size_t results_size( ACTUAL.GetSize());
      if( !util::IsDefined( CUTOFF))
      {
        return storage::Vector< size_t>( results_size, size_t( 0));
      }
      storage::Vector< size_t> cutoff_side( results_size, size_t( 0));
      for( size_t i( 0); i < results_size; ++i)
      {
        if( ACTUAL( i) >= CUTOFF)
        {
          cutoff_side( i) = 1;
        }
      }

      return cutoff_side;
    }

    //! @brief Get the next feature for a thread to make predictions for
    //! @param LAST_FEATURE last feature that this thread computed
    //! @return the next feature for a thread to make predictions for
    size_t ScoreDatasetNeuralNetworkInputSensitivity::GetNextFeatureForPrediction( const size_t &LAST_FEATURE) const
    {
      m_FeatureNumberMutex.Lock();
      const size_t dataset_size( m_DatasetPtr->GetSize());

      if( util::IsDefined( LAST_FEATURE))
      {
        ++m_FeaturesToComputeSize;
      }

      // only write out status from the first thread
      // determine progress percent
      const size_t percent( float( m_FeaturesToComputeSize) * 100.0f / dataset_size);

      // determine number of stars in the status bar
      const size_t number_stars( percent / 5);

      const std::string status
      (
        "["
        + std::string( number_stars, '*')
        + std::string( 20 - number_stars, ' ')
        + "] "
        + util::Format()( percent) + "% "
        + util::Format()( m_FeaturesToComputeSize) + " / " + util::Format()( size_t( dataset_size))
        + " features predicted with"
      );

      if( m_FeatureNumber == dataset_size)
      {
        m_FeatureNumberMutex.Unlock();
        return dataset_size;
      }

      const size_t number_to_return( m_FeatureNumber);

      ++m_FeatureNumber;
      util::GetLogger().LogStatus( status);
      m_FeatureNumberMutex.Unlock();
      return number_to_return;
    }

    //! @brief run a thread to compute the kernel for all features with ID = THREAD_ID % n_threads
    //! @param THREAD_ID id of the thread to run (0-indexed)
    void ScoreDatasetNeuralNetworkInputSensitivity::RunThread( const size_t &THREAD_ID) const
    {
      const size_t feature_size( m_DatasetPtr->GetFeatureSize());
      const size_t result_size( m_DatasetPtr->GetResultSize());
      const size_t dataset_size( m_DatasetPtr->GetSize());

      const storage::Vector< util::OwnPtr< NeuralNetwork> > &models( *m_Models);
      const size_t n_models( models.GetSize());

      // get a reference to the features and results
      linal::MatrixConstReference< float> features_ref( m_DatasetPtr->GetFeaturesReference());

      // store all weight effect matrices
      storage::Vector< linal::Matrix< float> > weight_effects
      (
        n_models,
        linal::Matrix< float>( feature_size, result_size)
      );

      linal::MatrixConstReference< float> dataset_reference( m_DatasetPtr->GetFeaturesReference());
      linal::MatrixConstReference< float> results_reference( m_DatasetPtr->GetResultsReference());

      linal::Matrix< float> model_predictions( result_size, n_models);

      // make predictions.  Threads work together to generate all predictions
      for
      (
        size_t feature_number( GetNextFeatureForPrediction( util::GetUndefined< size_t>()));
        feature_number < dataset_size;
        feature_number = GetNextFeatureForPrediction( feature_number)
      )
      {
        // create a feature data set containing just that feature
        linal::Matrix< float> features_matrix( size_t( 1), feature_size, dataset_reference[ feature_number]);
        FeatureDataSet< float> original_feature( features_matrix);
        // for each model
        for( size_t model_number( 0); model_number < n_models; ++model_number)
        {
          // get a reference to the model
          const NeuralNetwork &model( *models( model_number));
          linal::Matrix< float> &model_weight_effects( weight_effects( model_number));
          original_feature.DeScale();
          original_feature.Rescale( *model.GetRescaleInput());

          // get the rescaled vector
          linal::VectorReference< float> original_feature_vector( original_feature.GetRawMatrix().GetRow( 0));

          // get the result and the input sensitivity
          storage::Pair< linal::Vector< float>, linal::Matrix< float> >
            result_input_sensitivity( model.ComputeResultInputSensitivity( original_feature_vector));

          // compute the original results with this model
          const FeatureDataSet< float> original_results_fds( model( original_feature));
          const linal::VectorConstReference< float> results( result_input_sensitivity.First());
          for( size_t result_number( 0); result_number < result_size; ++result_number)
          {
            model_predictions( result_number, model_number) = results( result_number);
            model_weight_effects = result_input_sensitivity.Second();
          }
        }

        m_Scorer.Score
        (
          weight_effects,
          PartitionModels( m_PredictionClassifications, feature_number),
          m_AveDescriptorScores( THREAD_ID)
        );
      }
    }

  } // namespace model
} // namespace bcl

