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
#include "model/bcl_model_score_dataset_input_sensitivity.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_feature_data_reference.h"
#include "model/bcl_model_score_dataset_neural_network_input_sensitivity.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_unary_function_job_with_data.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_time.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetInputSensitivity::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetInputSensitivity())
    );

    // Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
    // gigantic models repetitively
    storage::Map< std::string, RetrieveInterface::t_Container> ScoreDatasetInputSensitivity::s_Models =
      storage::Map< std::string, RetrieveInterface::t_Container>();

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetInputSensitivity
    ScoreDatasetInputSensitivity *ScoreDatasetInputSensitivity::Clone() const
    {
      return new ScoreDatasetInputSensitivity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetInputSensitivity::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetInputSensitivity>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetInputSensitivity::GetAlias() const
    {
      static const std::string s_Name( "InputSensitivity");
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
    linal::Vector< float> ScoreDatasetInputSensitivity::Score( const descriptor::Dataset &DATASET) const
    {
      const size_t feature_size( DATASET.GetFeatureSize());
      const size_t result_size( DATASET.GetResultSize());
      const size_t dataset_size( DATASET.GetSize());

      // handle empty dataset
      if( dataset_size == size_t( 0))
      {
        return linal::Vector< float>( feature_size, float( 0.0));
      }

      // get a reference to the features and results
      linal::MatrixConstReference< float> features_ref( DATASET.GetFeaturesReference());
      linal::MatrixConstReference< float> results_ref( DATASET.GetResultsReference());

      // standard deviations scaled by delta
      linal::Vector< float> scaled_sd;
      math::RunningAverageSD< linal::Vector< float> > result_std;
      {
        // setup running averages / standard deviations
        // initialize them all with a feature vector so that the running statistic will be for vectors of the right size
        math::RunningAverageSD< linal::Vector< float> > ave_std;

        // compute averages and variances
        for( size_t counter( 0); counter < dataset_size; ++counter)
        {
          // add the feature to the averages for the actives
          ave_std += features_ref.GetRow( counter);
          result_std += results_ref.GetRow( counter);
        }

        // compute standard deviation for each feature
        scaled_sd = ave_std.GetSampleStandardDeviation();
        scaled_sd *= m_Delta;
        m_ResultStd = result_std.GetSampleStandardDeviation();
        for( size_t counter( 0); counter < result_size; ++counter)
        {
          if( m_ResultStd( counter) < float( 0.001))
          {
            m_ResultStd( counter) = float( 0.001);
          }
        }
      }
      m_DescriptorStd = util::ToSiPtr( scaled_sd);
      m_DatasetPtr    = util::ToSiPtr( DATASET);
      m_FeatureNumber = 0;

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
      m_Predictions =
          storage::Vector< linal::Matrix< float> >
          (
            dataset_size,
            linal::Matrix< float>( result_size, m_Models->GetSize())
          );
      m_ChoseFeatures = false;

      util::ShPtrVector< sched::JobInterface> jobs( n_threads);
      storage::Vector< size_t> thread_ids( n_threads, size_t( 0));

      // initialize data for each thread
      for( size_t thread_number( 0); thread_number < n_threads; ++thread_number)
      {
        thread_ids( thread_number) = thread_number;
        jobs( thread_number) =
          util::ShPtr< sched::JobInterface>
          (
            new sched::UnaryFunctionJobWithData< const size_t, void, ScoreDatasetInputSensitivity>
            (
              0,
              *this,
              &ScoreDatasetInputSensitivity::RunThread,
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

      // return input-sensitivity for every feature column in dataset
      return m_AveDescriptorScores( 0).GetAverage();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool ScoreDatasetInputSensitivity::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // get the models if they are not already in memory
      RetrieveInterface::t_Container &models( s_Models[ GetString()]);

      if( m_Key.empty())
      {
        // get all the keys for this storage
        storage::Vector< std::string> keys( m_ModelStorage->GetAllKeys());

        BCL_MessageStd( "Found " + util::Format()( keys.GetSize()) + " models for input sensitivity");

        if( models.IsEmpty() && !keys.IsEmpty())
        {
          // no models already loaded
          RetrieveInterface::t_Container interfaces( m_ModelStorage->RetrieveEnsemble( keys));
          models = RetrieveInterface::t_Container( interfaces.Begin(), interfaces.End());
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
    io::Serializer ScoreDatasetInputSensitivity::GetSerializer() const
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
        "delta",
        "change of each feature used to calculate the derivative; will be scaled by the standard deviation of the feature",
        io::Serialization::GetAgentWithRange( &m_Delta, 0.0, 200.0),
        "0.5"
      );

      parameters.AddInitializer
      (
        "weights",
        "Change the weighting of various measures the derivative scores",
        io::Serialization::GetAgent( &m_Scorer)
      );

      parameters.AddInitializer
      (
        "feature limit",
        "limit the # of features that are used.   Maximum # of features to compute input sensitivity for, 0 to use them all. "
        "If the final objective function was classification-based, selects features with the highest value of: "
        " # models correct * # models incorrect / ( # models correct ^ 2 + # models incorrect ^ 2) "
        " For regression targets, selects features with the highest result RMSD standard deviation",
        io::Serialization::GetAgent( &m_MaxFeatures),
        "0"
      );

      return parameters;
    }

    //! @brief Get the next feature for a thread to make predictions for
    //! @param LAST_FEATURE last feature that this thread computed
    //! @return the next feature for a thread to make predictions for
    size_t ScoreDatasetInputSensitivity::GetNextFeatureForPrediction( const size_t &LAST_FEATURE) const
    {
      m_FeatureNumberMutex.Lock();
      const size_t dataset_size( m_DatasetPtr->GetSize());

      if( util::IsDefined( LAST_FEATURE))
      {
        BCL_Assert
        (
          !m_ChoseFeatures,
          "Should not have already chosen features while other threads are still predicting!"
        );
        m_FeaturesToCompute.PushBack( std::make_pair( 0.0, LAST_FEATURE));
        ++m_FeaturesToComputeSize;
      }
      else if( m_ChoseFeatures)
      {
        m_FeatureNumberMutex.Unlock();
        // thread never had to predict anything; all predictions already made by other threads
        return dataset_size;
      }

      if( m_FeatureNumber == dataset_size)
      {
        if( m_FeaturesToComputeSize == dataset_size)
        {
          ChooseInputSensitivityFeatures();
          m_FeatureNumber = 0;
          m_ChoseFeatures = true;
        }
        m_FeatureNumberMutex.Unlock();
        return dataset_size;
      }

      // only write out status from the first thread
      // determine progress percent
      const size_t percent( float( m_FeatureNumber) * 100.0f / dataset_size);

      // determine number of stars in the status bar
      const size_t number_stars( percent / 5);

      const std::string status
      (
        "["
        + std::string( number_stars, '*')
        + std::string( 20 - number_stars, ' ')
        + "] "
        + util::Format()( percent) + "% "
        + util::Format()( m_FeatureNumber) + " / " + util::Format()( size_t( dataset_size))
        + " features predicted with"
      );

      const size_t number_to_return( m_FeatureNumber);

      ++m_FeatureNumber;
      util::GetLogger().LogStatus( status);
      m_FeatureNumberMutex.Unlock();
      return number_to_return;
    }

    //! @brief Get the next feature for a thread to compute sensitivity for
    //! @return the next feature for a thread to use input sensitivity with
    size_t ScoreDatasetInputSensitivity::GetNextFeatureForInputSensitivity() const
    {
      m_FeatureNumberMutex.Lock();
      const size_t dataset_size( m_DatasetPtr->GetSize());

      // check for whether this thread is done
      if( !m_FeaturesToComputeSize)
      {
        m_FeatureNumberMutex.Unlock();
        return dataset_size;
      }

      const size_t n_to_choose( m_MaxFeatures ? std::min( m_MaxFeatures, dataset_size) : dataset_size);

      // only write out status from the first thread
      // determine progress percent
      const size_t percent( float( m_FeatureNumber) * 100.0f / n_to_choose);

      // determine number of stars in the status bar
      const size_t number_stars( percent / 5);

      const std::string status
      (
        "["
        + std::string( number_stars, '*')
        + std::string( 20 - number_stars, ' ')
        + "] "
        + util::Format()( percent) + "% "
        + util::Format()( m_FeatureNumber) + " / " + util::Format()( size_t( n_to_choose))
        + " features tested"
      );

      const size_t number_to_return( m_FeaturesToCompute.FirstElement().second);
      m_FeaturesToCompute.PopFront();
      --m_FeaturesToComputeSize;

      ++m_FeatureNumber;
      util::GetLogger().LogStatus( status);
      m_FeatureNumberMutex.Unlock();
      return number_to_return;
    }

    //! @brief choose the features on which to compute input sensitivity
    void ScoreDatasetInputSensitivity::ChooseInputSensitivityFeatures() const
    {
      const size_t result_size( m_DatasetPtr->GetResultSize());
      const size_t dataset_size( m_DatasetPtr->GetSize());
      const size_t n_models( m_Models->GetSize());

      linal::Matrix< float> model_predictions( dataset_size, result_size);
      m_PredictionClassifications = storage::Vector< linal::Matrix< char> >( n_models);
      FeatureDataReference< float> experimental_results_fds( m_DatasetPtr->GetResultsReference());
      if( m_DatasetPtr->GetIdsPtr().IsDefined())
      {
        m_Objective->SetData( experimental_results_fds, *m_DatasetPtr->GetIdsPtr());
      }
      else
      {
        m_Objective->SetData( experimental_results_fds);
      }
      FeatureDataReference< float> model_predictions_fds( model_predictions);

      // iterate over each model's predictions
      for( size_t model( 0); model < n_models; ++model)
      {
        // translate the predictions from m_Predictions into the predictions for this model
        for( size_t datum( 0); datum < dataset_size; ++datum)
        {
          for( size_t result( 0); result < result_size; ++result)
          {
            model_predictions( datum, result) = m_Predictions( datum)( result, model);
          }
        }
        // have the objective function compute the classifications
        m_PredictionClassifications( model)
          = m_Objective->GetFeaturePredictionClassifications( experimental_results_fds, model_predictions_fds);
      }
      m_Scorer.InitializeBalancing( m_PredictionClassifications, m_DatasetPtr->GetFeatureSize());

      if( !m_MaxFeatures || m_FeaturesToComputeSize <= m_MaxFeatures)
      {
        return;
      }

      // score each feature
      for
      (
        storage::List< std::pair< float, size_t> >::iterator
          itr( m_FeaturesToCompute.Begin()), itr_end( m_FeaturesToCompute.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Vector< std::string> model_classifications
        (
          ScoreDatasetNeuralNetworkInputSensitivity::PartitionModels( m_PredictionClassifications, itr->second)
        );

        // compute a score based off of the number of models that got this feature wrong
        math::RunningAverage< float> ave_result_reduced_mass;
        for( size_t result( 0); result < result_size; ++result)
        {
          size_t number_correct( 0);
          for( size_t model_number( 0); model_number < n_models; ++model_number)
          {
            if( isupper( m_PredictionClassifications( model_number)( itr->second, result)))
            {
              ++number_correct;
            }
          }
          // compute reduced mass
          ave_result_reduced_mass +=
            float( 4 * number_correct * ( n_models - number_correct))
            / float( math::Sqr( number_correct) + math::Sqr( n_models - number_correct));
        }
        itr->first = ave_result_reduced_mass.GetAverage();
      }

      // descending sort
      m_FeaturesToCompute.Sort( std::greater< std::pair< float, size_t> >());

      // find the given location
      storage::List< std::pair< float, size_t> >::iterator itr( m_FeaturesToCompute.Begin());
      std::advance( itr, m_MaxFeatures);

      // get the value
      float score( itr->first);

      // find the surrounding range with the same value
      storage::List< std::pair< float, size_t> >::iterator itr_start( itr), itr_last( itr);
      for
      (
        storage::List< std::pair< float, size_t> >::iterator itr_begin( m_FeaturesToCompute.Begin());
        itr_start != itr_begin && itr_start->first == score;
        --itr_start
      )
      {
      }
      if( itr_start->first != score)
      {
        ++itr_start;
      }
      for
      (
        storage::List< std::pair< float, size_t> >::iterator itr_end( m_FeaturesToCompute.End());
        itr_last != itr_end && itr_last->first == score;
        ++itr_last
      )
      {
      }
      storage::List< std::pair< float, size_t> >::iterator itr_last_prev( itr_last);
      --itr_last_prev;
      if( itr_last != m_FeaturesToCompute.End())
      {
        m_FeaturesToCompute.Remove( itr_last, m_FeaturesToCompute.End());
        itr_last = m_FeaturesToCompute.End();
      }
      if( itr != itr_last_prev)
      {
        // create a vector with the cases that had the same scores
        std::vector< std::pair< float, size_t> > same_scores( itr_start, itr_last);
        std::random_shuffle( same_scores.begin(), same_scores.end());
        // erase the list following itr_start
        m_FeaturesToCompute.Remove( itr_start, m_FeaturesToCompute.End());
        m_FeaturesToComputeSize = m_FeaturesToCompute.GetSize();
        size_t n_to_append( m_MaxFeatures - m_FeaturesToComputeSize);
        // append the elements from the vector that were selected
        m_FeaturesToCompute.Append
        (
          same_scores.begin(),
          same_scores.begin() + n_to_append
        );
      }
      m_FeaturesToComputeSize = m_MaxFeatures;
    }

    //! @brief run a thread to compute the kernel for all features with ID = THREAD_ID % n_threads
    //! @param THREAD_ID id of the thread to run (0-indexed)
    void ScoreDatasetInputSensitivity::RunThread( const size_t &THREAD_ID) const
    {
      const size_t feature_size( m_DatasetPtr->GetFeatureSize());
      const size_t result_size( m_DatasetPtr->GetResultSize());
      const size_t dataset_size( m_DatasetPtr->GetSize());

      const RetrieveInterface::t_Container &models( *m_Models);
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
        linal::Matrix< float> &model_predictions( m_Predictions( feature_number));
        // for each model
        for( size_t model_number( 0); model_number < n_models; ++model_number)
        {
          // get a reference to the model
          const Interface &model( *models( model_number));

          // compute the original results with this model
          const FeatureDataSet< float> original_results_fds( model( original_feature));
          const linal::VectorConstReference< float> original_results( original_results_fds.GetMatrix().GetRow( 0));
          for( size_t result_number( 0); result_number < result_size; ++result_number)
          {
            model_predictions( result_number, model_number) = original_results( result_number);
          }
        }
      }

      // barrier until all predictions were made
      while( !m_ChoseFeatures)
      {
        util::Time::Delay( util::Time( 0, 1));
      }

      // compute input sensitivities on the selected features
      for
      (
        size_t feature_number( GetNextFeatureForInputSensitivity());
        feature_number < dataset_size;
        feature_number = GetNextFeatureForInputSensitivity()
      )
      {
        // create a feature data set containing just that feature
        linal::Matrix< float> features_matrix( size_t( 1), feature_size, dataset_reference[ feature_number]);
        linal::Matrix< float> &model_predictions( m_Predictions( feature_number));

        // compute the original results with this model
        const linal::Matrix< float> original_results( model_predictions.Transposed());

        // for each model
        for( size_t model_number( 0); model_number < n_models; ++model_number)
        {
          linal::Matrix< float> &model_weight_effects( weight_effects( model_number));

          // get a reference to the model
          const Interface &model( *models( model_number));

          FeatureDataReference< float> features = linal::MatrixConstReference< float>( features_matrix);
          linal::VectorConstReference< float> original_results_row( original_results.GetRow( model_number));
          for( size_t descriptor_col( 0); descriptor_col < feature_size; ++descriptor_col)
          {
            const float descriptor_std( ( *m_DescriptorStd)( descriptor_col));
            features_matrix( 0, descriptor_col) += descriptor_std;

            // test the perturbed matrix
            FeatureDataSet< float> perturbed_results( model( features));
            float *itr_effect( model_weight_effects[ descriptor_col]);
            for
            (
              const float *itr( original_results_row.Begin()),
                          *itr_end( original_results_row.End()),
                          *itr_perturbed( perturbed_results.GetMatrix().Begin());
              itr != itr_end;
              ++itr, ++itr_perturbed, ++itr_effect
            )
            {
              *itr_effect = *itr_perturbed - *itr;
            }
            features_matrix( 0, descriptor_col) -= descriptor_std;
          }
        }

        m_Scorer.Score
        (
          weight_effects,
          ScoreDatasetNeuralNetworkInputSensitivity::PartitionModels( m_PredictionClassifications, feature_number),
          m_AveDescriptorScores( THREAD_ID)
        );
      }
    }

  } // namespace model
} // namespace bcl

