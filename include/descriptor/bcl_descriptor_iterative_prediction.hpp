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

// include header of this class
#include "bcl_descriptor_iterative_prediction.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_feature_data_set.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, accepts bool of whether to auto-compute mean
    template< typename t_DataType>
    IterativePrediction< t_DataType>::IterativePrediction( const bool &MEAN) :
      Prediction< t_DataType>( MEAN),
      m_TertiaryPredictor( false),
      m_Iterations( 2),
      m_TertiaryPredictionLabel( "Prediction"),
      m_TertiaryPredictionMeanLabel( "PredictionMean"),
      m_LowLevelPredictionLabel( m_TertiaryPredictionLabel),
      m_LowLevelPredictionMeanLabel( m_TertiaryPredictionMeanLabel)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    IterativePrediction< t_DataType> *IterativePrediction< t_DataType>::Clone() const
    {
      return new IterativePrediction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &IterativePrediction< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &IterativePrediction< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "IterativePrediction"), s_mean_name( "IterativePredictionMean");
      return this->DoesPerformMean() ? s_mean_name : s_name;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType>
    void IterativePrediction< t_DataType>::SetObjectHook()
    {
      m_LastRoundsPredictions = linal::Matrix< float>();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool IterativePrediction< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( !Prediction< t_DataType>::ReadInitializerSuccessHook( LABEL, ERR_STREAM))
      {
        return false;
      }

      // ensure that the number of models is the same in both storages
      if( this->GetNumberOfModels() != m_InitialModelsStorage->GetSize())
      {
        ERR_STREAM << "Iterative prediction only works if both the predictor and the internal prediction storage have "
                   << "the same number of models!";
        return false;
      }

      // ensure that the number of outputs is also the same
      util::ObjectDataLabel result_label_replaced( m_InitialModelsStorage->RetrieveResultDescriptor());
      util::Implementation< Base< t_DataType, float> > result_descriptor( result_label_replaced, ERR_STREAM);
      if( !result_descriptor.IsDefined())
      {
        ERR_STREAM << "Could not create result descriptor from " << m_InitialModelsStorage->GetString();
        return false;
      }
      if( this->GetModelsNumberOutputs() != result_descriptor->GetSizeOfFeatures())
      {
        ERR_STREAM << "Iterative prediction requires the primary predictor and the internal predictor to have the "
                   << "same size!";
        return false;
      }

      // setup the labels for cache retrieval
      {
        // low-level label
        util::ObjectDataLabel lowlevel_storage( "storage", m_InitialModelsStorage.GetLabel());
        storage::Vector< util::ObjectDataLabel> lowlevel_storage_vec( 1, lowlevel_storage);
        m_LowLevelPredictionLabel = util::ObjectDataLabel( "", "Prediction", lowlevel_storage_vec);
        m_LowLevelPredictionMeanLabel = util::ObjectDataLabel( "", "PredictionMean", lowlevel_storage_vec);
      }

      {
        // tertiary label
        util::ObjectDataLabel tertiary_storage( "storage", *this->GetLabel().FindName( "storage", false));
        storage::Vector< util::ObjectDataLabel> tertiary_storage_vec( 1, tertiary_storage);
        m_TertiaryPredictionLabel = util::ObjectDataLabel( "", "Prediction", tertiary_storage_vec);
        m_TertiaryPredictionMeanLabel = util::ObjectDataLabel( "", "PredictionMean", tertiary_storage_vec);
        // the tertiary predictor object is what will actually perform all the calculations
        if( !m_TertiaryPredictor.TryRead( tertiary_storage, ERR_STREAM))
        {
          return false;
        }
      }

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer IterativePrediction< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters( Prediction< t_DataType>::GetSerializer());
      parameters.SetClassDescription
      (
        this->DoesPerformMean()
        ? "computes the mean prediction of pre-trained machine learning model(s) on the given object, "
          "iteratively substituting the values returned by the sub-model with the values returned by this model"
        : "uses pre-trained machine learning model(s) on the given object"
      );

      parameters.AddInitializer
      (
        "lower layer",
        "storage for the predictor that this descriptor's predictions should overwrite.  Note that using this descriptor "
        "causes the cache to be modified; such that the original predictions of this model cannot be returned",
        io::Serialization::GetAgent( &m_InitialModelsStorage)
      );
      parameters.AddInitializer
      (
        "iterations",
        "number of iterations to perform",
        io::Serialization::GetAgent( &m_Iterations)
      );

      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > IterativePrediction< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_TertiaryPredictor, &m_TertiaryPredictor + 1);
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void IterativePrediction< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // load the models, if they are not already loaded
      if( m_LastRoundsPredictions.GetNumberRows() == size_t( 0))
      {
        BCL_MessageStd( "Recomputing iterative prediction");
        Iterate( ITR);
      }

      STORAGE.CopyValues( m_LastRoundsPredictions.GetRow( ITR.GetPosition()));
    }

    namespace
    {
      //! @brief compute prediction mean given the complete predictions vector
      //! @param STORAGE vector containing one spot for each predicted output for each predicted element
      //! @param PREDICTIONS the actual predictions
      //! @param N_FEATURES number of features represented in the vector
      void ComputePredictionMeanFromPredictions
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::VectorConstInterface< float> &PREDICTIONS,
        const size_t &N_FEATURES
      )
      {
        linal::VectorConstReference< float> predictions_ref( PREDICTIONS);
        float *storage_ptr( STORAGE.Begin());
        // PREDICTIONS is of size N_OUTPUTS * N_FEATURES * N_MODELS
        // STORAGE is of size N_OUTPUTS * N_FEATURES
        const size_t n_models( PREDICTIONS.GetSize() / STORAGE.GetSize());
        const size_t outputs_per_model( STORAGE.GetSize() / N_FEATURES);
        math::RunningAverage< linal::Vector< float> > running_average;
        size_t pos( 0);
        for( size_t feature_number( 0); feature_number < N_FEATURES; ++feature_number)
        {
          for( size_t model( 0); model < n_models; pos += outputs_per_model, ++model)
          {
            running_average += predictions_ref.CreateSubVectorConstReference( outputs_per_model, pos);
          }
          std::copy( running_average.GetAverage().Begin(), running_average.GetAverage().End(), storage_ptr);
          storage_ptr += outputs_per_model;
          running_average.Reset();
        }
      }
    }

    //! @brief Iteratively perform the prediction
    //! @param ITR iterator to use
    template< typename t_DataType>
    void IterativePrediction< t_DataType>::Iterate( const Iterator< t_DataType> &ITR)
    {
      // first, save all cache entries that use the lower-level prediction
      CacheMap lowlevel_cache_entries( this->GetCurrentObject()->ExtractRelatedCacheEntries( m_LowLevelPredictionLabel));
      lowlevel_cache_entries.Merge( this->GetCurrentObject()->ExtractRelatedCacheEntries( m_LowLevelPredictionMeanLabel));
      this->GetCurrentObject()->MergeCache( lowlevel_cache_entries);

      // first, compute the predictions for the upper layer using the lower layer's inputs directly
      // Note that because Predictions are always cached, it is sufficient to call the predictor's method with any
      // iterator for this object
      m_TertiaryPredictor( ITR);

      // Save the predictions from the lower layer; these will be re-inserted at the end
      const bool lowlevel_predictor_uses_all( this->GetCurrentObject()->IsCached( m_LowLevelPredictionLabel));
      const bool lowlevel_predictor_uses_mean( this->GetCurrentObject()->IsCached( m_LowLevelPredictionMeanLabel));

      BCL_Assert
      (
        lowlevel_predictor_uses_all || lowlevel_predictor_uses_mean,
        "Tertiary model did not use the given lower-level model storage"
      );

      const size_t n_outputs( Prediction< t_DataType>::GetModelsNumberOutputs());

      // save the predictions from the tertiary layer; these will change every time through
      util::SiPtr< const linal::Vector< float> > tertiary_all_predictions_ptr
      (
        this->GetCurrentObject()->FindInCache( m_TertiaryPredictionLabel)
      );

      BCL_Assert
      (
        tertiary_all_predictions_ptr.IsDefined(),
        "Tertiary predictor was not cached!!! Cache = "
        + util::Format()( *( this->GetCurrentObject()->GetCacheMap()))
      );
      linal::Vector< float> prev_tertiary_predictions( *tertiary_all_predictions_ptr);
      const size_t n_features( ITR.GetSize());
      linal::Vector< float> prev_tertiary_prediction_means( n_outputs * n_features);
      const size_t n_models( Prediction< t_DataType>::GetNumberOfModels());
      ComputePredictionMeanFromPredictions( prev_tertiary_prediction_means, prev_tertiary_predictions, n_features);
      const linal::Vector< float> original_tertiary_predictions( prev_tertiary_predictions);

      for( size_t iteration( 0); iteration < m_Iterations; ++iteration)
      {
        if( lowlevel_predictor_uses_all)
        {
          this->GetCurrentObject()->RemoveFromCache( m_LowLevelPredictionLabel, true);
          this->GetCurrentObject()->Cache( m_LowLevelPredictionLabel, prev_tertiary_predictions);
        }
        if( lowlevel_predictor_uses_mean)
        {
          // remove all descriptors computed using the low level predictors from the cache
          this->GetCurrentObject()->RemoveFromCache( m_LowLevelPredictionMeanLabel, true);
          this->GetCurrentObject()->Cache( m_LowLevelPredictionMeanLabel, prev_tertiary_prediction_means);
        }

        // remove the tertiary prediction descriptors associated with the high-level predictor from the cache
        // to force it to be recomputed using the new low-level prediction values
        this->GetCurrentObject()->RemoveFromCache( m_TertiaryPredictionLabel, false);

        // call set-object to force recomputation of whether each of the descriptors is already in the cache
        m_TertiaryPredictor.SetObject( *this->GetCurrentObject());

        // call the tertiary predictor with the iterator to cause this descriptor to be recached
        m_TertiaryPredictor( ITR);

        // save the predictions from the tertiary layer; these will change every time through
        tertiary_all_predictions_ptr = this->GetCurrentObject()->FindInCache( m_TertiaryPredictionLabel);
        BCL_Assert
        (
          tertiary_all_predictions_ptr.IsDefined(),
          "Tertiary predictor was not cached!!! Cache contained: "
          + util::Format()( *( this->GetCurrentObject()->GetCacheMap()))
        );
        prev_tertiary_predictions = *tertiary_all_predictions_ptr;
        ComputePredictionMeanFromPredictions( prev_tertiary_prediction_means, prev_tertiary_predictions, n_features);
      }

      // restore the low-level predictions within the cache to their proper values
      this->GetCurrentObject()->RemoveFromCache( m_LowLevelPredictionMeanLabel, true);
      this->GetCurrentObject()->RemoveFromCache( m_LowLevelPredictionLabel, true);
      this->GetCurrentObject()->MergeCache( lowlevel_cache_entries);

      // restore the tertiary predictions to the cache with their non-iterative values
      this->GetCurrentObject()->Cache( m_TertiaryPredictionLabel, original_tertiary_predictions);
      if( !this->GetCurrentObject()->IsCached( m_TertiaryPredictionMeanLabel))
      {
        linal::Vector< float> original_tertiary_prediction_means( prev_tertiary_prediction_means);
        ComputePredictionMeanFromPredictions( original_tertiary_prediction_means, original_tertiary_predictions, n_features);
        this->GetCurrentObject()->Cache( m_TertiaryPredictionMeanLabel, original_tertiary_prediction_means);
      }

      // set the tertiary predictor again to force it to recompute which descriptors are already in the cache
      m_TertiaryPredictor.SetObject( *this->GetCurrentObject());

      // store the last round's prediction values for use in the main recalculate function
      if( this->DoesPerformMean())
      {
        m_LastRoundsPredictions = linal::Matrix< float>( n_features, n_outputs, prev_tertiary_prediction_means.Begin());
      }
      else
      {
        m_LastRoundsPredictions = linal::Matrix< float>( n_features, n_outputs * n_models, prev_tertiary_predictions.Begin());
      }
    }

  } // namespace descriptor
} // namespace bcl
