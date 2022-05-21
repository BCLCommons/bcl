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
#include "bcl_descriptor_prediction.h"
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
    Prediction< t_DataType>::Prediction( const bool &MEAN) :
      m_NumberOutputs( 0),
      m_ComputeMean( MEAN),
      m_NumberModels( 0)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    Prediction< t_DataType> *Prediction< t_DataType>::Clone() const
    {
      return new Prediction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &Prediction< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the type of this descriptor
    //! @return the type of this descriptor (should ignore dimension setting)
    template< typename t_DataType>
    Type Prediction< t_DataType>::GetType() const
    {
      return m_Type;
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &Prediction< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "Prediction"), s_mean_name( "PredictionMean");
      return m_ComputeMean ? s_mean_name : s_name;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief load the models: called the 1st time that recalculate is called
    template< typename t_DataType>
    void Prediction< t_DataType>::LoadModels()
    {
      // get all the keys for this storage
      m_Models = m_ModelStorage->RetrieveEnsemble();

      // load all the properties
      storage::List< util::ObjectDataLabel> labels( m_ModelStorage->RetrieveEnsembleDescriptors());

      m_DescriptorIDs.AllocateMemory( labels.GetSize());
      m_Properties.AllocateMemory( labels.GetSize());
      for
      (
        storage::List< util::ObjectDataLabel>::const_iterator itr( labels.Begin()), itr_end( labels.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::List< util::ObjectDataLabel>::const_iterator itr_descriptor( labels.Begin());
        size_t descriptor_id_cnt( 0);
        for( ; itr_descriptor != itr; ++itr_descriptor, ++descriptor_id_cnt)
        {
          if( *itr_descriptor == *itr)
          {
            break;
          }
        }
        m_DescriptorIDs.PushBack( descriptor_id_cnt);
        // check whether the descriptor set was used before
        if( itr_descriptor == itr)
        {
          s_DescriptorsMutex.Lock();
          // descriptor set has not been seen before
          // create a small molecule properties combined to hold all the descriptors for this model
          Combine< t_DataType, float> descriptors;
          std::stringstream err_stream;
          if( !descriptors.TryRead( *itr, err_stream))
          {
            BCL_Exit
            (
              "Failed to read descriptor for model for " + this->GetString() + "; error was: " + err_stream.str(),
              -1
            );
          }
          // add the descriptor combine to the properties vector
          m_Properties.PushBack( descriptors);
          s_DescriptorsMutex.Unlock();
          // set the type up properly
          m_Properties.LastElement().SetDimension( m_Type.GetDimension());
        }
      }
      BCL_MessageVrb
      (
        "Found " + util::Format()( m_Properties.GetSize()) + " descriptor sets to calculate for predictions"
      );
      // if set object was already called on this prediction object, then we must also call set object on the internal
      // descriptors
      if( this->GetCurrentObject().IsDefined())
      {
        for
        (
          iterate::Generic< Base< t_DataType, float> > itr_properties( m_Properties.Begin(), m_Properties.End());
          itr_properties.NotAtEnd();
          ++itr_properties
        )
        {
          itr_properties->SetObject( *this->GetCurrentObject());
        }
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool Prediction< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // get all the keys for this storage
      storage::Vector< std::string> keys( m_ModelStorage->GetAllKeys());
      m_NumberModels = keys.GetSize();

      BCL_MessageVrb( "Found " + util::Format()( m_NumberModels) + " models for predictions");

      // load the result label
      util::ObjectDataLabel result_descriptor( m_ModelStorage->RetrieveResultDescriptor());

      // determine the # of outputs it would  yield
      if( result_descriptor.IsEmpty())
      {
        ERR_STREAM << "Result descriptor was not given for prediction of " << LABEL.ToString();
        return false;
      }
      util::Implementation< Base< t_DataType, float> > result( result_descriptor, ERR_STREAM);

      if( !result.IsDefined())
      {
        ERR_STREAM << "Result descriptor was invalid in prediction directory";
        return false;
      }

      // # of outputs = # of outputs per model * # of models
      m_NumberOutputs = m_ComputeMean ? result->GetSizeOfFeatures() : result->GetSizeOfFeatures() * m_NumberModels;
      m_Type = result->GetType();

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer Prediction< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_ComputeMean
        ? "computes the mean prediction of pre-trained machine learning model(s) on the given object"
        : "uses pre-trained machine learning model(s) on the given object"
      );

      parameters.AddInitializer
      (
        "storage",
        "type of storage for models",
        io::Serialization::GetAgent( &m_ModelStorage)
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
    iterate::Generic< Base< t_DataType, float> > Prediction< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( m_Properties.Begin(), m_Properties.End());
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void Prediction< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // load the models, if they are not already loaded
      if( m_Models.IsEmpty())
      {
        LoadModels();
      }

      // create all the unique descriptor sets
      storage::Vector< model::FeatureDataSet< float> > descriptors_vect;
      descriptors_vect.AllocateMemory( m_Properties.GetSize());

      for
      (
        typename storage::Vector< Combine< t_DataType, float> >::iterator
          itr( m_Properties.Begin()), itr_end( m_Properties.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine the # of features to use for this descriptor
        const size_t feature_size( itr->GetSizeOfFeatures());

        // calculate this unique descriptor set
        // calculate the properties
        const linal::Vector< float> descriptors( ( *itr)( ITR));

        if( feature_size != descriptors.GetSize())
        {
          // missing features; push back an undefined dataset
          descriptors_vect.PushBack( model::FeatureDataSet< float>( 1, feature_size, util::GetUndefined< float>()));
        }
        else
        {
          // convert into a feature dataset
          const linal::Matrix< float> descriptors_mat( size_t( 1), feature_size, descriptors.Begin());

          // get the descriptors for this model, and store them in a feature data set
          descriptors_vect.PushBack( model::FeatureDataSet< float>( descriptors_mat));
        }
      }

      // iterate over the cached models
      const size_t number_models( m_Models.GetSize());
      if( m_ComputeMean)
      {
        math::RunningAverage< linal::Vector< float> > ave_prediction;
        for( size_t model_number( 0); model_number < number_models; ++model_number)
        {
          // get a reference to this model
          const model::Interface &model( *m_Models( model_number));

          // get the feature for this model
          const model::FeatureDataSet< float> &features_this_model( descriptors_vect( m_DescriptorIDs( model_number)));

          // make the predictions with the model
          model::FeatureDataSet< float> prediction_fds( model( features_this_model));

          // get the feature reference
          ave_prediction += prediction_fds( 0);
        }
        STORAGE.CopyValues( ave_prediction.GetAverage());
      }
      else
      {
        size_t storage_pos( 0), features_per_model( std::max( size_t( 1), m_NumberOutputs) / number_models);
        for( size_t model_number( 0); model_number < number_models; ++model_number, storage_pos += features_per_model)
        {
          // get a reference to this model
          const model::Interface &model( *m_Models( model_number));

          // get the feature for this model
          const model::FeatureDataSet< float> &features_this_model( descriptors_vect( m_DescriptorIDs( model_number)));

          // make the predictions with the model
          model::FeatureDataSet< float> prediction_fds( model( features_this_model));

          // get the feature reference
          STORAGE.CreateSubVectorReference( features_per_model, storage_pos).CopyValues( prediction_fds( 0));
        }
      }
    }

  } // namespace descriptor
} // namespace bcl
