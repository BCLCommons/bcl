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
#include "bcl_descriptor_kohonen_map_info.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_operations.h"
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

    //! @brief default constructor
    template< typename t_DataType>
    KohonenMapInfo< t_DataType>::KohonenMapInfo() :
      m_NumberOutputs( 0),
      m_Positions( false),
      m_Distances( true),
      m_NumberModels( 0)
    {
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    KohonenMapInfo< t_DataType> *KohonenMapInfo< t_DataType>::Clone() const
    {
      return new KohonenMapInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &KohonenMapInfo< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the type of this descriptor
    //! @return the type of this descriptor (should ignore dimension setting)
    template< typename t_DataType>
    Type KohonenMapInfo< t_DataType>::GetType() const
    {
      return m_Type;
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &KohonenMapInfo< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "KohonenMapInfo");
      return s_name;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief load the models: called the 1st time that recalculate is called
    template< typename t_DataType>
    void KohonenMapInfo< t_DataType>::LoadModels()
    {
      // get all the keys for this storage
      storage::Vector< std::string> keys( m_ModelStorage->GetAllKeys());
      s_ModelsMutex.Lock();
      // get the models if they are not already in memory
      storage::Vector< util::OwnPtr< model::KohonenNetworkAverage> > &models( s_Models[ m_ModelStorage.GetString()]);
      if( models.GetSize() == size_t( 0))
      {
        // no models already loaded
        if( keys.GetSize() != size_t( 0))
        {
          BCL_MessageStd( "Loading models");
          models.AllocateMemory( keys.GetSize());
          // retrieve the models
          for
          (
            storage::Vector< std::string>::const_iterator itr( keys.Begin()), itr_end( keys.End());
            itr != itr_end;
            ++itr
          )
          {
            // add the model associated with this key
            models.PushBack( m_ModelStorage->Retrieve( *itr));
            BCL_Assert( models.LastElement().IsDefined(), "Found model in storage that is not a kohonen network!");
          }
        }
      }
      BCL_Assert( models.GetSize() == keys.GetSize(), "Number of models changed in storage!");
      m_Models = util::ToSiPtr( models);
      s_ModelsMutex.Unlock();
      m_NumberModels = m_Models->GetSize();
      m_NumberOutputs = 0;
      for
      (
        storage::Vector< util::OwnPtr< model::KohonenNetworkAverage> >::const_iterator
          itr( m_Models->Begin()), itr_end( m_Models->End());
        itr != itr_end;
        ++itr
      )
      {
        m_NumberOutputs += m_Distances ? ( *itr)->GetCodeBook().GetSize() : 0;
        m_NumberOutputs += m_Positions ? ( *itr)->GetCodeBook().FirstElement().GetPosition().GetSize() : 0;
      }

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
    bool KohonenMapInfo< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      // get all the keys for this storage
      storage::Vector< std::string> keys( m_ModelStorage->GetAllKeys());
      m_NumberModels = keys.GetSize();

      BCL_MessageVrb
      (
        "Found " + util::Format()( keys.GetSize()) + " models for predictions"
      );

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
      if( !m_Distances && !m_Positions)
      {
        ERR_STREAM << "Must compute one or both of distances or positions";
        return false;
      }
      // # of outputs = # of outputs per model * # of models
      m_Type = result->GetType();
      LoadModels();

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer KohonenMapInfo< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes distances to all members of one or more kohonen maps. May alternatively/additionally compute position "
        " of an element on the map. "
        "If both are computed, output will be in the format: position in kohonen map 1, "
        "distances to nodes of kohonen map 1, position in kohonen map 2, etc."
      );

      parameters.AddInitializer
      (
        "distances",
        "Whether to compute distances to all nodes of the kohonen map",
        io::Serialization::GetAgent( &m_Distances),
        "True"
      );
      parameters.AddInitializer
      (
        "positions",
        "Whether to return positions of best matching node in each kohonen map",
        io::Serialization::GetAgent( &m_Positions),
        "False"
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
    iterate::Generic< Base< t_DataType, float> > KohonenMapInfo< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( m_Properties.Begin(), m_Properties.End());
    }

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void KohonenMapInfo< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // load the models, if they are not already loaded
      if( !m_Models.IsDefined())
      {
        LoadModels();
      }
      if( m_Models->IsEmpty())
      {
        return;
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
        const linal::VectorConstInterface< float> &descriptors( ( *itr)( ITR));

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
      const size_t number_models( m_Models->GetSize());
      size_t storage_pos( 0);
      for( size_t model_number( 0); model_number < number_models; ++model_number)
      {
        // get a reference to this model
        const model::KohonenNetworkAverage &model( *( *m_Models)( model_number));
        const storage::Vector< model::KohonenNode> &nodes( model.GetCodeBook());

        // get the feature for this model
        model::FeatureDataSet< float> &features_this_model( descriptors_vect( m_DescriptorIDs( model_number)));
        model.Rescale( features_this_model);

        if( !m_Distances)
        {
          // no need to compute all distances. Use the optimized GetIndexOfWinningNode
          const size_t indx( model.GetIndexOfWinningNode( features_this_model( 0)));
          const linal::VectorConstInterface< float> &position( nodes( indx).GetPosition());
          STORAGE.CreateSubVectorReference( position.GetSize(), storage_pos).CopyValues( position);
          storage_pos += position.GetSize();
        }
        else
        {
          double best_distance( std::numeric_limits< double>::max());
          size_t best_distance_index( 0);
          const linal::VectorConstReference< float> descriptor( features_this_model( 0));
          for( size_t node_id( 0), n_nodes( nodes.GetSize()); node_id < n_nodes; ++node_id)
          {
            const double dist( linal::Distance( descriptor, nodes( node_id).GetFeatureVector()));
            STORAGE( storage_pos++) = dist;
            if( dist < best_distance)
            {
              best_distance = dist;
              best_distance_index = node_id;
            }
          }
          if( m_Positions)
          {
            const linal::VectorConstInterface< float> &position( nodes( best_distance_index).GetPosition());
            STORAGE.CreateSubVectorReference( position.GetSize(), storage_pos).CopyValues( position);
            storage_pos += position.GetSize();
          }
        }
        features_this_model.DeScale();
      }
    }

  } // namespace descriptor
} // namespace bcl
