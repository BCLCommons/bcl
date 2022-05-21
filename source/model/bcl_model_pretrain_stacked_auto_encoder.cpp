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
#include "model/bcl_model_pretrain_stacked_auto_encoder.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_interface_store_in_file.h"
#include "model/bcl_model_neural_network_selective_backpropagation_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    const util::SiPtr< const util::ObjectInterface> PretrainStackedAutoEncoder::s_PretrainInstance
    (
      util::Enumerated< PretrainNeuralNetworkInterface>::AddInstance( new PretrainStackedAutoEncoder())
    );

    //! @brief default constructor
    PretrainStackedAutoEncoder::PretrainStackedAutoEncoder() :
      m_Approximator( true),
      m_TrainingData(),
      m_EncodedFeatures(),
      m_CurrentAutoEncoderModel(),
      m_AutoEncoderArchitecture()
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
    PretrainStackedAutoEncoder::PretrainStackedAutoEncoder
    (
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA
    ) :
      m_Approximator( true),
      m_TrainingData(),
      m_EncodedFeatures(),
      m_CurrentAutoEncoderModel(),
      m_AutoEncoderArchitecture()
    {
      // set and rescale training data set
      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new PretrainStackedAutoEncoder copied from this instance
    PretrainStackedAutoEncoder *PretrainStackedAutoEncoder::Clone() const
    {
      return new PretrainStackedAutoEncoder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PretrainStackedAutoEncoder::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PretrainStackedAutoEncoder::GetAlias() const
    {
      static const std::string s_Name( "AutoEncoder");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void PretrainStackedAutoEncoder::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      DATA->GetFeatures().DeScale();
      DATA->GetFeatures().Rescale( m_Approximator.GetTransferFunction()->GetDynamicOutputRange(), RescaleFeatureDataSet::e_MinMax);

      // setup dataset for autoencoding, features and results are the same
      m_TrainingData = util::ShPtr< descriptor::Dataset>
      (
        new descriptor::Dataset( DATA->GetFeaturesPtr(), DATA->GetFeaturesPtr())
      );

    } // SetTrainingData

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> PretrainStackedAutoEncoder::GetCurrentModel() const
    {
      return m_CurrentAutoEncoderModel;
    }

    //! @brief the main operation, pretrains a neural network
    //! @param DATA the data for use in pretraining
    //! @param OBJECTIVE ShPtr to the objective function for the network
    util::ShPtr< NeuralNetwork> PretrainStackedAutoEncoder::PretrainNetwork
    (
      util::ShPtr< descriptor::Dataset> &DATA,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
    )
    {
      m_AutoEncoderArchitecture = m_Approximator.GetHiddenArchitecture();

      SetTrainingData( DATA);
      while( !m_AutoEncoderArchitecture.IsEmpty())
      {
        TrainNextLayer( OBJECTIVE);
      }

      if( !m_AutoEncoderStoragePath.empty())
      {
        // serialize encoder
        InterfaceStoreInFile model_storage( m_AutoEncoderStoragePath);

        std::string result_label( "Combined(0");
        for( size_t cnt( 1), num_input( m_CurrentAutoEncoderModel->GetNumberOutputs()); cnt < num_input; ++cnt)
        {
          result_label += ",";
          result_label += util::Format()( cnt);
        }
        result_label += ")";

        model_storage.StoreResultDescriptor( util::ObjectDataLabel( result_label));

        model_storage.Store( util::ShPtr< Interface>( m_CurrentAutoEncoderModel), this->GetFeatureCode().GetLabel());
      }

      // train output for approximator
      BCL_MessageStd( "Train output layer for usage as pretrainer");

      // make it compatible with approximator that surrounds this pretrainer
      m_AutoEncoderArchitecture.PushBack( DATA->GetResultSize());
      while( !m_AutoEncoderArchitecture.IsEmpty())
      {
        TrainNextLayer( OBJECTIVE);
      }

      return GetCurrentModel();
    }

    //! @brief conducts the next approximation step and stores the approximation
    void PretrainStackedAutoEncoder::TrainNextLayer( const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE)
    {
      BCL_MessageStd( "Construct autoencoder layer ...");

      // set up AutoEncoder (AE) architecture
      storage::Vector< size_t> ae_architecture;
      ae_architecture.PushBack( m_AutoEncoderArchitecture.FirstElement());

      // remove current layer from considering future ae layers
      m_AutoEncoderArchitecture.RemoveElements( 0, 1);

      util::ShPtr< descriptor::Dataset> encoded_train_data;
      if( m_CurrentAutoEncoderModel.IsDefined())
      {

        encoded_train_data = util::ShPtr< descriptor::Dataset>
        (
          new descriptor::Dataset( m_EncodedFeatures, m_EncodedFeatures)
        );
      }
      else
      {
        BCL_MessageStd( "No encode of training data!");
      }

      // in round 1 take training data, for subsequent rounds take encoded data to train
      util::ShPtr< descriptor::Dataset> encoded_training_data
      (
        encoded_train_data.IsDefined() ? encoded_train_data : m_TrainingData
      );

      // set up NN to train current AE layer
      ApproximatorNeuralNetwork autoencoder( m_Approximator);
      // set AE architecture for current layer
      autoencoder.SetHiddenArchitecture( ae_architecture);

      // train current AE layer
      util::ShPtr< NeuralNetwork> current_ae_layer( autoencoder.PretrainNetwork( encoded_training_data, OBJECTIVE));

      // construct encoder portion
      util::ShPtr< NeuralNetwork> current_ae_encoder( new NeuralNetwork( *current_ae_layer));
      current_ae_encoder->RemoveOutputLayer();

      m_EncodedFeatures = util::ShPtr< FeatureDataSet< float> >
      (
        new FeatureDataSet< float>( current_ae_encoder->PredictWithoutRescaling( *encoded_training_data->GetFeaturesPtr()))
      );

      // construct all layers to resamble current encoder
      if( m_CurrentAutoEncoderModel.IsDefined())
      {
        // encoder
        util::ShPtr< NeuralNetwork> previous_model( new NeuralNetwork( *m_CurrentAutoEncoderModel));
        previous_model->Append( current_ae_encoder);
        current_ae_encoder = previous_model;
      }

      // construct decoder portion
      util::ShPtr< NeuralNetwork> current_ae_decoder( new NeuralNetwork( *current_ae_layer));
      current_ae_decoder->RemoveInputLayer();

      // construct all layers to resamble current encoder
      if( m_CurrentAutoDecoderModel.IsDefined())
      {
        // decoder
        current_ae_decoder->Append( m_CurrentAutoDecoderModel.HardCopy());
      }

      // update cuurent encoder and decoder
      m_CurrentAutoEncoderModel = current_ae_encoder;
      m_CurrentAutoDecoderModel = current_ae_decoder;

      BCL_MessageStd( "Current Encoder: " + util::Format()( current_ae_encoder->GetArchitecture()));
      BCL_MessageStd( "Current Decoder: " + util::Format()( current_ae_decoder->GetArchitecture()));

      // test the ability to auto encode of the current model

//      util::ShPtr< NeuralNetwork> ae_model( new NeuralNetwork( *m_CurrentAutoEncoderModel));
//      ae_model->Append( m_CurrentAutoDecoderModel);
//
//      BCL_Debug( ae_model->GetArchitecture());
//
//      util::ShPtr< FeatureDataSet< float> > test_ae
//      (
//        new FeatureDataSet< float>( ae_model->PredictWithoutRescaling( *( m_TrainingData->GetFeaturesPtr())))
//      );
//
//      BCL_MessageStd( util::Format()( m_TrainingData->GetFeaturesPtr()->GetRawMatrix().GetRow( 0)).substr(0,300));
//      BCL_MessageStd( util::Format()( test_ae->GetRawMatrix().GetRow( 0)).substr(0,300));

      BCL_MessageStd( "Construct autoencoder layer ... done");
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PretrainStackedAutoEncoder::GetSerializer() const
    {
      io::Serializer parameters( m_Approximator.GetSerializer());
      parameters.SetClassDescription
      (
        "trains a stacked auto encoder based on neural networks (see http://en.wikipedia.org/wiki/Artificial_neural_network)"
      );

      parameters.AddInitializer
      (
        "model storage path",
        "output path that specifies the directory name where the autoencoder model is stored",
        io::Serialization::GetAgent( &m_AutoEncoderStoragePath),
        ""
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
