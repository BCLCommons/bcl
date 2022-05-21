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
#include "bcl_app_model_test_ann_with_dropout_resampling.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_neural_network.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "model/bcl_model_retrieve_interface.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> ModelTestANNWithDropoutResampling::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "TestANNWithDropoutResampling");
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ModelTestANNWithDropoutResampling::GetReadMe() const
    {
      util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, float> >::SetHaveDisplayedHelp();
      static io::FixedLineWidthWriter writer;
      static std::string s_read_me =
      "Application for testing an ANN using dropout at test-time to compute the distribution of values for each output.";
      return s_read_me;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ModelTestANNWithDropoutResampling::GetDescription() const
    {
      return "Test an ANN using dropout at test-time to compute the distribution of values for each output.";
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ModelTestANNWithDropoutResampling::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_FlagOutputBase);
      // flags
      sp_cmd->AddFlag( m_FlagDataSetRetriever);

      // flag for model storage
      sp_cmd->AddFlag( m_FlagModelInterfaceStorage);
      sp_cmd->AddFlag( m_FlagDropoutRates);
      sp_cmd->AddFlag( m_FlagNumberSamples);
      // flag for model storage
      sp_cmd->AddFlag( m_FlagIDCode);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

  ////////////////////
  // helper methods //
  ////////////////////

    //! @brief write out model predictions and experimental data
    //! @param DATA dataset with experimental data and basis of predictions
    //! @param PREDICTIONS predicted data
    //! @param FILENAME filename for written out predictions
    void ModelTestANNWithDropoutResampling::WritePredictions
    (
      const descriptor::Dataset &DATA,
      const linal::MatrixConstInterface< float> &PREDICTIONS,
      const linal::MatrixConstInterface< float> &STDS,
      const std::string &FILENAME
    ) const
    {
      // get predictions
      const model::FeatureDataSetInterface< float> &experimental( *( DATA.GetResultsPtr()));

      // write experimental and predicted data
      io::OFStream output;
      io::File::MustOpenOFStream( output, FILENAME);

      // output experimental, predicted, id
      output << "#Experimental\tStd\tPredicted\t";
      const linal::MatrixConstReference< char> ids( DATA.GetIdsReference());
      const size_t feature_size( experimental.GetFeatureSize());
      const size_t id_size( ids.GetNumberCols());
      if( id_size)
      {
        output << "\"" << m_FlagIDCode->GetFirstParameter()->GetValue() << "\"\n";
      }
      else
      {
        output << '\n';
      }
      for( size_t data_id( 0), data_set_size( experimental.GetNumberFeatures()); data_id < data_set_size; ++data_id)
      {
        for( size_t row_element_id( 0); row_element_id < feature_size; ++row_element_id)
        {
          // get the actual result
          output << experimental[ data_id][ row_element_id] << '\t';
        }
        for( size_t row_element_id( 0); row_element_id < feature_size; ++row_element_id)
        {
          // get the actual result
          output << STDS( data_id, row_element_id) << '\t';
        }
        for( size_t row_element_id( 0); row_element_id < feature_size; ++row_element_id)
        {
          output << PREDICTIONS( data_id, row_element_id) << '\t';
        }
        if( id_size)
        {
          output << std::string( ids.GetRow( data_id).Begin(), id_size) << '\n';
        }
        else
        {
          output << '\n';
        }
      }
      io::File::CloseClearFStream( output);
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ModelTestANNWithDropoutResampling::Main() const
    {
      // stopwatch
      util::Stopwatch timer;
      timer.Start();

      // storage for the final model::Interface
      util::Implementation< model::RetrieveInterface> model_storage
      (
        m_FlagModelInterfaceStorage->GetFirstParameter()->GetValue()
      );
      BCL_Assert
      (
        model_storage.IsDefined(),
        "Model storage was undefined.  Use -storage_model flag"
      );

      BCL_MessageStd( "Model storage initialized ... ");

      // store models from storage
      storage::Vector< std::string> model_ids;
      storage::List< util::ObjectDataLabel> descriptor_model_list;

      BCL_MessageStd( "Start reading descriptors and model ids ...");

      // store models from storage
      model_ids = model_storage->GetAllKeys();
      descriptor_model_list = model_storage->RetrieveEnsembleDescriptors();

      BCL_MessageStd( "Done reading descriptors and model ids");

      BCL_MessageStd
      (
        "Model storage initialized ... done: " + util::Format()( timer.GetTotalTime().GetSeconds()) + " [sec]"
      );

      // reset timer
      timer.Reset();
      timer.Start();

      storage::Vector< std::string>::const_iterator itr_model_ids( model_ids.Begin());
      storage::List< util::ObjectDataLabel>::const_iterator itr_descriptors( descriptor_model_list.Begin());

      // create
      util::Implementation< model::RetrieveDataSetBase> data_set_retriever
      (
        m_FlagDataSetRetriever->GetFirstParameter()->GetValue()
      );
      data_set_retriever->SelectResults( model_storage->RetrieveResultDescriptor());
      data_set_retriever->SelectIdsGivenFilenameFlag( *m_FlagIDCode);
      util::ObjectDataLabel descriptors;
      util::ShPtr< descriptor::Dataset> data_set;

      // feature dataset that contains
      math::RunningAverage< linal::Matrix< float> > ave_preds;
      math::RunningAverage< linal::Matrix< float> > ave_sds;
      math::RunningAverage< linal::Matrix< float> > ave_preds_real;

      model::CrossValidationInfo common_info( model_storage->RetrieveCommonCVInfo());
      common_info.SetIndependentDataset( data_set_retriever->GetLabel());

      // iterate over all models
      for
      (
        storage::Vector< std::string>::const_iterator
          itr_model_ids( model_ids.Begin()), itr_model_ids_end( model_ids.End());
        itr_model_ids != itr_model_ids_end;
        ++itr_model_ids, ++itr_descriptors
      )
      {
        // current model
        model::RetrieveInterface::t_ModelPtr model( model_storage->Retrieve( *itr_model_ids));
        BCL_Assert( model.IsDefined(), "Read in model from storage is undefined!");
        util::SiPtr< const model::NeuralNetwork> ann_cast( model.GetPointer());
        BCL_Assert( ann_cast.IsDefined(), "Read in ann from storage is undefined!");

        // if descriptor set is different as previous descriptor set then calculate independent data again
        if( descriptors.ToString() != itr_descriptors->ToString())
        {
          // set descriptor set
          descriptors = *itr_descriptors;

          // generated independent dataset
          data_set_retriever->SelectFeatures( descriptors);
          data_set = data_set_retriever->GenerateDataSet();
          BCL_MessageStd( "Number of retrieved feature/result pairs: " + util::Format()( data_set->GetSize()));
        }

        // generate predicted data
        // scale input data to avoid unnecessary regeneration of the dataset
        util::ShPtr< model::FeatureDataSet< float> > features_ptr( data_set->GetFeaturesPtr());

        auto predictions
        (
          ann_cast->TestWithDropout
          (
            *features_ptr,
            m_FlagDropoutRates->GetNumericalList< double>(),
            m_FlagNumberSamples->GetFirstParameter()->GetNumericalValue< size_t>()
          )
        );
        ave_preds += predictions.First().GetMatrix();
        ave_sds += predictions.Second().GetMatrix();
        ave_preds_real += predictions.Third().GetMatrix();
      }

      if( m_FlagOutputBase->GetFlag())
      {
        // prepare output
        std::string filename( m_FlagOutputBase->GetParameterList().FirstElement()->GetValue());

        // write predictions for individual model
        common_info.SetIndependentPredictionsFilename( filename);
        common_info.WritePredictions
        (
          data_set->GetResultsReference(),
          ave_preds_real.GetAverage(),
          data_set->GetIdsReference(),
          filename
        );
        common_info.WritePredictions
        (
          data_set->GetResultsReference(),
          ave_sds.GetAverage(),
          data_set->GetIdsReference(),
          filename + ".stds.txt"
        );
        common_info.WritePredictions
        (
          data_set->GetResultsReference(),
          ave_preds.GetAverage(),
          data_set->GetIdsReference(),
          filename + ".resampled.txt"
        );
        io::OFStream output;
        io::File::MustOpenOFStream( output, filename + ".info");
        output << common_info;
        io::File::CloseClearFStream( output);
      }

      // end
      return 0;
    } // Main

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief standard constructor
    ModelTestANNWithDropoutResampling::ModelTestANNWithDropoutResampling() :
      m_FlagOutputBase
      (
        new command::FlagStatic
        (
          "output",
          "flag for output base name (if -average is also given, this will be the actual filename instead)",
          command::Parameter( "filename", "flag for output of screening", "prediction.txt")
        )
      ),
      m_FlagDataSetRetriever
      (
        new command::FlagStatic
        (
          "retrieve_dataset",
          "method to retrieve the dataset to test the model(s) against",
          command::Parameter
          (
            "retrieve_dataset",
            "method to retrieve the dataset to test the model(s) against",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_FlagModelInterfaceStorage
      (
        new command::FlagStatic
        (
          "storage_model",
          "location of models to predict with",
          command::Parameter
          (
            "model_storage",
            "choice of model storage",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveInterface>()),
            "File"
          )
        )
      ),
      m_FlagDropoutRates
      (
        new command::FlagDynamic
        (
          "dropout",
          "rates of dropout, one for each layer in the model (except output layer)",
          command::Parameter
          (
            "rate",
            "rate of dropout for the layer",
            command::ParameterCheckRanged< double>( 0.0, 1.0),
            "0.0"
          ),
          1,
          100000
        )
      ),
      m_FlagNumberSamples
      (
        new command::FlagStatic
        (
          "number_samples",
          "Number of resamplings to perform",
          command::Parameter
          (
            "number",
            "number of times to resample the network with a different dropout mask",
            command::ParameterCheckRanged< size_t>( 2.0, 10000.0),
            "100"
          )
        )
      ),
      m_FlagIDCode
      (
        new command::FlagStatic
        (
          "id_labels",
          "label or file containing the id label; necessary if ids for each row are desired",
          command::Parameter( "id_labels", "", "")
        )
      )
    {
    }

    const ApplicationType ModelTestANNWithDropoutResampling::ModelTestANNWithDropoutResampling_Instance
    (
      GetAppGroups().AddAppToGroup( new ModelTestANNWithDropoutResampling(), GetAppGroups().e_Model)
    );

  } // namespace app
} // namespace bcl
