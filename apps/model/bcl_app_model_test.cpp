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
#include "bcl_app_model_test.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
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
    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ModelTest::GetDescription() const
    {
      return "Test any machine learning model, including ANNs, SVMs, and many more";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ModelTest::GetReadMe() const
    {
      static const std::string readme =
          "ModelTest is an application for testing any machine learning model generated with the BCL. "
          "Usually, models built with the launch.py script or with model:Train are tested/applied using "
          "the descriptor framework, e.g. bcl.exe molecule:Properties -add 'PredictionMean(help)'; however, "
          "model:Test is an alternative way to access this functionality. The decision regarding which "
          "workflow to use is typically dependent on the task.";
      return readme;
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> ModelTest::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "TestModel");
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ModelTest::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_FlagOutputBase);
      sp_cmd->AddFlag( m_FlagAverageOnly);
      // flags
      sp_cmd->AddFlag( m_FlagDataSetRetriever);

      // flag for model storage
      sp_cmd->AddFlag( m_FlagModelInterfaceStorage);
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
    void ModelTest::WritePredictions
    (
      const descriptor::Dataset &DATA,
      const linal::MatrixConstInterface< float> &PREDICTIONS,
      const std::string &FILENAME
    ) const
    {
      // get predictions
      const model::FeatureDataSetInterface< float> &experimental( *( DATA.GetResultsPtr()));

      // write experimental and predicted data
      io::OFStream output;
      io::File::MustOpenOFStream( output, FILENAME);
      if( !m_FlagIDCode->GetFlag())
      {
        // create matrix with twice the number of cols of the result vector
        linal::Matrix< float> matrix
        (
          experimental.GetNumberFeatures(), size_t( 2) * experimental.GetFeatureSize(), float( 0.0)
        );

        // off set to access indices of matrix for result values
        const size_t offset( experimental.GetFeatureSize());

        for( size_t data_id( 0), data_set_size( matrix.GetNumberRows()); data_id < data_set_size; ++data_id)
        {
          float *matrix_row( matrix[ data_id]);
          for( size_t row_element_id( 0), row_size( experimental.GetFeatureSize()); row_element_id < row_size; ++row_element_id)
          {
            // get the actual result
            matrix_row[ row_element_id] = experimental[ data_id][ row_element_id];
            matrix_row[ row_element_id + offset] = PREDICTIONS( data_id, row_element_id);
          }
        }

        // output the matrix directly
        output << matrix;
      }
      else
      {
        // output experimental, predicted, id
        output << "Experimental\tPredicted\t\"" << m_FlagIDCode->GetFirstParameter()->GetValue() << "\"\n";
        const linal::MatrixConstReference< char> ids( DATA.GetIdsReference());
        const size_t feature_size( experimental.GetFeatureSize());
        const size_t id_size( ids.GetNumberCols());
        for( size_t data_id( 0), data_set_size( experimental.GetNumberFeatures()); data_id < data_set_size; ++data_id)
        {
          for( size_t row_element_id( 0); row_element_id < feature_size; ++row_element_id)
          {
            // get the actual result
            output << experimental[ data_id][ row_element_id] << '\t';
          }
          for( size_t row_element_id( 0); row_element_id < feature_size; ++row_element_id)
          {
            output << PREDICTIONS( data_id, row_element_id) << '\t';
          }
          output << std::string( ids.GetRow( data_id).Begin(), id_size) << '\n';
        }
      }
      io::File::CloseClearFStream( output);
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ModelTest::Main() const
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
      linal::Matrix< float> avg_predictions;

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
        model->Rescale( *features_ptr);

        const model::FeatureDataSet< float> predictions( model->operator ()( *data_set->GetFeaturesPtr()));

        // initialize or add to average predictions matrix
        if( avg_predictions.GetNumberRows() > 0)
        {
          // add prediction of model to average predictions
          avg_predictions += predictions.GetMatrix();
        }
        else
        {
          avg_predictions = predictions.GetMatrix();
        }

        if( !m_FlagAverageOnly->GetFlag())
        {
          BCL_MessageStd
          (
            "predictions with model id " + util::Format()( *itr_model_ids)
            + " #" + util::Format()( predictions.GetNumberFeatures())
            + " predicted columns: " + util::Format()( predictions.GetFeatureSize())
          );

          const std::string filename
          (
            m_FlagOutputBase->GetParameterList().FirstElement()->GetValue()
            + "."
            + util::Format()( *itr_model_ids) + ".dat"
          );

          // write predictions for individual model
          WritePredictions( *data_set, predictions.GetRawMatrix(), filename);
        }
      }

      if( m_FlagOutputBase->GetFlag() && ( model_ids.GetSize() > 1 || m_FlagAverageOnly->GetFlag()))
      {
        // prepare output
        std::string filename( m_FlagOutputBase->GetParameterList().FirstElement()->GetValue());
        if( !m_FlagAverageOnly->GetFlag())
        {
          filename += ".avg.txt";
        }

        // obtain matrix with average predictions by dividing by the number of evaluated models
        avg_predictions /= float( model_ids.GetSize());

        // write predictions for individual model
        common_info.SetIndependentPredictionsFilename( filename);
        common_info.WritePredictions
        (
          data_set->GetResultsReference(),
          avg_predictions,
          data_set->GetIdsReference(),
          filename
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
    ModelTest::ModelTest() :
      m_FlagOutputBase
      (
        new command::FlagStatic
        (
          "output",
          "flag for output base name (if -average is also given, this will be the actual filename instead)",
          command::Parameter( "filename", "flag for output of screening", "prediction.txt")
        )
      ),
      m_FlagAverageOnly
      (
        new command::FlagStatic
        (
          "average",
          "set this flag to only write out the average file (to -output), rather than files for each model"
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

    const ApplicationType ModelTest::ModelTest_Instance
    (
      GetAppGroups().AddAppToGroup( new ModelTest(), GetAppGroups().e_Model)
    );

  } // namespace app
} // namespace bcl
