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
#include "bcl_app_model_train.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_binary_function_bind_second.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_approximator_base.h"
#include "model/bcl_model_feature_data_reference.h"
#include "model/bcl_model_meta_data_storage_interface.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "model/bcl_model_store_interface.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_convergence_result.h"
#include "opti/bcl_opti_criterion_elapsed_time.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_result_threshold.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "opti/bcl_opti_printer_default.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ModelTrain::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // parameters

      // filename for iterate object
      sp_cmd->AddParameter( m_Iterator);
      // whether to suppress output objective function results while training
      sp_cmd->AddFlag( m_SuppressProgressOutput);
      // maximum time, in minutes; the training is stopped at the next iteration thereafter
      sp_cmd->AddFlag( m_FlagTerminateAfterTimeMinutes);
      // add flag for max iterations, a common termination criteria
      sp_cmd->AddFlag( m_FlagMaxNumberIterations);
      // add flag for max_iterations_since_last_improvement
      sp_cmd->AddFlag( m_FlagMaxIterationsWithoutImprovement);
      // add flag for result averaging over rounds to smooth out noisy objective functions
      sp_cmd->AddFlag( m_FlagTrackerResultAverageRounds);

      // add flag for final objective function
      sp_cmd->AddFlag( m_FlagFinalObjectiveFunction);

      // flag for writing out a serialized iterate instance for continued training
      sp_cmd->AddFlag( m_FlagContinuedTraining);
      sp_cmd->AddFlag( m_FlagReadApproximator);

      // training, monitoring, independent
      sp_cmd->AddFlag( m_FlagTrainingDataSet);
      sp_cmd->AddFlag( m_FlagMonitoringDataSet);
      sp_cmd->AddFlag( m_FlagIndependentDataSet);

      // code object files
      sp_cmd->AddFlag( m_FlagFeatureCode);
      sp_cmd->AddFlag( m_FlagResultCode);
      sp_cmd->AddFlag( m_FlagIdCode);

      // flag for printing final predictions
      sp_cmd->AddFlag( m_FlagPrintTrainingPredictions);
      sp_cmd->AddFlag( m_FlagPrintMonitoringPredictions);
      sp_cmd->AddFlag( m_FlagPrintIndependentPredictions);

      // flag for model::Interface storage
      sp_cmd->AddFlag( m_FlagModelInterfaceStorage);

      // flag for model::MetaData storage used in descriptor selection
      sp_cmd->AddFlag( m_FlagMetaDataStorage);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief the Main function
    //! @return error code - 0 for success
    int ModelTrain::Main() const
    {
      // s_instances for serialization
      math::BinaryFunctionBindSecond< float, float, bool>::s_Instance.IsDefined();
      model::FeatureDataSet< float>::s_Instance.IsDefined();
      opti::CriterionConvergenceResult< util::ShPtr< model::Interface>, float>::s_Instance.IsDefined();
      opti::CriterionResultThreshold< util::ShPtr< model::Interface>, float>::s_Instance.IsDefined();
      storage::Pair< util::ShPtr< model::Interface>, float>::s_Instance.IsDefined();
      storage::Vector< int>::s_Instance.IsDefined();
      storage::Vector< size_t>::s_Instance.IsDefined();

      BCL_Assert
      (
        m_FlagTerminateAfterTimeMinutes->GetFlag()
        || m_FlagMaxNumberIterations->GetFlag()
        || m_FlagMaxIterationsWithoutImprovement->GetFlag(),
        "At least one termination criteria must be given (-max_minutes | -max_iterations | -max_unimproved_iterations)"
      );

      // stopwatch
      util::Stopwatch timer;
      timer.Start();

      // initialize monitor data
      util::Implementation< model::RetrieveDataSetBase> monitoring_retriever
      (
        util::ObjectDataLabel( m_FlagMonitoringDataSet->GetFirstParameter()->GetValue())
      );
      monitoring_retriever->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);
      monitoring_retriever->SelectResultsGivenFilenameFlag( *m_FlagResultCode);
      monitoring_retriever->SelectIdsGivenFilenameFlag( *m_FlagIdCode);

      // storage for the final model::Interface
      util::Implementation< model::StoreInterface> model_storage;
      if( m_FlagModelInterfaceStorage->GetFlag())
      {
        // store the result code; asserts if it cannot be created or is different from the existing result code
        model_storage =
          util::Implementation< model::StoreInterface>( m_FlagModelInterfaceStorage->GetFirstParameter()->GetValue());
        model_storage->StoreResultDescriptor( monitoring_retriever->GetResultCodeWithSizes().GetLabel());
      }

      // create sh-ptrs to the different datasets
      util::ShPtr< descriptor::Dataset> monitoring_data_set, training_data_set, independent_data_set;

      // generate the monitoring dataset
      monitoring_data_set = monitoring_retriever->GenerateDataSet();

      // check whether monitoring dataset is not empty
      BCL_Assert( !monitoring_data_set->IsEmpty(), "Monitoring dataset is empty!");

      // initialize independent data
      util::Implementation< model::RetrieveDataSetBase> independent_retriever
      (
        util::ObjectDataLabel( m_FlagIndependentDataSet->GetFirstParameter()->GetValue())
      );

      // set up the feature result codes
      independent_retriever->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);
      independent_retriever->SelectResultsGivenFilenameFlag( *m_FlagResultCode);
      independent_retriever->SelectIdsGivenFilenameFlag( *m_FlagIdCode);

      if( independent_retriever.GetLabel() == monitoring_retriever.GetLabel())
      {
        independent_data_set = monitoring_data_set;
      }
      else
      {
        independent_data_set = independent_retriever->GenerateDataSet();
      }

      // check whether monitoring dataset is not empty
      BCL_Assert( !independent_data_set->IsEmpty(), "Independent dataset is empty!");

      // initialize training data
      util::Implementation< model::RetrieveDataSetBase> training_retriever
      (
        util::ObjectDataLabel( m_FlagTrainingDataSet->GetFirstParameter()->GetValue())
      );

      // set up the feature result codes
      training_retriever->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);
      training_retriever->SelectResultsGivenFilenameFlag( *m_FlagResultCode);
      training_retriever->SelectIdsGivenFilenameFlag( *m_FlagIdCode);
      if( training_retriever.GetLabel() == monitoring_retriever.GetLabel())
      {
        training_data_set = monitoring_data_set;
      }
      else if( training_retriever.GetLabel() == independent_retriever.GetLabel())
      {
        training_data_set = independent_data_set;
      }
      else
      {
        training_data_set = training_retriever->GenerateDataSet();
      }

      // check whether monitoring dataset is not empty
      BCL_Assert( !training_data_set->IsEmpty(), "Training dataset is empty!");

      BCL_MessageStd
      (
        "Created independent data set with " + util::Format()( independent_data_set->GetSize())
        + " points\nCreated monitoring data set with " + util::Format()( monitoring_data_set->GetSize())
        + " points\nCreated training data set with " + util::Format()( training_data_set->GetSize())
        + " points"
      );

      BCL_MessageStd
      (
        "\nNumber of descriptor values: " + util::Format()( independent_data_set->GetFeatureSize())
        + "\nNumber of result values:  " + util::Format()( independent_data_set->GetResultSize())
      );

      BCL_MessageDbg( "read in dataset finished ... done");

      BCL_MessageStd
      (
        "constructing datasets finished total: " + util::Format()( timer.GetTotalTime().GetSeconds()) + " [sec]"
      );

      // reset timer
      timer.Reset();
      timer.Start();

      util::ShPtr< model::ApproximatorBase> sp_iterate;

      if( m_FlagReadApproximator->GetFlag())
      {
        if( m_Iterator->GetWasSetInCommandLine())
        {
          BCL_MessageCrt( "Iterate parameter was ignored -- reading value from given file");
        }
        // initialize input stream
        io::IFStream read;
        io::File::MustOpenIFStream( read, m_FlagReadApproximator->GetFirstParameter()->GetValue());
        io::Serialize::Read( sp_iterate, read);
        io::File::CloseClearFStream( read);
        sp_iterate->SetTrainingContinued( true);
      }
      else if( m_Iterator->GetWasSetInCommandLine())
      {
        // READ IN ITERATE FUNCTION
        util::Implementation< model::ApproximatorBase> iterate_impl( m_Iterator->GetValue());
        // convert to a shared pointer because approximator requires it
        sp_iterate = ConvertToShPtr( iterate_impl);
      }

      BCL_Assert
      (
        sp_iterate.IsDefined(),
        "Undefined model iterator!"
      );

      // set id, feature, and result code objects
      sp_iterate->SetIds( training_retriever->GetIdCodeWithSizes());
      sp_iterate->SetFeatures( training_retriever->GetFeatureLabelsWithSizes());
      sp_iterate->SetResults( training_retriever->GetResultCodeWithSizes());

      // set termination criterion.  This is done here in case the model does extensive work in SetTrainingData that
      // could significantly cut into the training time.  For example, neural network pre-training.  We always want to
      // respect the max-minutes flag, if possible

      // set termination
      sp_iterate->SetCriterion( CreateTerminationCriterion());

      // invoke dataset preparation
      sp_iterate->SetTrainingData( training_data_set);

      BCL_MessageDbg( "read in iterate ... done");

      // READ IN OBJECTIVE FUNCTION
      util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function( sp_iterate->GetObjectiveFunction());

      sp_objective_function->SetData
      (
        monitoring_data_set,
        sp_iterate->GetRescaleFeatureDataSet()
      );

      BCL_MessageDbg( "read in objective function ... done");

      // FINAL OBJECTIVE FUNCTION
      util::ShPtr< model::ObjectiveFunctionWrapper> sp_final_obj_function;

      // if final objective function is set by flag
      if( m_FlagFinalObjectiveFunction->GetFlag())
      {
        sp_final_obj_function = util::ShPtr< model::ObjectiveFunctionWrapper>
        (
          new model::ObjectiveFunctionWrapper
          (
            util::Implementation< model::ObjectiveFunctionInterface>
            (
              util::ObjectDataLabel( m_FlagFinalObjectiveFunction->GetFirstParameter()->GetValue())
            )
          )
        );
      }

      BCL_MessageDbg( "read in final objective function ... done");

      // Create printer
      sp_iterate->SetPrinter
      (
        opti::PrinterDefault< util::ShPtr< model::Interface>, float>
        (
          false,     // never write model
          !m_SuppressProgressOutput->GetFlag() // write result unless user has specifically requested otherwise
        )
      );
      sp_iterate->GetTracker().SetAveragingWindowSize
      (
        m_FlagTrackerResultAverageRounds->GetFirstParameter()->GetNumericalValue< size_t>()
      );

      BCL_MessageDbg( "build approximator finished ... done");

      // call minimizer till criterion met
      sp_iterate->Approximate();

      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > result;
      if( m_FlagMaxNumberIterations->GetFirstParameter()->GetNumericalValue< size_t>() == 0)
      {
        result = util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
        (
          new storage::Pair< util::ShPtr< model::Interface>, float>( sp_iterate->GetCurrentModel(), 0.0)
        );
      }
      else
      {
        result = sp_iterate->GetTracker().GetBest();
      }

      // reference to current model
      const model::Interface &model( *result->First());

      BCL_MessageStd
      (
        "approximator training total: " + util::Format()( timer.GetTotalTime().GetSeconds()) + " [sec]"
      );
      // reset timer
      timer.Reset();
      timer.Start();

      // apply objective function on final model
      float score( result->Second());

      const util::ShPtr< model::Interface> &model_sp( result->First());

      BCL_MessageStd( "monitoring dataset score:" + util::Format()( score));

      // if the user wants the iterate written out, write it out
      if( m_FlagContinuedTraining->GetFlag())
      {
        io::OFStream output;
        if( io::File::TryOpenOFStream( output, m_FlagContinuedTraining->GetFirstParameter()->GetValue()))
        {
          io::Serialize::Write( sp_iterate, output);
          io::File::CloseClearFStream( output);
        }
        else
        {
          BCL_MessageCrt
          (
            "Could not open: " + m_FlagContinuedTraining->GetFirstParameter()->GetValue() + " for writing"
          );
        }
      }

      // create all the cross-validation information
      model::CrossValidationInfo cv_info
      (
        score,
        sp_final_obj_function.IsDefined()
        ? sp_final_obj_function->GetImprovementType()
        : opti::e_LargerEqualIsBetter,
        independent_retriever.GetLabel(),
        monitoring_retriever.GetLabel(),
        training_retriever.GetLabel(),
        sp_final_obj_function.IsDefined()
        ? sp_final_obj_function->GetImplementation().GetLabel()
        : util::ObjectDataLabel(),
        sp_iterate->GetLabel(),
        m_FlagPrintIndependentPredictions->GetFlag()
        ? m_FlagPrintIndependentPredictions->GetFirstParameter()->GetValue()
        : "",
        independent_retriever->GetIdCodeWithSizes(),
        model_sp->GetNumberOutputs()
      );

      if
      (
        m_FlagFinalObjectiveFunction->GetFlag()
        || m_FlagMetaDataStorage->GetFlag()
        || m_FlagModelInterfaceStorage->GetFlag()
        || m_FlagIndependentDataSet->GetFlag()
      )
      {
        if( sp_iterate->GetRescaleFeatureDataSet().IsDefined())
        {
          // rescale the dataset to obviate the need to copy and then rescale the dataset internally
          independent_data_set->GetFeatures().Rescale( *sp_iterate->GetRescaleFeatureDataSet());
        }

        // create independent predictions
        model::FeatureDataSet< float> independent_predictions
        (
          model.PredictWithoutRescaling( *independent_data_set->GetFeaturesPtr()).DeScale()
        );

        if( m_FlagFinalObjectiveFunction->GetFlag())
        {
          // set independent data for final objective function
          sp_final_obj_function->SetData
          (
            independent_data_set,
            sp_iterate->GetRescaleFeatureDataSet()
          );

          // the final score is determined by final objective function on independent dataset
          score = sp_final_obj_function->Evaluate( independent_predictions);

          BCL_MessageStd
          (
            "independent data: final score ("
            + sp_final_obj_function->GetImplementation().GetAlias()
            + "): " + util::Format()( score)
            + " time: " + util::Format()( timer.GetTotalTime().GetSeconds()) + " [sec]"
          );
        }

        if( m_FlagModelInterfaceStorage->GetFlag())
        {
          const std::string model_key
          (
            model_storage->Store
            (
              model_sp,
              score,
              training_retriever->GetFeatureLabelsWithSizes().GetLabel(),
              m_Iterator->GetValue(),
              m_FlagFinalObjectiveFunction->GetFirstParameter()->GetValue(),
              cv_info
            )
          );

          BCL_MessageStd( "Final model written (key = " + model_key + ") ... done");
        }

        if( m_FlagMetaDataStorage->GetFlag())
        {
           // storage for the final model::Interface
          util::Implementation< model::MetaDataStorageInterface> meta_storage
          (
            m_FlagMetaDataStorage->GetFirstParameter()->GetValue()
          );

          meta_storage->Store
          (
            score,
            training_retriever->GetFeatureLabelsWithSizes().GetLabel(),
            m_Iterator->GetValue(),
            m_FlagFinalObjectiveFunction->GetFirstParameter()->GetValue()
          );

          BCL_MessageStd( "Descriptor selection data written ... done");
        }

        if( m_FlagPrintIndependentPredictions->GetFlag())
        {
          cv_info.WritePredictions
          (
            independent_data_set->GetResults(),
            independent_predictions,
            independent_data_set->GetIds(),
            m_FlagPrintIndependentPredictions->GetFirstParameter()->GetValue()
          );
          BCL_MessageStd( "Independent Predictions written ... done");
        }
      }

      if( m_FlagPrintTrainingPredictions->GetFlag())
      {
        model::FeatureDataSet< float> predictions( model( *( training_data_set->GetFeaturesPtr())));
        cv_info.WritePredictions
        (
          training_data_set->GetResults(),
          predictions,
          training_data_set->GetIds(),
          m_FlagPrintTrainingPredictions->GetFirstParameter()->GetValue()
        );
        BCL_MessageStd( "Training Predictions written ... done");
      }

      if( m_FlagPrintMonitoringPredictions->GetFlag())
      {
        model::FeatureDataSet< float> predictions( model( *( monitoring_data_set->GetFeaturesPtr())));
        cv_info.WritePredictions
        (
          monitoring_data_set->GetResults(),
          predictions,
          monitoring_data_set->GetIds(),
          m_FlagPrintMonitoringPredictions->GetFirstParameter()->GetValue()
        );
        BCL_MessageStd( "Monitoring Predictions written ... done");
      }

      // end
      return 0;
    } // Main

    //! @brief Create termination criterion from command line flag
    opti::CriterionCombine< util::ShPtr< model::Interface>, float> ModelTrain::CreateTerminationCriterion() const
    {
      // Create terminate
      opti::CriterionCombine< util::ShPtr< model::Interface>, float> termination_criteria;

      if( m_FlagMaxNumberIterations->GetFlag())
      {
        // just use the max number iterations termination criteria
        termination_criteria.InsertCriteria
        (
          opti::CriterionNumberIterations< util::ShPtr< model::Interface>, float>
          (
            m_FlagMaxNumberIterations->GetFirstParameter()->GetNumericalValue< size_t>()
          )
        );
      }
      if( m_FlagTerminateAfterTimeMinutes->GetFlag())
      {
        const double max_time_in_minutes( m_FlagTerminateAfterTimeMinutes->GetFirstParameter()->GetNumericalValue< double>());

        termination_criteria.InsertCriteria
        (
          opti::CriterionElapsedTime< util::ShPtr< model::Interface>, float>
          (
            util::Time( size_t( max_time_in_minutes * 60.0), size_t( 0))
          )
        );
      }
      if( m_FlagMaxIterationsWithoutImprovement->GetFlag())
      {
        const size_t max_unimproved_iterations
        (
          m_FlagMaxIterationsWithoutImprovement->GetFirstParameter()->GetNumericalValue< double>()
        );
        termination_criteria.InsertCriteria
        (
          opti::CriterionUnimproved< util::ShPtr< model::Interface>, float>( max_unimproved_iterations)
        );
      }
      return termination_criteria;
    }

    //! @brief standard constructor
    ModelTrain::ModelTrain() :
      m_Iterator
      (
        new command::Parameter
        (
          "iterator",
          "data label of a specific iterate used as training algorithm in opti::Approximator",
          command::ParameterCheckSerializable( util::Implementation< model::ApproximatorBase>()),
          "LinearRegression"
        )
      ),
      m_SuppressProgressOutput
      (
        new command::FlagStatic
        (
          "suppress_progress_output",
          "by default, objective function evaluations are shown as the model is training; set this flag to suppress that output"
        )
      ),
      m_FlagTerminateAfterTimeMinutes
      (
        new command::FlagStatic
        (
          "max_minutes",
          "maximum # of minutes to train; if # of iterations is reached first, iteration will stop then",
          command::Parameter
          (
            "minutes",
            "",
            command::ParameterCheckRanged< double>( 0.5, double( 60 * 24 * 365)),
            "1440"
          )
        )
      ),
      m_FlagMaxNumberIterations
      (
        new command::FlagStatic
        (
          "max_iterations",
          "maximum number of iterations",
          command::Parameter
          (
            "iterations",
            "",
            command::ParameterCheckRanged< size_t>(),
            util::Format()( std::numeric_limits< size_t>::max())
          )
        )
      ),
      m_FlagMaxIterationsWithoutImprovement
      (
        new command::FlagStatic
        (
          "max_unimproved_iterations",
          "maximum number of iterations that can pass between improvement steps without stopping the training",
          command::Parameter
          (
            "iterations",
            "",
            command::ParameterCheckRanged< size_t>(),
            util::Format()( std::numeric_limits< size_t>::max())
          )
        )
      ),
      m_FlagTrackerResultAverageRounds
      (
        new command::FlagStatic
        (
          "result_averaging_window",
          "Window size for computing the current average result, helps smooth noisy objective functions. Example values:"
          "0 - Choose the last model; 1 - Choose the best model on the monitoring dataset; 2 - Choose the best model "
          "based on the last two iteration's objective functions on the monitoring dataset according to 1/3 Last round "
          "objective function + 2/3 This round objective function. Higher values consider additional previous "
          "rounds according to triangularly-weighted average ",
          command::Parameter
          (
            "window_size",
            "",
            command::ParameterCheckRanged< size_t>(),
            "1"
          )
        )
      ),
      m_FlagFinalObjectiveFunction
      (
        new command::FlagStatic
        (
          "final_objective_function",
          "data label for an objective function the evaluates the final model",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::ObjectiveFunctionInterface>()),
            "RMSD"
          )
        )
      ),
      m_FlagContinuedTraining
      (
        new command::FlagStatic
        (
          "continued_training",
          "write out a serialized iterate instance for later continuation of training (using the -continue flag)",
          command::Parameter
          (
            "file_name",
            "file name of written out iterate instance at the conclusion of training",
            ""
          )
        )
      ),
      m_FlagReadApproximator
      (
        new command::FlagStatic
        (
          "continue",
          "read in an iterator for training",
          command::Parameter
          (
            "file_name",
            "file name of written out iterate instance",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_FlagTrainingDataSet
      (
        new command::FlagStatic
        (
          "training",
          "source of dataset used to train the model",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_FlagMonitoringDataSet
      (
        new command::FlagStatic
        (
          "monitoring",
          "source of dataset used when deciding whether model has improved",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_FlagIndependentDataSet
      (
        new command::FlagStatic
        (
          "independent",
          "source of dataset used when evaluating the final objective function",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_FlagFeatureCode
      (
        new command::FlagStatic
        (
          "feature_labels",
          "label or file containing the label for the feature descriptors",
          command::Parameter( "feature_labels", "", "")
        )
      ),
      m_FlagResultCode
      (
        new command::FlagStatic
        (
          "result_labels",
          "label or file containing the label for the result descriptors",
          command::Parameter( "result_labels", "", "")
        )
      ),
      m_FlagIdCode
      (
        new command::FlagStatic
        (
          "id_labels",
          "label or file containing the label for the id descriptors",
          command::Parameter( "id_labels", "", "")
        )
      ),
      m_FlagPrintTrainingPredictions
      (
        new command::FlagStatic
        (
          "print_training_predictions",
          "file to print the predicted values for training data to (1st column = actual, 2nd column = predicted)",
          command::Parameter( "output_filename", "", "training_pred.txt")
        )
      ),
      m_FlagPrintMonitoringPredictions
      (
        new command::FlagStatic
        (
          "print_monitoring_predictions",
          "file to print the predicted values for monitoring data to (1st column = actual, 2nd column = predicted)",
          command::Parameter( "output_filename", "", "monitor_pred.txt")
        )
      ),
      m_FlagPrintIndependentPredictions
      (
        new command::FlagStatic
        (
          "print_independent_predictions",
          "file to print the predicted values for independent data to (1st column = actual, 2nd column = predicted)",
          command::Parameter( "output_filename", "", "ind_pred.txt")
        )
      ),
      m_FlagModelInterfaceStorage
      (
        new command::FlagStatic
        (
          "storage_model",
          "choice of model storage",
          command::Parameter
          (
            "model_storage",
            "choice of model storage",
            command::ParameterCheckSerializable( util::Implementation< model::StoreInterface>()),
            "File"
          )
        )
      ),
      m_FlagMetaDataStorage
      (
        new command::FlagStatic
        (
          "storage_descriptor_selection",
          "choice of storage for meta data (used in descriptor selection)",
          command::Parameter
          (
            "storage",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::MetaDataStorageInterface>()),
            "File"
          )
        )
      )
    {
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> ModelTrain::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "TrainModel");
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ModelTrain::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL model:Train, terms of use, appropriate citation, installation "
        "procedures, BCL model:Train execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL model:Train?\n"
        "BCL model:Train is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons. BCL model:Train trains a "
        " machine learning model. Available are implementations for artificial neural networks (ANN),"
        "support vector machines (SVM), kohonen networks (KN) and decision trees (DT). Additionally, GPU acceleration is"
        "available for ANNs and SVMs to leverage the power of your graphics card. "
        "The training is applied on a given data set which is split into monitoring, independent, and training "
        "partitions. While training the machine learning technique is done iteratively, an specified objective function is "
        "evaluated to determine model performance improvement. BCL model:Train is set up to perform cross-validation."
        "The application is highly flexible in terms of trained model evaluation and serialization of model predictions"
        "for each of the monitoring, independent, and training partitions as well as writing out the model itself."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL model:Train.\n"
        "When using BCL model:Train in a publication, please cite the following publication describing the application's "
        "development:\n"
        "\n"
        "Butkiewicz, M.; Lowe, E.W., Jr.; Mueller, R.; Mendenhall, J.L.; Teixeira, P.L.; Weaver, C.D.; Meiler, J.\n"
        "'Benchmarking Ligand-Based Virtual High-Throughput Screening with the PubChem Database'. Molecules 2013, 18, 735-756.\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL model:Train.\n"
        "Running BCL model:Train consists of the following steps.\n"
        "\n"
        "1) Generate a data set (.bin file) with BCL descriptor:GenerateDataSet\n"
        "\n"
        "2) Run BCL model:Train on one or multiple <experimental predicted value file>s\n"
        "\n"
        "3) Obtain a gnuplot plot of your chosen quality measures or a correlation plot.\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing: <bcl.exe> model:Train -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        "For more general information about the product, type <bcl.exe>  model:Train -help\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL model:Train.\n"
        "BCL model:Train is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE\n"
        "Provide two data set for your active and inactive molecules with BCL descriptor:GenerateDataSet. "
        "The result will be two files: actives.bin and inactives.bin\n"
        "\n"
        "Execute the following commandline to train a neural network:\n"
        "<bcl.exe> model:Train "
        "  'NeuralNetwork( transfer function = Sigmoid, weight update = Simple(eta=0.1,alpha=0.1),"
        " objective function = RMSD, steps per update=1, hidden architecture(8))' \n"
        "  -max_minutes 1\n"
        "  -max_iterations 100\n"
        "  -feature_labels features.txt\n"
        "  -opencl Disable\n"
        "  -final_objective_function 'RMSD' \n"
        "  -training 'Balanced( Subset(number chunks=3,chunks=\"[0,3)-[0]-[1]\",filename=\"actives.bin\")),"
        "Subset(number chunks=3,chunks=\"[0,3)-[0]-[1]\",filename=\"inactives.bin\"))))'\n"
        "  -monitoring 'Combined(Subset(number chunks=3,chunks=\"[1]\", filename = \"actives.bin\")),"
        "Subset(number chunks=3,chunks=\"[1]\", filename = \"inactives.bin\"))))'\n"
        "  -independent 'Combined(Subset(number chunks=3,chunks=\"[0]\",filename=\"actives.bin\"),"
        "Subset(number chunks=3,chunks=\"[0]\",filename=\"inactives.bin\"))'\n"
        "  -print_independent_predictions ./independent.txt.gz \n"
        "  -storage_model 'File(directory=./models/)' \n"
        "  -logger File ./log.txt\n"
        "\n"
        "To get further information on any of the flag parameter having a bracket notation '()' you can place the key word "
        "'help' in the brackets and obtain more specific information about the flag parameter."
        "Examples:\n"
        "  - get all available options of a NeuralNetwork with : <bcl.exe> NeuralNetwork(help)\n"
        "  - get all available options of the Balanced option for the -training flag: <bcl.exe> -training Balanced(help)\n"
        + DefaultSectionSeparator()
      );

      return readme;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &ModelTrain::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::ModelTrain: Sophisticated Machine Learning \n\n"
        "As the center-piece application in the application suite BCL::ChemInfo, BCL::ModelTrain paves the road to"
        "sophisticated machine learning modeling. It is an important step for creating quantitative structure activity "
        "relationship (QSAR) models.\n\n"
        "The aim of this project is provide an flexible infrastructure to process a given data set in various "
        "formats (eg. csv, sdf) and allow for machine learning. Designed to handle large high-throughput screening data"
        " sets, protein related data set, or plain numeric data sets BCL::ModelTrain allows cross-validation training "
        "through data set partitioning flexibility. A large number of machine learning algorithms implementations "
        "is easily accessible:\n"
        "\n"
        "Artificial neural networks\n"
        "Support vector machines with extension for regression\n"
        "Kohonen maps (self-organizing maps)\n"
        "Decision tree\n"
        "kappa-Nearest Neighbor\n"
        "Restricted Boltzmann Machines\n"
        "Auto Encoders\n"
        "\n"
        "For a subset of machine learning algorithms a GPU version is available that leverage the vast computational"
        "power of your graphics card(s) and allows for significant performance speed-ups!\n"
        "\n"
        "!bcl_model_train_workflow.png!\n"
        "\nFig. 1:\n"
        "\n\n"
        "BCL::ModelTrain is placed in context of the BCL::ChemInfo suite workflow and a subset of "
        "machine learning algorithms are highlighted.\n"
        ""
      );
      return s_web_text;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ModelTrain::GetDescription() const
    {
      return "Train any machine learning model, including ANNs, SVMs, and many more";
    }

    const ApplicationType ModelTrain::ModelTrain_Instance
    (
      GetAppGroups().AddAppToGroup( new ModelTrain(), GetAppGroups().e_Model)
    );

  } // namespace app
} // namespace bcl
