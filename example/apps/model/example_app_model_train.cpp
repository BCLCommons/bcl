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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "model/bcl_app_model_train.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_directory.h"
#include "sched/bcl_sched_schedulers.h"
#include "util/bcl_util_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_train_model.cpp
  //! @brief application example for application app::TrainModel.
  //!
  //! @author butkiem1
  //! @date Mar 03, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppModelTrain :
    public ExampleInterface
  {
  public:

    ExampleAppModelTrain *Clone() const
    {
      return new ExampleAppModelTrain( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    void SetTrainingFlagsBinFiles
    (
      ApplicationExampleHelper &APP_HELPER,
      const std::string &ACTIVES_FILENAME,
      const std::string &INACTIVES_FILENAME,
      const std::string &DESCRIPTOR_LABELS_FILENAME
    ) const
    {
      APP_HELPER.SetFlag
      (
        "training",
        "Combined"
        "( "
        " Subset( filename=" + ACTIVES_FILENAME + ",number chunks=10,chunks=\"[0,10)-[1]-[0]\"),"
        " Subset( filename=" + INACTIVES_FILENAME + ",number chunks=10,chunks=\"[0,10)-[1]-[0]\")"
        ")"
      );
      // specify independent data set
      APP_HELPER.SetFlag
      (
        "independent",
        "Combined"
        "( "
        " Subset( filename=" + ACTIVES_FILENAME + ",number chunks=10,chunks=[0]),"
        " Subset( filename=" + INACTIVES_FILENAME + ",number chunks=10,chunks=[0])"
        ")"
      );
      // specify monitoring data set
      APP_HELPER.SetFlag
      (
        "monitoring",
        "Combined"
        "( "
        " Subset( filename=" + ACTIVES_FILENAME + ",number chunks=10,chunks=[1]),"
        " Subset( filename=" + INACTIVES_FILENAME + ",number chunks=10,chunks=[1])"
        ")"
      );

      // specify maximal number of approximator iterations
      APP_HELPER.SetFlag( "max_iterations", "3");

      // specify descriptor labels
      APP_HELPER.SetFlag( "feature_labels", DESCRIPTOR_LABELS_FILENAME);

      // specify objective function to evaluate final final model
      APP_HELPER.SetFlag
      (
        "final_objective_function",
        "EnrichmentAverage( cutoff = 4, enrichment max = 0.05, step size = 0.005, parity = 1)"
      );

      // use scheduler with pthreads
      if( sched::Scheduler( "PThread").IsDefined())
      {
        APP_HELPER.SetFlag( "scheduler", storage::Vector< std::string>::Create( "PThread", "6"));
      }
    }

    void SetTrainingFlagsBinFilesDatFile
    (
      ApplicationExampleHelper &APP_HELPER,
      const std::string &TRAINING_DATA_FILENAME
    ) const
    {
      APP_HELPER.SetFlag
      (
        "training",
        "File"
        "( "
        "  filename=" + TRAINING_DATA_FILENAME + ",number chunks=10,chunks=\"[0,10)-[1]-[0]\""
        ")"
      );
      // specify independent data set
      APP_HELPER.SetFlag
      (
        "independent",
        "File"
        "( "
        "  filename=" + TRAINING_DATA_FILENAME + ",number chunks=10,chunks=[0]"
        ")"
      );
      // specify monitoring data set
      APP_HELPER.SetFlag
      (
        "monitoring",
        "File"
        "( "
        "  filename=" + TRAINING_DATA_FILENAME + ",number chunks=10,chunks=[1]"
        ")"
      );

      // specify maximal number of approximator iterations
      APP_HELPER.SetFlag( "max_iterations", "3");

      // specify objective function to evaluate final final model
      APP_HELPER.SetFlag
      (
        "final_objective_function",
        "EnrichmentAverage( cutoff = 4, enrichment max = 0.05, step size = 0.005, parity = 1)"
      );

      // use scheduler with pthreads
      APP_HELPER.SetFlag( "scheduler", storage::Vector< std::string>::Create( "PThread", "6"));
    }

    int Run() const
    {

    ////////////////
    // unit tests //
    ////////////////

      // unit tests of application class functions, if applicable, go here

    ////////////////////
    // initialization //
    ////////////////////

      // create a helper to run this application
      ApplicationExampleHelper train_model( app::ModelTrain::ModelTrain_Instance);

    ///////////
    // files //
    ///////////

      // create filenames for any input/output files used/generated by the test

      // path to bin file for active compounds
      const std::string actives_bin
      (
        train_model.GetThisApplicationExampleInputPath() + "aid891_actives.bin"
      );

      // path to bin file for inactive compounds
      const std::string inactives_bin
      (
        train_model.GetThisApplicationExampleInputPath() + "aid891_inactives.bin"
      );

      // path to dat file with all datapoints of interest
      const std::string data_points_dat
      (
        train_model.GetThisApplicationExampleInputPath() + "descriptors_10.dat"
      );

      // file containing descriptors stored in bin file for each molecule
      const std::string filename_descriptor_labels
      (
        train_model.GetThisApplicationExampleInputPath() + "code_input.object"
      );

    ///////////////////////
    // integration tests //
    ///////////////////////

    /////////////////
    // SVM example //
    /////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "SupportVectorMachine"
        "( "
        " objective function=EnrichmentAverage( cutoff = 4, enrichment max = 0.05, step size = 0.005, parity = 1),"
        " kernel = RBF( gamma=0.4 ),"
        " iterations=100,"
        " cost=0.1,"
        " gap_threshold=0.1"
        ")"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    //////////////////////////////////////////////
    // ANN Simple propagation non-batch example //
    //////////////////////////////////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "NeuralNetwork"
        "("
        " transfer function = Sigmoid,"
        " weight update = Simple(eta=0.1,alpha=0.5),"
        " objective function = RMSD,"
        " steps per update=1,"
        " hidden architecture(8)"
        ")"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    /////////////////
    // KNN example //
    /////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "KappaNearestNeighbor( objective function=RMSD, min kappa=1, max kappa=25)"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    /////////////////////////////
    // Kohonen network example //
    /////////////////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "Kohonen"
        "( "
        " objective function=EnrichmentAverage( cutoff = 4, enrichment max = 0.05, step size = 0.005, parity = 1), "
        " map dimensions(10.0,10.0),"
        " steps per update=0,"
        " length=10,"
        " radius=4,"
        " neighbor kernel = Gaussian"
        ")"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    //////////////////////////////////////////
    // KNN - kappa nearest neighbor example //
    //////////////////////////////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "KappaNearestNeighbor"
        "( "
        "  objective function=EnrichmentAverage( cutoff = 4, enrichment max = 0.05, step size = 0.005, parity = 1), "
        "  min kappa=1,"
        "  max kappa=5"
        ")"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    ///////////////////////////
    // Decision tree example //
    ///////////////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "DecisionTree( objective function=RMSD, partitioner=InformationGain, activity cutoff = 4)"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    //////////////////////////////////////////////
    // Decision tree example using storage flag //
    //////////////////////////////////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "DecisionTree( objective function=RMSD, partitioner=InformationGain, activity cutoff = 4)"
      );

      SetTrainingFlagsBinFiles( train_model, actives_bin, inactives_bin, filename_descriptor_labels);

      train_model.SetFlag
      (
        "storage_model",
        "File(directory=\"" + AddExampleOutputPathToFilename( "bcl::apps", "") + "\",prefix=\"\",write_descriptors=1)"
      );

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    ////////////////////////////////////////////////////////////////////////////////////////////
    // Neural network example using a dat file for training monitoring and independet dataset //
    ////////////////////////////////////////////////////////////////////////////////////////////

      train_model.ResetFlagsAndParameters();

      // specify machine learning technique
      train_model.AddParameter
      (
        "NeuralNetwork"
        "("
        " transfer function = Sigmoid,"
        " weight update = Simple(eta=0.1,alpha=0.5),"
        " objective function = RMSD,"
        " steps per update=1,"
        " hidden architecture(8)"
        ")"
      );

      SetTrainingFlagsBinFilesDatFile( train_model, data_points_dat);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( train_model.RunCommand(), 0);

    /////////////////
    // operations  //
    /////////////////

      // check that flags are needed
      BCL_ExampleCheck( train_model.CheckCommandString( true), true);

      // clean up example storage directory
      io::Directory dir( AddExampleOutputPathToFilename( "bcl::apps", ""));
      dir.Remove( true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppModelTrain

  const ExampleClass::EnumType ExampleAppModelTrain::s_Instance
  (
    GetExamples().AddEnum( ExampleAppModelTrain())
  );

} // namespace bcl
