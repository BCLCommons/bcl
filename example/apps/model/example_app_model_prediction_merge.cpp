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
#include "model/bcl_app_model_prediction_merge.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_matrix.h"
#include "model/bcl_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_train_model_prediction_merge.cpp
  //! @brief application example for application app::TrainModelPredictionMerge.
  //!
  //! @author butkiem1
  //! @date Mar 11, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppModelPredictionMerge :
    public ExampleInterface
  {
  public:

    ExampleAppModelPredictionMerge *Clone() const
    {
      return new ExampleAppModelPredictionMerge( *this);
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
      ApplicationExampleHelper train_model_pred_merge
      (
        app::ModelPredictionMerge::ModelPredictionMerge_Instance
      );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing command line without flags and parameters
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      train_model_pred_merge.ResetFlagsAndParameters();

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model_pred_merge.CheckCommandString( false), false);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case where prediction matrices are appended
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // specify output directory
      const std::string input_directory_name( AddExampleInputPathToFilename( ExampleInterface::e_Model, ""));

      // specify output directory
      const std::string output_directory_name( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), ""));

      // create a command line
      train_model_pred_merge.ResetFlagsAndParameters();

      storage::Vector< std::string> input_flag;
      input_flag.PushBack( input_directory_name + "prediction_merge_1.txt");
      input_flag.PushBack( input_directory_name + "prediction_merge_2.txt");

      // flags
      train_model_pred_merge.SetFlag( "input", input_flag);

      train_model_pred_merge.SetFlag( "output", output_directory_name + "prediction_merge_result.txt");
      train_model_pred_merge.SetFlag( "modus");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model_pred_merge.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( train_model_pred_merge.RunCommand(), 0))
      {
        BCL_MessageStd( "command line executed successfully!");
      }

      // check whether result file was constructed
      io::IFStream result_file_stream;
      BCL_ExampleMustOpenInputFile( result_file_stream, output_directory_name + "prediction_merge_result.txt");

      // fill bcl object with content of file
      storage::Vector< storage::Vector< std::string> > result_matrix;
      result_matrix = util::SplittedStringLineListFromIStream( result_file_stream, "\n,");
      // close file stream
      io::File::CloseClearFStream( result_file_stream);

      // check the read in matrix for correctness
      BCL_ExampleCheck( result_matrix( 0).GetSize(), size_t( 2));
      BCL_ExampleCheck( result_matrix.GetSize(), size_t( 10));

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case where prediction matrices are merged
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      train_model_pred_merge.ResetFlagsAndParameters();

      input_flag.Reset();
      input_flag.PushBack( input_directory_name + "prediction_merge_1.txt");
      input_flag.PushBack( input_directory_name + "prediction_merge_2.txt");

      // flags
      train_model_pred_merge.SetFlag( "input", input_flag);

      train_model_pred_merge.SetFlag( "output", output_directory_name + "prediction_merge_result.txt");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model_pred_merge.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( train_model_pred_merge.RunCommand(), 0))
      {
        BCL_MessageStd( "command line executed successfully!");
      }

      // check whether result file was constructed
      BCL_ExampleMustOpenInputFile( result_file_stream, output_directory_name + "prediction_merge_result.txt");
      io::File::CloseClearFStream( result_file_stream);

      // check whether merge attempt result is the same matrix as matrix in comparison file
      BCL_ExampleCheck
      (
        io::File::FilesMatch
        (
          output_directory_name + "prediction_merge_result.txt",
          output_directory_name + "prediction_merge_result_compare.txt"
        ),
        true
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppModelPredictionMerge

  const ExampleClass::EnumType ExampleAppModelPredictionMerge::s_Instance
  (
    GetExamples().AddEnum( ExampleAppModelPredictionMerge())
  );

} // namespace bcl
