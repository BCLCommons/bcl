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
#include "model/bcl_app_model_compute_statistics.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_directory_entry.h"
#include "model/bcl_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_model_compute_statistics.cpp
  //! @brief application example for application app::ModelComputeStatistics
  //!
  //! @author butkiem1
  //! @date Feb 23, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppModelComputeStatistics :
    public ExampleInterface
  {
  public:

    ExampleAppModelComputeStatistics *Clone() const
    {
      return new ExampleAppModelComputeStatistics( *this);
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
      ApplicationExampleHelper comp_jury_stat
      (
        app::ModelComputeStatistics::ModelComputeStatistics_Instance
      );

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing command line without flags and parameters
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      comp_jury_stat.ResetFlagsAndParameters();

      // check an invalid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( comp_jury_stat.CheckCommandString( false), false);

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case with standard parameters and evaluate a roc curve
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // specify output directory
      std::string input_directory_name( AddExampleInputPathToFilename( ExampleInterface::e_Model, ""));

      storage::Vector< std::string> input_flag;
      input_flag.PushBack( input_directory_name + "prediction_merge_1.txt");
      input_flag.PushBack( input_directory_name + "prediction_merge_2.txt");

      // specify output directory
      const std::string output_directory_name( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "/compute_jury_statistics/"));

      // create output directory
      io::Directory output_directory( io::Directory::MkDir( output_directory_name));

      // check if directory exists
      BCL_ExampleCheck( output_directory.DoesExist(), true);

      // create a command line
      comp_jury_stat.ResetFlagsAndParameters();

      // flags
      comp_jury_stat.SetFlag( "input", input_flag( 0));
      comp_jury_stat.SetFlag( "potency_cutoff", "1.0");
      comp_jury_stat.SetFlag( "table_name", "compute_jury_statistics_table.txt");
      comp_jury_stat.SetFlag( "output_directory", output_directory_name);
      comp_jury_stat.SetFlag( "image_format", "png");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      if( BCL_ExampleCheck( comp_jury_stat.CheckCommandString( true), true))
      {
        // run a valid set of flags, check that the return status is 0
        if( BCL_ExampleCheck( comp_jury_stat.RunCommand(), 0))
        {
          BCL_MessageStd( "command line executed successfully!");
        }
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case where a correlation analysis is evaluated
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      comp_jury_stat.ResetFlagsAndParameters();

      // flags
      comp_jury_stat.SetFlag( "input", input_flag( 0));
      comp_jury_stat.SetFlag( "correlation");
      comp_jury_stat.SetFlag( "table_name", "compute_jury_statistics_table.txt");
      comp_jury_stat.SetFlag( "output_directory", output_directory_name);
      comp_jury_stat.SetFlag( "image_format", "png");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      if( BCL_ExampleCheck( comp_jury_stat.CheckCommandString( true), true))
      {
        // run a valid set of flags, check that the return status is 0
        if( BCL_ExampleCheck( comp_jury_stat.RunCommand(), 0))
        {
          BCL_MessageStd( "command line executed successfully!");
        }
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case multiple predictions are passed as input and a consensus analyis is performed
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      comp_jury_stat.ResetFlagsAndParameters();

      // flags
      comp_jury_stat.SetFlag( "input", input_flag);
      comp_jury_stat.SetFlag( "potency_cutoff", "1.0");
      comp_jury_stat.SetFlag( "table_name", "compute_jury_statistics_table.txt");
      comp_jury_stat.SetFlag( "output_directory", output_directory_name);
      comp_jury_stat.SetFlag( "image_format", "png");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      if( BCL_ExampleCheck( comp_jury_stat.CheckCommandString( true), true))
      {
        // run a valid set of flags, check that the return status is 0
        if( BCL_ExampleCheck( comp_jury_stat.RunCommand(), 0))
        {
          BCL_MessageStd( "command line executed successfully!");
        }
      }

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // cleanup
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // clean up generated output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).Remove( true), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppModelComputeJuryStatistics

  const ExampleClass::EnumType ExampleAppModelComputeStatistics::s_Instance
  (
    GetExamples().AddEnum( ExampleAppModelComputeStatistics())
  );

} // namespace bcl
