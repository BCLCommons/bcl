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
#include "internal/biology/bcl_app_loop_template.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_loop_template.cpp
  //! @brief this example tests the implementation of the application LoopTemplate, which creates a loop
  //! conformation library from a given set of pdbs
  //!
  //! @author fischea
  //! @date Dec 16, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppLoopTemplate :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleAppLoopTemplate
    ExampleAppLoopTemplate *Clone() const
    {
      return new ExampleAppLoopTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create an application instance for testing
      const app::ApplicationType app_loop_template( "LoopTemplate");
      BCL_ExampleAssert( app_loop_template.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // application helper for calling the object with command line arguments
      ApplicationExampleHelper loop_template_helper( app_loop_template);

      // get paths and filenames
      const std::string pdb_list( AddExampleInputPathToFilename( e_Fold, "loop_templates.ls"));
      const std::string storage_path( AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), ""));
      const std::string library_output_file
      (
        AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), "loop_library.ls")
      );
      const std::string statistics_output_file
      (
        AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), "loop_statistics.tbl")
      );

      // test using a single pdb file
      loop_template_helper.SetFlag( "pdb_list", pdb_list);
      loop_template_helper.SetFlag( "library_output_path", library_output_file);
      loop_template_helper.SetFlag( "statistics_output_path", statistics_output_file);
      const std::string correct_library_file( AddExampleInputPathToFilename( e_Fold, "loop_library.ls"));
      const std::string correct_statistics_file( AddExampleInputPathToFilename( e_Fold, "loop_statistics.tbl"));

      // check the command line and run the application
      BCL_ExampleAssert( loop_template_helper.CheckCommandString( true), true);
      BCL_ExampleAssert( loop_template_helper.RunCommand(), 0);

      // check for correctness of the created loop template library
      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( library_output_file, correct_library_file, 0.1),
          true,
          "Loop template library was generated incorrectly"
        )
      )
      {
        // if the example check was passed remove the output file
        remove( library_output_file.c_str());
      }

      // check for correctness of the created loop statistics
      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( statistics_output_file, correct_statistics_file, 0.1),
          true,
          "Loop template statistics were generated incorrectly"
        )
      )
      {
        // if the example check was passed remove the output file
        remove( statistics_output_file.c_str());
      }

      return 0;
    }

  }; // class ExampleAppLoopTemplate

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppLoopTemplate::s_Instance
  (
     GetExamples().AddEnum( ExampleAppLoopTemplate())
  );

} // namespace bcl
