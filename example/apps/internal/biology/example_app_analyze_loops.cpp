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
#include "internal/biology/bcl_app_analyze_loops.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "fold/bcl_fold.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_analyze_loops.cpp
  //! @brief this example tests the implementation of the application AnalyzeLoops, which analyzes the geometry of loop
  //! regions in protein models.
  //!
  //! @author fischea
  //! @date Dec 29, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppAnalyzeLoops :
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
    //! @return pointer to a new ExampleAppAnalyzeLoops
    ExampleAppAnalyzeLoops *Clone() const
    {
      return new ExampleAppAnalyzeLoops( *this);
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
      const app::ApplicationType app_loop( "AnalyzeLoops");
      BCL_ExampleAssert( app_loop.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // application helper for calling the object with command line arguments
      ApplicationExampleHelper loop_helper( app_loop);

      // create different test cases
      const std::string input_pdb_1_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      const std::string output_prefix
      (
        AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), "app_analyze_loop_test_")
      );

      // check test case 1
      loop_helper.SetFlag( "input_pdb", input_pdb_1_filename);
      loop_helper.SetFlag( "output_prefix", output_prefix);
      BCL_ExampleAssert( loop_helper.CheckCommandString( true), true);
      BCL_ExampleAssert( loop_helper.RunCommand(), 0);

      // // get paths and file names
      // const std::string input_pdb_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      // const std::string loop_template_library( AddExampleInputPathToFilename( e_Fold, "loop_library.ls"));
      // const std::string output_prefix
      // (
      //  AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), "app_loop_test_")
      // );
      // const std::string number_models( "3");

      // // test building five models for the test PDB
      // loop_helper.SetFlag( "input_pdb", input_pdb_filename);
      // loop_helper.SetFlag( "loop_library", loop_template_library);
      // loop_helper.SetFlag( "output_prefix", output_prefix);
      // loop_helper.SetFlag( "num_models", number_models);

      // // check the command line and run the application
      // BCL_ExampleAssert( loop_helper.CheckCommandString( true), true);
      // BCL_ExampleAssert( loop_helper.RunCommand(), 0);

      // // TODO check correctness of the created models

      // // remove output files
      // remove( ( output_prefix + "0.pdb").c_str());
      // remove( ( output_prefix + "1.pdb").c_str());
      // remove( ( output_prefix + "2.pdb").c_str());

      return 0;
    }

  }; // class ExampleAppAnalyzeLoops

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppAnalyzeLoops::s_Instance
  (
     GetExamples().AddEnum( ExampleAppAnalyzeLoops())
  );

} // namespace bcl
