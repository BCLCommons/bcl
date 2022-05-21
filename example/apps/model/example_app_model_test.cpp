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
#include "model/bcl_app_model_test.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
//#include "io/bcl_io_directory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_model_test.cpp
  //! @brief application example for application model::Test
  //!
  //! @author butkiem1
  //! @date Jun 22, 2013
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppModelTest :
    public ExampleInterface
  {
  public:

    ExampleAppModelTest *Clone() const
    {
      return new ExampleAppModelTest( *this);
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
      ApplicationExampleHelper test_model( app::ModelTest::ModelTest_Instance);

    ///////////
    // files //
    ///////////

    ///////////////////////
    // integration tests //
    ///////////////////////

    /////////////////////////////
    // load in model from file //
    /////////////////////////////

      // reset all flags and parameters for that command line
      test_model.ResetFlagsAndParameters();

      std::string example_path( AddExampleInputPathToFilename( e_Model, ""));

      // set dataset retriever
      test_model.SetFlag
      (
        "retrieve_dataset",
        "Subset(filename=" + example_path + "dataset_aid891_100act_100inact.bin)"
      );

      // set model storage where the appropriate model is stored
      test_model.SetFlag
      (
        "storage_model",
        "File(directory=" + example_path + ",prefix=model,key=000011)"
      );

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( test_model.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      BCL_ExampleCheck( test_model.RunCommand(), 0);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppModelTest

  const ExampleClass::EnumType ExampleAppModelTest::s_Instance
  (
    GetExamples().AddEnum( ExampleAppModelTest())
  );

} // namespace bcl
