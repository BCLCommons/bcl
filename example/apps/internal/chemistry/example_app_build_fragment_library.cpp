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
#include "internal/chemistry/bcl_app_build_fragment_library.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_app_build_fragment_library.cpp
  //! @brief this example tests the implementation of the application BuildFragmentLibrary, which fragments a given
  //! ensemble of molecules
  //!
  //! @see @link example_app_build_fragment_library.cpp @endlink
  //! @author kothiwsk
  //! @date Jan 28, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppBuildFragmentLibrary :
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
    ExampleAppBuildFragmentLibrary *Clone() const
    {
      return new ExampleAppBuildFragmentLibrary( *this);
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

      ApplicationExampleHelper build_fragment_library_helper( app::BuildFragmentLibrary::BuildFragmentLibrary_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      std::string extensions( ".sdf");
      // determine whether file compression can be used
      if( io::GetStreamBufferClasses().GetCompressionFromExtension( "bz2").IsDefined())
      {
        extensions += ".bz2";
      }

      const std::string input_ensemble( AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));

      const std::string output_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "mGluR5_app_fragments.sdf"));

      // create a command line that will generate descriptors from a stored sdf
      build_fragment_library_helper.ResetFlagsAndParameters();

      build_fragment_library_helper.SetFlag( "molecules_filenames", input_ensemble);
      build_fragment_library_helper.SetFlag( "output", output_filename);
      build_fragment_library_helper.SetFlag( "constitutions");
      build_fragment_library_helper.SetFlag( "max_rot", "1000");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( build_fragment_library_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( build_fragment_library_helper.RunCommand(), 0))
      {
        std::string correct_output( output_filename + ".correct");
        // check if the generated file is right
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_filename, correct_output),
            true,
            "AppBuildFragmentLibrary"
          )
        )
        {
          remove( output_filename.c_str());
        }

      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      }

  }; // class ExampleAppBuildFragmentLibrary

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppBuildFragmentLibrary::s_Instance
  (
     GetExamples().AddEnum( ExampleAppBuildFragmentLibrary())
  );

} // namespace bcl
