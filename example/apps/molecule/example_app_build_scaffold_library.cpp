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
#include "internal/chemistry/bcl_app_build_scaffold_library.h"

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
  //! @example bcl_app_build_scaffold_library.cpp
  //! @brief this example tests the implementation of the application BuildScaffoldLibrary, which generates the common
  //! substructures (i.e. scaffolds) of pairs of molecules
  //!
  //! @see @link example_app_build_scaffold_library.cpp @endlink
  //! @author geanesar
  //! @date Apr 18, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppBuildScaffoldLibrary :
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
    ExampleAppBuildScaffoldLibrary *Clone() const
    {
      return new ExampleAppBuildScaffoldLibrary( *this);
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

      ApplicationExampleHelper build_scaffold_library_helper( app::BuildScaffoldLibrary::BuildScaffoldLibrary_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // IO files
      const std::string input_ensemble( AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      const std::string output_basename( "mGluR5_scaffold_library.sdf");
      const std::string output_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_basename));

      // Set flags
      build_scaffold_library_helper.ResetFlagsAndParameters();
      build_scaffold_library_helper.SetFlag( "ensemble", input_ensemble);
      build_scaffold_library_helper.SetFlag( "output", output_filename);
      build_scaffold_library_helper.SetFlag( "sampling_fraction", "0.8");

      // Min size of 7 will get rid of benzene rings
      build_scaffold_library_helper.SetFlag( "min_size", "7");
      build_scaffold_library_helper.SetFlag( "ignore_scaffolds_with_open_rings", "1");
      build_scaffold_library_helper.SetFlag( "bond_coloring_scheme", "BondOrderOrAromaticWithRingness");

      BCL_ExampleCheck( build_scaffold_library_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( build_scaffold_library_helper.RunCommand(), 0))
      {
        std::string correct_output( output_filename + ".correct");

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_filename, correct_output),
            true,
            this->GetClassIdentifier() + " correct output"
          )
        )
        {
          remove( output_filename.c_str());
        }
      }

      return 0;
    }

  }; // class ExampleAppBuildScaffoldLibrary

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppBuildScaffoldLibrary::s_Instance
  (
     GetExamples().AddEnum( ExampleAppBuildScaffoldLibrary())
  );

} // namespace bcl
