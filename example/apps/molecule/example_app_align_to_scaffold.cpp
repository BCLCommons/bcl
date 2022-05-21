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
#include "molecule/bcl_app_align_to_scaffold.h"

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
  //! @example bcl_app_align_to_scaffold.cpp
  //! @brief this example tests the implementation of the application MoleculeAlignToScaffold, aligns emselbles of
  //! molecules to a given scaffold
  //!
  //! @see @link example_app_align_to_scaffold.cpp @endlink
  //! @author geanesar
  //! @date Apr 18, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppMoleculeAlignToScaffold :
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
    ExampleAppMoleculeAlignToScaffold *Clone() const
    {
      return new ExampleAppMoleculeAlignToScaffold( *this);
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

      ApplicationExampleHelper align_to_scaffold_helper( app::AlignToScaffold::AlignToScaffold_Instance);

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
      const std::string input_ensemble( AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf"));
      const std::string input_scaffold( AddExampleInputPathToFilename( e_Chemistry, "diazepam_base.sdf"));
      const std::string output_basename( "diazepam_aligned.sdf");
      const std::string output_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_basename));

    ////////////////////////////
    // test overall alignment //
    ////////////////////////////

      // Set Parameters
      align_to_scaffold_helper.ResetFlagsAndParameters();
      align_to_scaffold_helper.AddParameter( input_scaffold);
      align_to_scaffold_helper.AddParameter( input_ensemble);
      align_to_scaffold_helper.AddParameter( output_filename);

      BCL_ExampleCheck( align_to_scaffold_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( align_to_scaffold_helper.RunCommand(), 0))
      {
        std::string correct_output( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_basename + ".correct"));

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

    ////////////////////////////////
    // test alignment using atoms //
    ////////////////////////////////

      // Set Parameters
      align_to_scaffold_helper.ResetFlagsAndParameters();
      align_to_scaffold_helper.AddParameter( input_scaffold);
      align_to_scaffold_helper.AddParameter( input_ensemble);
      align_to_scaffold_helper.AddParameter( output_filename);
      align_to_scaffold_helper.SetFlag( "align_scaffold_atoms", storage::Vector< std::string>::Create( "0", "1", "2"));
      align_to_scaffold_helper.SetFlag( "align_ensemble_atoms", storage::Vector< std::string>::Create( "2", "8", "7"));

      BCL_ExampleCheck( align_to_scaffold_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( align_to_scaffold_helper.RunCommand(), 0))
      {
        // The output for this should be exactly the same as the previous one
        std::string correct_output( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_basename + ".correct"));

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_filename, correct_output),
            true,
            this->GetClassIdentifier() + " correct output for align_atoms"
          )
        )
        {
          remove( output_filename.c_str());
        }
      }

      return 0;
    }

  }; // class ExampleAppMoleculeAlignToScaffold

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppMoleculeAlignToScaffold::s_Instance
  (
     GetExamples().AddEnum( ExampleAppMoleculeAlignToScaffold())
  );

} // namespace bcl
