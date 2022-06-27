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
#include "example_application_example_helper.h"

// include the header of the class which this example is for
#include "molecule/bcl_app_molecule_reorder.h"

// bcl headers, sorted alphabetically
#include "chemistry/bcl_chemistry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_app_molecule_reorder.cpp
  //! @brief this example tests the implementation of the application MoleculeReorder, which changes the order of
  //! molecules (or atoms in those molecules)
  //!
  //! @see @link example_app_molecule_reorder.cpp @endlink
  //! @author geanesar
  //! @date Apr 18, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppMoleculeReorder :
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
    ExampleAppMoleculeReorder *Clone() const
    {
      return new ExampleAppMoleculeReorder( *this);
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

      ApplicationExampleHelper molecule_reorder_helper( app::MoleculeReorder::MoleculeReorder_Instance);

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

      const std::string output_reverse_basename( "mGluR5_actives_reversed.sdf");
      const std::string output_random_basename( "mGluR5_actives_randomized.sdf");
      const std::string output_atomorder_basename( "mGluR5_actives_atomorder.sdf");
      const std::string output_sort_basename( "mGluR5_actives_sort.sdf");

      const std::string output_reverse_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_reverse_basename));
      const std::string output_random_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_random_basename));
      const std::string output_atomorder_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_atomorder_basename));
      const std::string output_sort_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_sort_basename));

    ///////////////////
    // reverse order //
    ///////////////////

      // Set flags
      molecule_reorder_helper.ResetFlagsAndParameters();
      molecule_reorder_helper.SetFlag( "add_h");
      molecule_reorder_helper.SetFlag( "input_filenames", input_ensemble);
      molecule_reorder_helper.SetFlag( "output", output_reverse_filename);
      molecule_reorder_helper.SetFlag( "reverse");

      BCL_ExampleCheck( molecule_reorder_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( molecule_reorder_helper.RunCommand(), 0))
      {
        std::string correct_output( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_reverse_basename + ".correct"));

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_reverse_filename, correct_output),
            true,
            this->GetClassIdentifier() + " correct output reverse"
          )
        )
        {
          remove( output_reverse_filename.c_str());
        }
      }

    ////////////////////
    // sort molecules //
    ////////////////////

      // Set flags
      molecule_reorder_helper.ResetFlagsAndParameters();
      molecule_reorder_helper.SetFlag( "add_h");
      molecule_reorder_helper.SetFlag( "input_filenames", input_ensemble);
      molecule_reorder_helper.SetFlag( "output", output_sort_filename);
      molecule_reorder_helper.SetFlag( "sort", "NAtoms");

      BCL_ExampleCheck( molecule_reorder_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( molecule_reorder_helper.RunCommand(), 0))
      {
        std::string correct_output( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_sort_basename + ".correct"));

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_sort_filename, correct_output),
            true,
            this->GetClassIdentifier() + " correct output sort"
          )
        )
        {
          remove( output_sort_filename.c_str());
        }
      }

    ////////////////
    // randomized //
    ////////////////

      // Set flags
      molecule_reorder_helper.ResetFlagsAndParameters();
      molecule_reorder_helper.SetFlag( "add_h");
      molecule_reorder_helper.SetFlag( "input_filenames", input_ensemble);
      molecule_reorder_helper.SetFlag( "output", output_random_filename);
      molecule_reorder_helper.SetFlag( "randomize");

      BCL_ExampleCheck( molecule_reorder_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( molecule_reorder_helper.RunCommand(), 0))
      {
        std::string correct_output( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_random_basename + ".correct"));

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_random_filename, correct_output)
            || io::File::FilesMatch( output_random_filename, correct_output + ".mac"),
            true,
            this->GetClassIdentifier() + " correct output randomize"
          )
        )
        {
          remove( output_random_filename.c_str());
        }
      }

    ///////////////////////
    // atom order change //
    ///////////////////////

      // Set flags
      molecule_reorder_helper.ResetFlagsAndParameters();
      molecule_reorder_helper.SetFlag( "add_h");
      molecule_reorder_helper.SetFlag( "input_filenames", input_ensemble);
      molecule_reorder_helper.SetFlag( "output", output_atomorder_filename);

      molecule_reorder_helper.SetFlag( "atom_order", storage::Vector< std::string>::Create( "3", "0", "2", "6", "12", "15", "13", "5"));

      BCL_ExampleCheck( molecule_reorder_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( molecule_reorder_helper.RunCommand(), 0))
      {
        std::string correct_output( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_atomorder_basename + ".correct"));

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_atomorder_filename, correct_output),
            true,
            this->GetClassIdentifier() + " correct output atom_order"
          )
        )
        {
          remove( output_atomorder_filename.c_str());
        }
      }

      return 0;
    }

  }; // class ExampleAppMoleculeReorder

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppMoleculeReorder::s_Instance
  (
     GetExamples().AddEnum( ExampleAppMoleculeReorder())
  );

} // namespace bcl
