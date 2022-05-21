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
#include "molecule/bcl_app_molecule_unique.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry_fragment_feed_from_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_molecule_unique.cpp
  //!
  //! @author butkiem1
  //! @date Nov 28, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppMoleculeUnique :
    public ExampleInterface
  {
  public:

    ExampleAppMoleculeUnique *Clone() const
    {
      return new ExampleAppMoleculeUnique( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      ApplicationExampleHelper molecule_unique_helper( app::MoleculeUnique::MoleculesUnique_Instance);

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

      const std::string csd_input_filename
      (
        molecule_unique_helper.GetThisApplicationExampleInputPath() + "csd_1000mols" + extensions
      );
      const std::string output_filename_unique_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_1000mols_unique" + extensions)
      );

      // determine the number of molecules in the complete csd
      const size_t total_csd_size( chemistry::FragmentFeedFromFile::CountMolecules( csd_input_filename));

      const std::string output_filename_duplicates_mols_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_1000mols_duplicate" + extensions)
      );

      // create a command line that will generate descriptors from a stored sdf
      molecule_unique_helper.ResetFlagsAndParameters();

      molecule_unique_helper.SetFlag( "input_filenames", csd_input_filename);
      molecule_unique_helper.SetFlag( "output", output_filename_unique_test);
      molecule_unique_helper.SetFlag( "output_dupes", output_filename_duplicates_mols_test);
      molecule_unique_helper.SetFlag( "compare", "Constitutions");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_unique_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_unique_helper.RunCommand(), 0))
      {
        // get the file sizes
        const size_t number_mols_unique_test
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_unique_test)
        );
        const size_t number_mols_duplicates_test
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_duplicates_mols_test)
        );

        bool all_checks_passed( true);

        // check that all molecules were stored in one of the two files
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_unique_test + number_mols_duplicates_test,
            total_csd_size,
            "All molecules were put in one of the two files"
          );

        // check the count of molecules
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_unique_test,
            866,
            "number molecules with good atom types"
          );
      }

      // create a command line that will generate descriptors from a stored sdf
      molecule_unique_helper.ResetFlagsAndParameters();

      molecule_unique_helper.SetFlag( "input_filenames", csd_input_filename);
      molecule_unique_helper.SetFlag( "output", output_filename_unique_test);
      molecule_unique_helper.SetFlag( "output_dupes", output_filename_duplicates_mols_test);
      molecule_unique_helper.SetFlag( "compare", "Exact");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_unique_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_unique_helper.RunCommand(), 0))
      {
        // get the file sizes
        const size_t number_mols_unique_test_two
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_unique_test)
        );
        const size_t number_mols_duplicates_test_two
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_duplicates_mols_test)
        );

        bool all_checks_passed( true);

        // check that all molecules were stored in one of the two files
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_unique_test_two + number_mols_duplicates_test_two,
            total_csd_size,
            "All molecules were put in one of the two files"
          );

        // check the count of molecules
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_unique_test_two,
            914,
            "number molecules with good atom types"
          );
      }

      // create a command line that will generate descriptors from a stored sdf
      molecule_unique_helper.ResetFlagsAndParameters();

      molecule_unique_helper.SetFlag( "input_filenames", csd_input_filename);
      molecule_unique_helper.SetFlag( "output", output_filename_unique_test);
      molecule_unique_helper.SetFlag( "output_dupes", output_filename_duplicates_mols_test);
      molecule_unique_helper.SetFlag( "compare", "Exact");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_unique_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_unique_helper.RunCommand(), 0))
      {
        // get the file sizes
        const size_t number_mols_unique_test_two
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_unique_test)
        );
        const size_t number_mols_duplicates_test_two
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_duplicates_mols_test)
        );

        bool all_checks_passed( true);

        // check that all molecules were stored in one of the two files
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_unique_test_two + number_mols_duplicates_test_two,
            total_csd_size,
            "All molecules were put in one of the two files"
          );

        // check the count of molecules
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_unique_test_two,
            914,
            "number molecules with good atom types"
          );
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end Exampleexample_name

  const ExampleClass::EnumType ExampleAppMoleculeUnique::s_Instance
  (
    GetExamples().AddEnum( ExampleAppMoleculeUnique())
  );

} // namespace bcl
