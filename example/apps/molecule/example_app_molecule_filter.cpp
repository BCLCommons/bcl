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
#include "molecule/bcl_app_molecule_filter.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_feed_from_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_molecule_filter.cpp
  //!
  //! @author mendenjl
  //! @date May 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppMoleculeFilter :
    public ExampleInterface
  {
  public:

    ExampleAppMoleculeFilter *Clone() const
    {
      return new ExampleAppMoleculeFilter( *this);
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

      ApplicationExampleHelper molecule_filter_helper( app::MoleculeFilter::MoleculeFilter_Instance);

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
        molecule_filter_helper.GetThisApplicationExampleInputPath() + "csd_whole" + extensions
      );
      const std::string output_filename_good_atm_types_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_good_atm_types" + extensions)
      );

      // determine the number of molecules in the complete csd
      const size_t total_csd_size( chemistry::FragmentFeedFromFile::CountMolecules( csd_input_filename));

      const std::string output_filename_bad_atm_types_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_bad_atm_types" + extensions)
      );

      // create a command line that will generate descriptors from a stored sdf
      molecule_filter_helper.ResetFlagsAndParameters();

      molecule_filter_helper.SetFlag( "input_filenames", csd_input_filename);
      molecule_filter_helper.SetFlag( "output_matched", output_filename_good_atm_types_test);
      molecule_filter_helper.SetFlag( "output_unmatched", output_filename_bad_atm_types_test);
      molecule_filter_helper.SetFlag( "defined_atom_types");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_filter_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_filter_helper.RunCommand(), 0))
      {
        // get the file sizes
        const size_t number_mols_good_atom_types
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_good_atm_types_test)
        );
        const size_t number_mols_bad_atom_types
        (
          chemistry::FragmentFeedFromFile::CountMolecules( output_filename_bad_atm_types_test)
        );

        bool all_checks_passed( true);

        // check that all molecules were stored in one of the two files
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            number_mols_good_atom_types + number_mols_bad_atom_types,
            total_csd_size,
            "All molecules were put in one of the two files"
          );

        // check the count of molecules
        all_checks_passed &=
          BCL_ExampleIndirectCheckWithinAbsTolerance
          (
            number_mols_good_atom_types,
            54320,
            2,
            "number molecules with good atom types"
          );

        chemistry::FragmentFeed feed_good
        (
          storage::Vector< std::string>( 1, output_filename_good_atm_types_test),
          sdf::e_Maintain
        );
        chemistry::FragmentFeed feed_bad
        (
          storage::Vector< std::string>( 1, output_filename_bad_atm_types_test),
          sdf::e_Maintain
        );
        // check that all the molecules in the bad ensemble have bad atom types
        bool had_undefined( false);
        for( ; feed_bad.NotAtEnd(); ++feed_bad)
        {
          had_undefined = false;
          for
          (
            iterate::Generic< const chemistry::AtomConformationalInterface> itr( feed_bad->GetAtomsIterator());
            itr.NotAtEnd();
            ++itr
          )
          {
            if( !itr->GetAtomType()->IsGasteigerAtomType())
            {
              had_undefined = true;
              break;
            }
          }
          if( !had_undefined)
          {
            break;
          }
        }
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            had_undefined,
            true,
            "All molecules in the unmatched list should have had completely valid atom types"
          );
        // check that all the molecules in the good ensemble have good atom types
        had_undefined = false;
        for( ; feed_good.NotAtEnd() && !had_undefined; ++feed_good)
        {
          for
          (
            iterate::Generic< const chemistry::AtomConformationalInterface> itr( feed_good->GetAtomsIterator());
            itr.NotAtEnd();
            ++itr
          )
          {
            if( !itr->GetAtomType()->IsGasteigerAtomType())
            {
              had_undefined = true;
              break;
            }
          }
        }
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            had_undefined,
            false,
            "All molecules in the matched list should have had completely valid atom types"
          );
        if( all_checks_passed)
        {
          remove( output_filename_good_atm_types_test.c_str());
          remove( output_filename_bad_atm_types_test.c_str());
        }
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

  const ExampleClass::EnumType ExampleAppMoleculeFilter::s_Instance
  (
    GetExamples().AddEnum( ExampleAppMoleculeFilter())
  );

} // namespace bcl
