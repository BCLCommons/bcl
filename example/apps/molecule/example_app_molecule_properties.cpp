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
#include "molecule/bcl_app_molecule_properties.h"

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
  //! @example example_app_molecule_properties.cpp
  //!
  //! @author mendenjl
  //! @date Apr 20, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppMoleculeProperties :
    public ExampleInterface
  {
  public:

    ExampleAppMoleculeProperties *Clone() const
    {
      return new ExampleAppMoleculeProperties( *this);
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

      ApplicationExampleHelper molecule_properties_helper( app::MoleculeProperties::MoleculeProperties_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      const std::string csd_input_filename
      (
        AddExampleInputPathToFilename( e_Chemistry, "csd_first_1115_simple.sdf")
      );
      const std::string output_filename_props_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_100mols_with_props.sdf")
      );

      const std::string output_filename_hists_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_100mols.hist")
      );
      const std::string output_filename_table_test
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd_100mols.table")
      );

      // create a command line that will generate descriptors from a stored sdf
      molecule_properties_helper.ResetFlagsAndParameters();

      molecule_properties_helper.SetFlag( "input_filenames", csd_input_filename);
      molecule_properties_helper.SetFlag( "output", output_filename_props_test);
      molecule_properties_helper.SetFlag( "output_histogram", output_filename_hists_test);
      molecule_properties_helper.SetFlag( "output_table", output_filename_table_test);
      molecule_properties_helper.SetFlag( "input_start", "5");
      molecule_properties_helper.SetFlag( "input_max", "100");
      molecule_properties_helper.SetFlag
      (
        "numeric_histogram",
        storage::Vector< std::string>::Create( "Atom_EffectivePolarizability", "1", "0.25", "30")
      );
      molecule_properties_helper.SetFlag
      (
        "string_histogram",
        storage::Vector< std::string>::Create( "AtomTypes", "BondTypes")
      );
      molecule_properties_helper.SetFlag
      (
        "tabulate",
        storage::Vector< std::string>::Create( "Numeric(NAtoms)", "Numeric(NStereo)", "AtomTypes")
      );
      molecule_properties_helper.SetFlag
      (
        "statistics",
        storage::Vector< std::string>::Create
        (
          "Atom_Vcharge",
          "Atom_SigmaCharge",
          "Atom_PiCharge",
          "Atom_EffectivePolarizability",
          "TotalFormalCharge",
          "TopologicalPolarSurfaceArea",
          "Equal(Atom_Stereocenters,Constant(1))",  // R
          "Equal(Atom_Stereocenters,Constant(-1))", // S
          "Equal(Atom_Stereocenters,Constant(0))"   // achiral
        )
      );
      molecule_properties_helper.SetFlag( "add", "Add(Atom_SigmaCharge,Atom_Vcharge)");
      molecule_properties_helper.SetFlag
      (
        "rename",
        storage::Vector< std::string>::Create( "Add(Atom_SigmaCharge,Atom_Vcharge)", "AtomSigmaVCharge")
      );

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_properties_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_properties_helper.RunCommand(), 0))
      {
        // check that the individual files are correct
        if
        (
          BCL_ExampleCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance
            (
              output_filename_props_test,
              output_filename_props_test + ".correct",
              0.05
            ),
            true
          )
        )
        {
          remove( output_filename_props_test.c_str());
        }
        if
        (
          BCL_ExampleCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance
            (
              output_filename_hists_test,
              output_filename_hists_test + ".correct",
              0.4
            ),
            true
          )
        )
        {
          remove( output_filename_hists_test.c_str());
        }
        if
        (
          BCL_ExampleCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance
            (
              output_filename_table_test,
              output_filename_table_test + ".correct",
              0.05
            ),
            true
          )
        )
        {
          remove( output_filename_table_test.c_str());
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

  const ExampleClass::EnumType ExampleAppMoleculeProperties::s_Instance
  (
    GetExamples().AddEnum( ExampleAppMoleculeProperties())
  );

} // namespace bcl
