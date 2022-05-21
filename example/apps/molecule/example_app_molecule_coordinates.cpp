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
#include "molecule/bcl_app_molecule_coordinates.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_molecule_coordinates.cpp
  //!
  //! @author mendenjl
  //! @date May 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppMoleculeCoordinates :
    public ExampleInterface
  {
  public:

    ExampleAppMoleculeCoordinates *Clone() const
    {
      return new ExampleAppMoleculeCoordinates( *this);
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

      ApplicationExampleHelper molecule_coordinates_helper( app::MoleculeCoordinates::MoleculeCoordinates_Instance);

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
        molecule_coordinates_helper.GetApplicationExampleInputPath() + "/MoleculeFilter/csd_whole" + extensions
      );
      const std::string output_path
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "csd")
      );

      // create strings to hold all the output paths
      const std::string output_filename_bond_length_stats( output_path + ".bond_lengths.stats.txt");
      const std::string output_filename_bond_angle_stats( output_path + ".bond_angles.stats.txt");
      const std::string output_filename_bond_angle_histograms( output_path + ".bond_angles.histograms.txt");
      const std::string output_filename_dihedral_histograms( output_path + ".dihedral.histograms.txt");

      // create a command line that will generate descriptors from a stored sdf
      molecule_coordinates_helper.ResetFlagsAndParameters();

      molecule_coordinates_helper.SetFlag( "input_filenames", csd_input_filename);
      molecule_coordinates_helper.SetFlag( "output", output_path);
      molecule_coordinates_helper.SetFlag( "dihedral_bin_size", "15");
      molecule_coordinates_helper.SetFlag( "bond_angle_bin_size", "5");
      molecule_coordinates_helper.SetFlag( "statistics");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_coordinates_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_coordinates_helper.RunCommand(), 0))
      {
        bool all_checks_passed( true);

        // tolerance for histograms; 2 means that 2 items may be in one bin instead of another; use this to give some
        // wiggle room for numerical differences between 32/64 bit archs on such a large dataset
        const double histogram_tolerance( 2.0);

        // tolerance for angular statistics; differences < 0.1 degrees are irrelevant
        // no tolerance is needed for bond lengths
        const double angle_stats_tolerance( 0.1);

        // check all of the output files
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_filename_bond_length_stats, output_filename_bond_length_stats + ".correct"),
            true,
            "Bond length statistics"
          );
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance
            (
              output_filename_bond_angle_stats,
              output_filename_bond_angle_stats + ".correct",
              angle_stats_tolerance
            ),
            true,
            "Bond angle statistics"
          );
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance
            (
              output_filename_bond_angle_histograms,
              output_filename_bond_angle_histograms + ".correct",
              histogram_tolerance
            ),
            true,
            "Bond angle histograms"
          );
        all_checks_passed &=
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance
            (
              output_filename_dihedral_histograms,
              output_filename_dihedral_histograms + ".correct",
              histogram_tolerance
            ),
            true,
            "Dihedral angle histograms"
          );
        if( all_checks_passed)
        {
          remove( output_filename_bond_length_stats.c_str());
          remove( output_filename_bond_angle_stats.c_str());
          remove( output_filename_bond_angle_histograms.c_str());
          remove( output_filename_dihedral_histograms.c_str());
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

  const ExampleClass::EnumType ExampleAppMoleculeCoordinates::s_Instance
  (
    GetExamples().AddEnum( ExampleAppMoleculeCoordinates())
  );

} // namespace bcl
