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
#include "restraint/bcl_app_restraint_saxs_prep.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_format.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_restraint_saxs_prep.cpp
  //!
  //! @author putnamdk
  //! @date August 26, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppRestraintSaxsPrep :
    public ExampleInterface
  {
  public:

    ExampleAppRestraintSaxsPrep *Clone() const
    {
      return new ExampleAppRestraintSaxsPrep( *this);
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

      // get application
      const app::ApplicationType app_enum_restraint_saxs_prep( "restraint:SaxsPrep");
      BCL_ExampleAssert( app_enum_restraint_saxs_prep.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // Set Global variables

      // Read in experimental SAXS data
      const std::string saxs_file_name( AddExampleInputPathToFilename( e_Biology, "3ICL_SAXS.gnom"));
      const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "3ICL.pdb"));
      std::string maximum_dimension( "65.46");
      std::string sampling_rounds( "10");
      ApplicationExampleHelper restraint_saxs_prep_helper( app_enum_restraint_saxs_prep);

      // Test the simulate error function of the application
      {
        BCL_MessageStd( " Inside Error Simulation Test: ");

        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "3ICL_simulated_errors_correct.data")
        );

        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "3ICL_simulated_errors.data")
        );

        // application parameters and flags
        restraint_saxs_prep_helper.SetFlag( "exp_data", saxs_file_name);
        restraint_saxs_prep_helper.SetFlag( "simulate_errors");
        restraint_saxs_prep_helper.SetFlag( "output_file", data_file_name);

        // check the command line
        BCL_ExampleAssert( restraint_saxs_prep_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_prep_helper.RunCommand(), 0);

        // stream to read data files
        io::IFStream read;

        // read output file created by this example
        BCL_ExampleMustOpenInputFile( read, data_file_name);
        restraint::SasScatteringData scattering_data;
        scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        // read known correct data file
        BCL_ExampleMustOpenInputFile( read, correct_data_file_name);
        restraint::SasScatteringData correct_scattering_data;
        correct_scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        bool match
        (
          io::File::ObjectsMatchWithinAbsoluteTolerance
          (
            scattering_data.GetScatteringData(),
            correct_scattering_data.GetScatteringData(),
            0.001
          )
        );

        // check if data is the same
        BCL_ExampleIndirectCheck( match, true, "Error Simulation");
      }

      restraint_saxs_prep_helper.ResetFlagsAndParameters();

      // Test the reduce data function with min error as a selection criteria
      {
        BCL_MessageDbg( " Inside Reduce Data with min error as selection critera: ");

        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "3ICL_reduced_min_error_correct.data")
        );

        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "3ICL_reduced_min_error.data")
        );

        // application parameters and flags
        restraint_saxs_prep_helper.SetFlag( "exp_data", saxs_file_name);
        restraint_saxs_prep_helper.SetFlag( "reduce_data_min_error");
        restraint_saxs_prep_helper.SetFlag( "dmax", maximum_dimension);
        restraint_saxs_prep_helper.SetFlag( "output_file", data_file_name);

        // check the command line
        BCL_ExampleAssert( restraint_saxs_prep_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_prep_helper.RunCommand(), 0);

        // stream to read data files
        io::IFStream read;

        // read output file created by this example
        BCL_ExampleMustOpenInputFile( read, data_file_name);
        restraint::SasScatteringData scattering_data;
        scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        // read known correct data file
        BCL_ExampleMustOpenInputFile( read, correct_data_file_name);
        restraint::SasScatteringData correct_scattering_data;
        correct_scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        bool match
        (
          io::File::ObjectsMatchWithinAbsoluteTolerance
          (
            scattering_data.GetScatteringData(),
            correct_scattering_data.GetScatteringData(),
            0.001
          )
        );

        // check if data is the same
        BCL_ExampleIndirectCheck( match, true, "Data reduction with min error approximation");
      }

      restraint_saxs_prep_helper.ResetFlagsAndParameters();

      // Test the reduce data function with full shannon sampling and noisy data reduction criteria
      {
        BCL_MessageDbg( " Inside Reduce Data with shannon sampling and noisy data reduction criteria: ");

        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "3ICL_reduced_shannon_correct.data")
        );

        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "3ICL_reduced_shannon.data")
        );

        // application parameters and flags
        restraint_saxs_prep_helper.SetFlag( "exp_data", saxs_file_name);
        restraint_saxs_prep_helper.SetFlag( "pdb_file", pdb_file_name);
        restraint_saxs_prep_helper.SetFlag( "reduce_data");
        restraint_saxs_prep_helper.SetFlag( "dmax", maximum_dimension);
        restraint_saxs_prep_helper.SetFlag( "sampling_rounds", sampling_rounds);
        restraint_saxs_prep_helper.SetFlag( "output_file", data_file_name);

        // check the command line
        BCL_ExampleAssert( restraint_saxs_prep_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_prep_helper.RunCommand(), 0);

        // stream to read data files
        io::IFStream read;

        // read output file created by this example
        BCL_ExampleMustOpenInputFile( read, data_file_name);
        restraint::SasScatteringData scattering_data;
        scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        // read known correct data file
        BCL_ExampleMustOpenInputFile( read, correct_data_file_name);
        restraint::SasScatteringData correct_scattering_data;
        correct_scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        bool match
        (
          io::File::ObjectsMatchWithinAbsoluteTolerance
          (
            scattering_data.GetScatteringData(),
            correct_scattering_data.GetScatteringData(),
            0.001
          )
        );

        // check if data is the same
        BCL_ExampleIndirectCheck( match, true, "Data reduction with Shannon Sampling");
      }

      restraint_saxs_prep_helper.ResetFlagsAndParameters();
      // Test the Saxs Fit Read function with setting the error values to a constant 1.0
      {
        BCL_MessageDbg( " Inside Saxs Fit and Error Setting check: ");
        std::string constant_error( "1.0");

        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "3ICL_saxs_fit_constant_error_correct.data")
        );

        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "3ICL_saxs_fit_constant_error.data")
        );

        // application parameters and flags
        restraint_saxs_prep_helper.SetFlag( "exp_data", saxs_file_name);
        restraint_saxs_prep_helper.SetFlag( "gnome_fit");
        restraint_saxs_prep_helper.SetFlag( "set_error", constant_error);
        restraint_saxs_prep_helper.SetFlag( "output_file", data_file_name);

        // check the command line
        BCL_ExampleAssert( restraint_saxs_prep_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_prep_helper.RunCommand(), 0);

        // stream to read data files
        io::IFStream read;

        // read output file created by this example
        BCL_ExampleMustOpenInputFile( read, data_file_name);
        restraint::SasScatteringData scattering_data;
        scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        // read known correct data file
        BCL_ExampleMustOpenInputFile( read, correct_data_file_name);
        restraint::SasScatteringData correct_scattering_data;
        correct_scattering_data.ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        bool match
        (
          io::File::ObjectsMatchWithinAbsoluteTolerance
          (
            scattering_data.GetScatteringData(),
            correct_scattering_data.GetScatteringData(),
            0.001
          )
        );

        // check if data is the same
        BCL_ExampleIndirectCheck( match, true, "Error Set to constant value with Gnome Fit data");
      }
      restraint_saxs_prep_helper.ResetFlagsAndParameters();

      // Test the Loop simulation for regula falsi approximation
      {
        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "3ICL_simulated_analytic_loops_correct.pdb")
        );

        const std::string data_file_name
        (
           AddExampleOutputPathToFilename( assemble::ProteinModel(), "3ICL_simulated_analytic_loops.pdb")
        );

        // application parameters and flags
        restraint_saxs_prep_helper.SetFlag( "pdb_file", pdb_file_name);
        restraint_saxs_prep_helper.SetFlag( "use_analytic");
        restraint_saxs_prep_helper.SetFlag( "output_model", data_file_name);

        // check the command line
        BCL_ExampleAssert( restraint_saxs_prep_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_prep_helper.RunCommand(), 0);

        // read output file created by this example
        assemble::ProteinModel created_protein_model
        (
          pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename( data_file_name)
        );

        // read known correct data file
        assemble::ProteinModel correct_protein_model
        (
          pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename( correct_data_file_name)
        );

        bool match
        (
          io::File::ObjectsMatchWithinAbsoluteTolerance
          (
            created_protein_model.GetAtomCoordinates(),
            correct_protein_model.GetAtomCoordinates(),
            0.15
          )
        );

        BCL_ExampleIndirectCheck( match, true, "Analytic Normalization loop approximation");
      }

      restraint_saxs_prep_helper.ResetFlagsAndParameters();

      // Test the Loop simulation for triangular approximation
      {
        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "3ICL_simulated_triangular_loops_correct.pdb")
        );

        const std::string data_file_name
        (
           AddExampleOutputPathToFilename( assemble::ProteinModel(), "3ICL_simulated_triangular_loops.pdb")
        );

        // application parameters and flags
        restraint_saxs_prep_helper.SetFlag( "pdb_file", pdb_file_name);
        restraint_saxs_prep_helper.SetFlag( "output_model", data_file_name);

        // check the command line
        BCL_ExampleAssert( restraint_saxs_prep_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_prep_helper.RunCommand(), 0);

        // read output file created by this example
        assemble::ProteinModel created_protein_model
        (
          pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename( data_file_name)
        );

        // read known correct data file
        assemble::ProteinModel correct_protein_model
        (
          pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename( correct_data_file_name)
        );

        bool match
        (
          io::File::ObjectsMatchWithinAbsoluteTolerance
          (
            created_protein_model.GetAtomCoordinates(),
            correct_protein_model.GetAtomCoordinates(),
            0.15
          )
        );

        BCL_ExampleIndirectCheck( match, true, "Fast triangular loop approximation");
      }

      // reset all pdb factory flags, since this application changes them
      restraint_saxs_prep_helper.ResetFlagsAndParameters();
      pdb::Factory::ResetFlagDefaults();

      // end
      return 0;
    }

    // put string comparison check in io file

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppRestraintSaxsPrep

  const ExampleClass::EnumType ExampleAppRestraintSaxsPrep::s_Instance
  (
    GetExamples().AddEnum( ExampleAppRestraintSaxsPrep())
  );

} // namespace bcl
