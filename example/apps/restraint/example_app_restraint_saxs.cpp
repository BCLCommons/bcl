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
#include "restraint/bcl_app_restraint_saxs.h"

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
  //! @example example_app_restraint_saxs.cpp
  //!
  //! @author putnamdk, heinzes1
  //! @date March 16, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppRestraintSaxs :
    public ExampleInterface
  {
  public:

    ExampleAppRestraintSaxs *Clone() const
    {
      return new ExampleAppRestraintSaxs( *this);
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
      const app::ApplicationType app_enum_restraint_saxs( "restraint:Saxs");
      BCL_ExampleAssert( app_enum_restraint_saxs.IsDefined(), true);

      //  Files to be used in the analysis
      const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ENH.pdb"));
      const std::string saxs_file_name( AddExampleInputPathToFilename( e_Biology, "1ENH00.saxs"));
      const std::string sasa_file_name( AddExampleInputPathToFilename( e_Biology, "1ENH.area"));

    ////////////////
    // operations //
    ////////////////

      // restraint:Saxs 'SasDebye(consider loops=0, analytic=0)'
      // -pdb 1ENH.pdb
      // -exp_data 1ENH00.saxs
      // -aaclass AAComplete
      // -rmsd
      {
        BCL_MessageStd( " Inside Application Example: ");

        ApplicationExampleHelper restraint_saxs_helper( app_enum_restraint_saxs);
        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "1ENH_aacomplete_raw_saxs_correct.data")
        );
        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "1ENH_aacomplete_raw_saxs.data")
        );

        // application parameters and flags
        restraint_saxs_helper.SetFlag( "pdb", pdb_file_name);
        restraint_saxs_helper.SetFlag( "exp_data", saxs_file_name);
        //restraint_saxs_helper.SetFlag( "sasa_data", sasa_file_name);
        restraint_saxs_helper.SetFlag( "aaclass", "AAComplete");
        restraint_saxs_helper.SetFlag( "output_file", data_file_name);
        restraint_saxs_helper.SetFlag( "use_errors");

        // check the command line
        BCL_ExampleAssert( restraint_saxs_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_helper.RunCommand(), 0);

        // stream to read data files
        io::IFStream read;

        // read output file created by this example
        BCL_ExampleMustOpenInputFile( read, data_file_name);
        restraint::SasScatteringData scattering_data;
        scattering_data.ReadFromDataFile( read);
// scattering_data.ReadBCLModel( read);
        io::File::CloseClearFStream( read);

        // read known correct data file
        BCL_ExampleMustOpenInputFile( read, correct_data_file_name);
        restraint::SasScatteringData correct_scattering_data;
        correct_scattering_data.ReadFromDataFile( read);
// correct_scattering_data.ReadBCLModel( read);
        io::File::CloseClearFStream( read);

        // Turn both objects into strings
        std::string created_object( util::Format()( scattering_data.GetScatteringData()));
        std::string correct_object( util::Format()( correct_scattering_data.GetScatteringData()));

        // Construct standard string streams
        std::stringstream created_stream( created_object);
        std::stringstream correct_stream( correct_object);

        // Test stream within tolerance
        bool match( io::File::StreamsMatchWithinAbsoluteTolerance( created_stream, correct_stream, 0.001));

        // check if data is the same
        BCL_ExampleCheck( match, true);
      }

      // restraint:Saxs 'SasDebye(consider loops=0, analytic=0)'
      // -pdb 1ENH.pdb
      // -exp_data 1ENH00.saxs
      // -aaclass AABackBone
      // -rmsd
      {
        ApplicationExampleHelper restraint_saxs_helper( app_enum_restraint_saxs);
        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "1ENH_aabackbone_raw_saxs_correct.data")
        );
        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "1ENH_aabackbone_raw_saxs.data")
        );

        // application parameters and flags
        restraint_saxs_helper.SetFlag( "pdb", pdb_file_name);
        restraint_saxs_helper.SetFlag( "exp_data", saxs_file_name);
        restraint_saxs_helper.SetFlag( "aaclass", "AABackBone");
        restraint_saxs_helper.SetFlag( "output_file", data_file_name);
        restraint_saxs_helper.SetFlag( "apx_sc");

        // check the command line
        BCL_ExampleAssert( restraint_saxs_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_helper.RunCommand(), 0);

        // stream to read data files
        io::IFStream read;

        // read output file created by this example
        BCL_ExampleMustOpenInputFile( read, data_file_name);
        restraint::SasScatteringData scattering_data;
        scattering_data.ReadFromDataFile( read);
// scattering_data.ReadBCLModel( read);
        io::File::CloseClearFStream( read);

        // read known correct data file
        BCL_ExampleMustOpenInputFile( read, correct_data_file_name);
        restraint::SasScatteringData correct_scattering_data;
        correct_scattering_data.ReadFromDataFile( read);
// correct_scattering_data.ReadBCLModel( read);
        io::File::CloseClearFStream( read);

        // check if data is the same
        BCL_ExampleCheck( scattering_data.GetScatteringData(), correct_scattering_data.GetScatteringData());
      }

      // restraint:Saxs 'SasDebye(consider loops=1, analytic=0)'
      // -pdb 1ENH.pdb
      // -exp_data 1ENH00.saxs
      // -aaclass AABackBone
      // -rmsd
      // -min_sse_size 5 3 999
      {
        ApplicationExampleHelper restraint_saxs_helper( app_enum_restraint_saxs);
        // get input and output file names
        const std::string correct_data_file_name
        (
          AddExampleInputPathToFilename( e_Biology, "1ENH_aabackbone_minssesize_raw_saxs_correct.data")
        );
        const std::string data_file_name
        (
          AddExampleOutputPathToFilename( restraint::SasScatteringData(), "1ENH_aabackbone_minssesize_raw_saxs.data")
        );

        // application parameters and flags
        restraint_saxs_helper.SetFlag( "pdb", pdb_file_name);
        restraint_saxs_helper.SetFlag( "exp_data", saxs_file_name);
        restraint_saxs_helper.SetFlag( "aaclass", "AABackBone");
        restraint_saxs_helper.SetFlag( "min_sse_size", storage::Vector< std::string>::Create( "5", "3", "999"));
        restraint_saxs_helper.SetFlag( "output_file", data_file_name);
        restraint_saxs_helper.SetFlag( "apx_sc");
        restraint_saxs_helper.SetFlag( "apx_loops");

        // check the command line
        BCL_ExampleAssert( restraint_saxs_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( restraint_saxs_helper.RunCommand(), 0);

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

        // check if data is the same
        BCL_ExampleCheck( scattering_data.GetScatteringData(), correct_scattering_data.GetScatteringData());
      }

      // reset all pdb factory flags, since this application changes them
      pdb::Factory::ResetFlagDefaults();

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppRestraintSaxs

  const ExampleClass::EnumType ExampleAppRestraintSaxs::s_Instance
  (
    GetExamples().AddEnum( ExampleAppRestraintSaxs())
  );

} // namespace bcl
