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
//#include "release/bcl_app_fit_in_density.cpp"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_fit_in_density.cpp
  //!
  //! @author woetzen
  //! @date November 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppFitInDensity :
    public ExampleInterface
  {
  public:

    ExampleAppFitInDensity *Clone() const
    {
      return new ExampleAppFitInDensity( *this);
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

      const app::ApplicationType app_enum_fit_in_density( "FitInDensity");
      BCL_ExampleAssert( app_enum_fit_in_density.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // FitInDensity example/example_files/input/biology/1ubi.pdb example/example_files/input/biology/1ubi_res_6.6voxelsize_2.200Gaussian.mrc
      // -mrc_resolution 6.6
      // -hash_storage HashMap
      // -prefix example/example_files/output/biol/fit_in_density/
      // -protein_storage example/example_files/output/biol/fit_in_density/ Overwrite
      // -coordinatesystem Spherical
      // -atoms CA N C O
      {
        ApplicationExampleHelper fit_in_density_helper( app_enum_fit_in_density);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string mrc_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi_res_6.6voxelsize_2.200Gaussian.mrc"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "fit_in_density/")
        );

        // flags
        fit_in_density_helper.AddParameter( pdb_file_name);
        fit_in_density_helper.AddParameter( mrc_file_name);
        fit_in_density_helper.SetFlag( "mrc_resolution", "6.6");
        fit_in_density_helper.SetFlag( "hash_storage", "HashMap");
        fit_in_density_helper.SetFlag( "prefix", output_prefix);
        fit_in_density_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( output_prefix, "Overwrite"));
        fit_in_density_helper.SetFlag( "coordinatesystem", "Spherical");
        fit_in_density_helper.SetFlag( "atoms", storage::Vector< std::string>::Create( "CA", "N", "C", "O"));

        // check the command line
        BCL_ExampleAssert( fit_in_density_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( fit_in_density_helper.RunCommand(), 0);

        // output
        const std::string result_table_name( output_prefix + "result.table");
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, result_table_name);
        storage::Table< double> result;
        result.ReadFormatted( read);
        io::File::CloseClearFStream( read);

        // stored results
        const util::ShPtr< assemble::ProteinStorageFile> sp_storage
        (
          new assemble::ProteinStorageFile( output_prefix)
        );

        // some simple checks
        BCL_ExampleCheck( sp_storage->GetSize( "transformed"), 10);
        BCL_ExampleCheck( sp_storage->GetSize( "transformed_min"), 10);
        BCL_ExampleCheck( result.GetSize(), 10);
        BCL_ExampleCheckWithinTolerance( 1.0, result.Begin()->Second()[ "min_corr"], 0.002);
        BCL_ExampleCheck( result.Begin()->Second()[ "min_RMSD"] < 0.25, true);

        // cleanup
        io::Directory output_dir( output_prefix);
        BCL_ExampleCheck( output_dir.Remove( true), true);
      }

      // reset all pdb factory flags, since this application changes them
      pdb::Factory::ResetFlagDefaults();

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppFitInDensity

  const ExampleClass::EnumType ExampleAppFitInDensity::s_Instance
  (
    GetExamples().AddEnum( ExampleAppFitInDensity())
  );

} // namespace bcl
