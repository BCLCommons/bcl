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
//#include "release/bcl_app_fit_in_density_minimize.cpp"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "density/bcl_density_fit_protein_minimizers.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_protein_agreement_ccc.h"
#include "density/bcl_density_simulators.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_fit_in_density_minimize.cpp
  //!
  //! @author woetzen
  //! @date January 27, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppFitInDensityMinimize :
    public ExampleInterface
  {
  public:

    ExampleAppFitInDensityMinimize *Clone() const
    {
      return new ExampleAppFitInDensityMinimize( *this);
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

      const app::ApplicationType app_enum_fit_in_density( "FitInDensityMinimize");
      BCL_ExampleAssert( app_enum_fit_in_density.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // iterate over minimizers
      for( density::FitProteinMinimizers::const_iterator itr( density::GetFitProteinMinimizers().Begin()), itr_end( density::GetFitProteinMinimizers().End()); itr != itr_end; ++itr)
      {
        if( !( *itr)->IsDefined())
        {
          continue;
        }
        ApplicationExampleHelper fit_in_density_helper( app_enum_fit_in_density);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));
        const std::string mrc_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi_res_6.6voxelsize_2.200Gaussian.mrc"));
        util::ShPtr< io::Directory> output_directory
        (
          new io::Directory( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "fit_in_density_minimize"))
        );

        BCL_ExampleIndirectAssert( output_directory->Make(), true, "create directory for output: " + output_directory->GetPath());

        // flags
        fit_in_density_helper.AddParameter( pdb_file_name);
        fit_in_density_helper.AddParameter( mrc_file_name);
        fit_in_density_helper.SetFlag( "mrc_resolution", "6.6");
        fit_in_density_helper.SetFlag( "prefix", output_directory->GetPath() + '/');
        fit_in_density_helper.SetFlag( "approximator", itr->GetName());

        // check the command line
        BCL_ExampleAssert( fit_in_density_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( fit_in_density_helper.RunCommand(), 0);

        // output
        const io::DirectoryEntry result_file_name( output_directory, "transformed_min.pdb");
        BCL_ExampleIndirectAssert( result_file_name.DoesExist(), true, "check if " + result_file_name.GetFullName() + " was written");
        const assemble::ProteinModel fitted_structure( Proteins::GetModel( result_file_name.GetFullName(), biol::GetAAClasses().e_AAComplete));

        // cleanup
        BCL_ExampleCheck( output_directory->Remove( true), true);

        // checks for structure similarity
        const assemble::ProteinModel starting_structure( Proteins::GetModel( pdb_file_name, biol::GetAAClasses().e_AAComplete));

        // rmsd
        const double rmsd( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD_NoSuperimposition, fitted_structure, starting_structure, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)));

        // ccc
        const double resolution( 6.6), voxelsize( 2.2);
        const util::ShPtr< density::SimulateInterface> simulator
        (
          density::GetSimulators().CreateSimulator
          (
            density::GetSimulators().e_Gaussian, linal::Vector3D( voxelsize), resolution
          )
        );
        // density map of original model
        util::ShPtr< density::Map> sp_density_map( new density::Map());
        io::IFStream read;

        BCL_ExampleMustOpenBinaryInputFile( read, mrc_file_name);
        sp_density_map->ReadMRC( read);
        io::File::CloseClearFStream( read);
        density::ProteinAgreementCCC agreement_ccc( false, false);
        agreement_ccc.SetDensityMap( sp_density_map);
        agreement_ccc.SetSimulator( simulator);

        const double ccc( agreement_ccc( fitted_structure));
        BCL_MessageStd( "ccc: " + util::Format()( ccc) + " rmsd: " + util::Format()( rmsd));
        BCL_ExampleCheck( rmsd < 0.5, true);
        BCL_ExampleCheck( ccc < -0.7, true);
      }

      // reset all pdb factory flags, since this application changes them
      pdb::Factory::ResetFlagDefaults();

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppFitInDensityMinimize

  const ExampleClass::EnumType ExampleAppFitInDensityMinimize::s_Instance
  (
    GetExamples().AddEnum( ExampleAppFitInDensityMinimize())
  );

} // namespace bcl
