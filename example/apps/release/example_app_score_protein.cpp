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
//#include "release/bcl_app_score_protein.cpp"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "command/bcl_command_command.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_score_protein.cpp
  //!
  //! @author woetzen
  //! @date November 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppScoreProtein :
    public ExampleInterface
  {
  public:

    ExampleAppScoreProtein *Clone() const
    {
      return new ExampleAppScoreProtein( *this);
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

      const app::ApplicationType app_enum_score_protein( "ScoreProtein");
      BCL_ExampleAssert( app_enum_score_protein.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // ScoreProtein
      // -pdb 1ubi.pdb
      // -native 1ubi.pdb
      // -quality RMSD GDT_TS
      // -sspred JUFO PSIPRED
      // -sequence_data [example_path] 1ubi
      // -score_density_agreement 1ubi_res_6.6voxelsize_2.200Gaussian.mrc 6.6 GaussianSphere CCCScaled
      // -score_table_write scores.table
      {
        ApplicationExampleHelper score_protein_helper( app_enum_score_protein);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string ss_pred_path( AddExampleInputPathToFilename( e_Biology, ""));
        const std::string ss_prefix( "1ubi");
        const std::string mrc_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi_res_6.6voxelsize_2.200Gaussian.mrc"));
        const std::string output_table
        (
          AddExampleOutputPathToFilename( score::GetNamespaceIdentifier(), "scores.table")
        );

        // flags
        score_protein_helper.SetFlag( "pdb", pdb_file_name);
        score_protein_helper.SetFlag( "native", pdb_file_name);
        score_protein_helper.SetFlag( "quality", storage::Vector< std::string>::Create( "RMSD", "GDT_TS"));
        score_protein_helper.SetFlag( "sspred", storage::Vector< std::string>::Create( "JUFO", "PSIPRED"));
        score_protein_helper.SetFlag( "sequence_data", storage::Vector< std::string>::Create( ss_pred_path, ss_prefix));
        score_protein_helper.SetFlag( "score_density_agreement", storage::Vector< std::string>::Create( mrc_file_name, "6.6", "GaussianSphere", "CCCScaled"));
        score_protein_helper.SetFlag( "score_table_write", output_table);

        // check the command line
        BCL_ExampleAssert( score_protein_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( score_protein_helper.RunCommand(), 0);

        // output
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, output_table);
        storage::Table< double> result;
        result.ReadFormatted( read);
        io::File::CloseClearFStream( read);

        // some simple checks
        BCL_ExampleCheck( result.GetSize(), 1);
      }

      // ScoreProtein
      // -pdblist cluster_20pdbs.ls
      // -native 1ubi.pdb
      // -quality RMSD GDT_TS
      // -sspred JUFO PSIPRED
      // -sequence_data [example_path] 1ubi
      // -score_density_agreement 1ubi_res_6.6voxelsize_2.200Gaussian.mrc 6.6 GaussianSphere CCCScaled
      // -score_table_write scores20.table
      {
        ApplicationExampleHelper score_protein_helper( app_enum_score_protein);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string ss_pred_path( AddExampleInputPathToFilename( e_Biology, ""));
        const std::string ss_prefix( "1ubi");
        const std::string pdb_list( AddExampleInputPathToFilename( e_Cluster, "cluster_20pdbs.ls"));
        const std::string mrc_file_name( AddExampleInputPathToFilename( e_Biology, "1ubi_res_6.6voxelsize_2.200Gaussian.mrc"));
        const std::string output_table
        (
          AddExampleOutputPathToFilename( score::GetNamespaceIdentifier(), "scores20.table")
        );

        // flags
        score_protein_helper.SetFlag( "pdblist", pdb_list);
        score_protein_helper.SetFlag( "native", pdb_file_name);
        score_protein_helper.SetFlag( "quality", storage::Vector< std::string>::Create( "RMSD", "GDT_TS"));
        score_protein_helper.SetFlag( "sspred", storage::Vector< std::string>::Create( "JUFO", "PSIPRED"));
        score_protein_helper.SetFlag( "sequence_data", storage::Vector< std::string>::Create( ss_pred_path, ss_prefix));
        score_protein_helper.SetFlag( "score_density_agreement", storage::Vector< std::string>::Create( mrc_file_name, "6.6", "GaussianSphere", "CCCScaled"));
        score_protein_helper.SetFlag( "score_table_write", output_table);

        // check the command line
        BCL_ExampleAssert( score_protein_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( score_protein_helper.RunCommand(), 0);

        // output
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, output_table);
        storage::Table< double> result;
        result.ReadFormatted( read);
        io::File::CloseClearFStream( read);

        // some simple checks
        BCL_ExampleCheck( result.GetSize(), 20);
      }

      // reset all pdb factory flags, since this application changes them
      pdb::Factory::ResetFlagDefaults();

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppScoreProtein

  const ExampleClass::EnumType ExampleAppScoreProtein::s_Instance
  (
    GetExamples().AddEnum( ExampleAppScoreProtein())
  );

} // namespace bcl
