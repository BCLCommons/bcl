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
//#include "release/bcl_app_fold.cpp"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "fold/bcl_fold.h"
#include "io/bcl_io_directory_entry.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_fold.cpp
  //!
  //! @author woetzen, weinerbe
  //! @date June 5, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppFold :
    public ExampleInterface
  {
  public:

    ExampleAppFold *Clone() const
    {
      return new ExampleAppFold( *this);
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

      const app::ApplicationType app_enum_fold( "Fold");
      BCL_ExampleAssert( app_enum_fold.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // example with native pool
      // Fold
      // -native example/example_files/input/biology/1ubi.pdb
      // -use_native_pool
      // -diable_sse_resize
      // -mc_number_iterations 5000 400
      // -sspred JUFO PSIPRED
      // -sequence_data example/example_files/input/biology/ 1ubi
      // -prefix 1ubi_test_
      // -protein_storage example/example_files/output/fold/ Overwrite
      {
        ApplicationExampleHelper fold_helper( app_enum_fold);

        // get paths and filenames
        const std::string native_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string storage_path( AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), ""));
        const std::string prefix( "1ubi_test_");
        const std::string output_pdb
        (
          AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), prefix + "final_0000.pdb")
        );
        const std::string ss_path( AddExampleInputPathToFilename( e_Biology, ""));

        // test using native pool, default scores, and single stage
        fold_helper.SetFlag( "native", native_pdb);
        fold_helper.SetFlag( "use_native_pool");
        fold_helper.SetFlag( "mc_number_iterations", storage::Vector< std::string>::Create( "5000", "400"));
        fold_helper.SetFlag( "sspred", storage::Vector< std::string>::Create( "JUFO", "PSIPRED"));
        fold_helper.SetFlag( "sequence_data", storage::Vector< std::string>::Create( ss_path, "1ubi"));
        fold_helper.SetFlag( "prefix", prefix);
        fold_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( storage_path, "Overwrite"));

        // check the command line and run
        BCL_ExampleAssert( fold_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fold_helper.RunCommand(), 0);

        // check the output pdb
        const size_t expected_helices( 1);
        const size_t expected_strands( 4);
        assemble::ProteinModel protein_model( Proteins::GetModel( output_pdb));
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().HELIX).GetSize(), expected_helices);
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().STRAND).GetSize(), expected_strands);

        // clean up
        io::DirectoryEntry output_file( output_pdb);
        output_file.Remove();
      }

      // example with predicted pool and multiple stages
      // Fold
      // -native example/example_files/input/biology/1ubi.pdb
      // -quality RMSD GDT_TS
      // -pool example/example_files/input/biology/1ubiA.SSPredMC_PSIPRED_JUFO.pool
      // -pool_separate
      // -pool_min_sse_lengths 5 3
      // -function_cache
      // -sspred JUFO PSIPRED
      // -sequence_data example/example_files/input/biology/ 1ubi
      // -stages_read example/example_files/input/fold/fold_stages.txt
      // -prefix 1ubi_test_
      // -protein_storage example/example_files/output/fold/ Overwrite
      {
        ApplicationExampleHelper fold_helper( app_enum_fold);

        // get paths and filenames
        const std::string native_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string storage_path( AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), ""));
        const std::string prefix( "1ubi_test_");
        const std::string output_pdb
        (
          AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), prefix + "final_0000.pdb")
        );
        const std::string ss_path( AddExampleInputPathToFilename( e_Biology, ""));
        const std::string pool( AddExampleInputPathToFilename( e_Biology, "1ubiA.SSPredMC_PSIPRED_JUFO.pool"));
        const std::string stages( AddExampleInputPathToFilename( e_Fold, "fold_stages.txt"));

        // test using predicted pool and multiple stages
        fold_helper.SetFlag( "native", native_pdb);
        fold_helper.SetFlag( "quality", storage::Vector< std::string>::Create( "RMSD", "GDT_TS"));
        fold_helper.SetFlag( "pool", pool);
        fold_helper.SetFlag( "pool_separate");
        fold_helper.SetFlag( "pool_min_sse_lengths", storage::Vector< std::string>::Create( "5", "3"));
        fold_helper.SetFlag( "sspred", storage::Vector< std::string>::Create( "JUFO", "PSIPRED"));
        fold_helper.SetFlag( "sequence_data", storage::Vector< std::string>::Create( ss_path, "1ubi"));
        fold_helper.SetFlag( "stages_read", stages);
        fold_helper.SetFlag( "prefix", prefix);
        fold_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( storage_path, "Overwrite"));

        // check the command line and run
        BCL_ExampleAssert( fold_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fold_helper.RunCommand(), 0);

        // check the output pdb
        const size_t expected_helices( 1);
        const size_t expected_strands( 4);
        assemble::ProteinModel protein_model( Proteins::GetModel( output_pdb));
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().HELIX).GetSize() >= expected_helices, true);
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().STRAND).GetSize() >= expected_strands, true);

        // clean up
        io::DirectoryEntry output_file( output_pdb);
        output_file.Remove();
      }

      // EMFold example
      // Fold
      // -protocols EM
      // -nmodels 1
      // -fasta example/example_files/input/biology/1IE9A.fasta
      // -pool example/example_files/output/assemble/1IE9.pool
      // -mc_number_iterations 2000 400
      // -mc_temperature_fraction 0.25 0.05
      // -protein_storage example/example_files/output/fold/emfold/ Overwrite
      // -body_restraint example/example_files/input/biology/1IE9.pdb
      // -print_body_assignment
      // -score_weightset_read example/example_files/input/fold/assembly_weightset.table
      // -score_density_connectivity example/example_files/output/fold/1IE9.mrc
      // -write_minimization improved
      // -sspred JUFO PSIPRED
      // -sequence_data example/example_files/input/biology/ 1IE9
      {
        ApplicationExampleHelper fold_helper( app_enum_fold);
        // get input and output file names
        const std::string pdb_file_name( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
        const std::string mrc_file_name( AddExampleInputPathToFilename( e_Fold, "1IE9.mrc"));
        const std::string fasta_file_name( AddExampleInputPathToFilename( e_Biology, "1IE9A.fasta"));
        const std::string pool_file_name( AddExampleInputPathToFilename( e_Fold, "1IE9.pool"));
        const std::string weightset_file_name( AddExampleInputPathToFilename( e_Fold, "emfold_assembly.weightset"));
        const std::string output_prefix
        (
          AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), "emfold/")
        );

        // flags
        fold_helper.SetFlag( "protocols", "EM");
        fold_helper.SetFlag( "nmodels", "1");
        fold_helper.SetFlag( "fasta", fasta_file_name);
        fold_helper.SetFlag( "pool", pool_file_name);
        fold_helper.SetFlag( "mc_number_iterations", storage::Vector< std::string>::Create( "1000", "200"));
        fold_helper.SetFlag( "mc_temperature_fraction", storage::Vector< std::string>::Create( "0.25", "0.05"));
        fold_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( output_prefix, "Overwrite"));
        fold_helper.SetFlag( "body_restraint", pdb_file_name);
        fold_helper.SetFlag( "print_body_assignment");
        fold_helper.SetFlag( "score_weightset_read", weightset_file_name);
        fold_helper.SetFlag( "score_density_connectivity", mrc_file_name);
        fold_helper.SetFlag( "write_minimization", "improved");
        fold_helper.SetFlag( "sspred", storage::Vector< std::string>::Create( "JUFO", "PSIPRED"));
        fold_helper.SetFlag
        (
          "sequence_data",
          storage::Vector< std::string>::Create( AddExampleInputPathToFilename( e_Biology, ""), "1IE9")
        );

        // check the command line
        BCL_ExampleAssert( fold_helper.CheckCommandString( true), true);

        // run
        BCL_ExampleAssert( fold_helper.RunCommand(), 0);

        // output
        const std::string output_file( output_prefix + "final_0000.pdb");
        const assemble::ProteinModel correct_model( Proteins::GetModel( pdb_file_name));
        const assemble::ProteinModel final_model( Proteins::GetModel( output_file));

        const double gdt4A( assemble::Quality::Calculate( quality::GetMeasures().e_GDT_4A, final_model, correct_model, storage::Set< biol::AtomType>( biol::GetAtomTypes().CA)));
        BCL_MessageStd( "gdt4A of final model to correct structure: " + util::Format()( gdt4A));

        // check # of sses
        if( BCL_ExampleCheck( final_model.GetNumberSSEs() >= 15, true))
        {
          // cleanup
          io::Directory output_dir( output_prefix);
          BCL_ExampleCheck( output_dir.Remove( true), true);
        }
      }

      // MP-Fold example
      // -native example/example_files/input/biology/4A2N.pdb
      // -use_native_pool
      // -protocols Membrane
      // -membrane
      // -tm_helices example/example_files/input/biology/4A2N.SSPredHighest_OCTOPUS.pool
      // -diable_sse_resize
      // -mc_number_iterations 5000 400
      // -prefix 4A2N_test_
      // -protein_storage example/example_files/output/fold/ Overwrite
      {
        ApplicationExampleHelper fold_helper( app_enum_fold);

        // get paths and filenames
        const std::string native_pdb( AddExampleInputPathToFilename( e_Biology, "4A2N.pdb"));
        const std::string tm_pool( AddExampleInputPathToFilename( e_Biology, "4A2N.SSPredHighest_OCTOPUS.pool"));
        const std::string storage_path( AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), ""));
        const std::string prefix( "4A2N_test_");
        const std::string output_pdb
        (
          AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), prefix + "final_0000.pdb")
        );
        const std::string ss_path( AddExampleInputPathToFilename( e_Biology, ""));

        // test using native pool, default scores, and single stage
        fold_helper.SetFlag( "native", native_pdb);
        fold_helper.SetFlag( "use_native_pool");
        fold_helper.SetFlag( "mc_number_iterations", storage::Vector< std::string>::Create( "10000", "800"));
        fold_helper.SetFlag( "membrane");
        fold_helper.SetFlag( "tm_helices", tm_pool);
        fold_helper.SetFlag( "protocols", "Membrane");
        fold_helper.SetFlag( "prefix", prefix);
        fold_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( storage_path, "Overwrite"));

        // check the command line and run
        BCL_ExampleAssert( fold_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fold_helper.RunCommand(), 0);

        // check the output pdb
        const size_t expected_helices( 5);
        assemble::ProteinModel protein_model( Proteins::GetModel( output_pdb));
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().HELIX).GetSize(), expected_helices);

        // clean up
        io::DirectoryEntry output_file( output_pdb);
        output_file.Remove();
      }

      // MP-Fold multimer example
      // -native example/example_files/input/biology/1J4N.pdb
      // -use_native_pool
      // -protocols Membrane Multimer
      // -membrane
      // -tm_helices example/example_files/input/biology/1J4N.SSPredMC_OCTOPUS.pool
      // -symmetry C4
      // -diable_sse_resize
      // -mc_number_iterations 3000 400
      // -prefix 1J4N_test_
      // -protein_storage example/example_files/output/fold/ Overwrite
      {
        ApplicationExampleHelper fold_helper( app_enum_fold);

        // get paths and filenames
        const std::string native_pdb( AddExampleInputPathToFilename( e_Biology, "1J4N.pdb"));
        const std::string tm_pool( AddExampleInputPathToFilename( e_Biology, "1J4N.SSPredMC_OCTOPUS.pool"));
        const std::string storage_path( AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), ""));
        const std::string prefix( "1J4N_test_");
        const std::string output_pdb
        (
          AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), prefix + "final_0000.pdb")
        );
        const std::string ss_path( AddExampleInputPathToFilename( e_Biology, ""));

        // test using native pool, default scores, and single stage
        fold_helper.SetFlag( "native", native_pdb);
        fold_helper.SetFlag( "use_native_pool");
        fold_helper.SetFlag( "mc_number_iterations", storage::Vector< std::string>::Create( "3000", "400"));
        fold_helper.SetFlag( "membrane");
        fold_helper.SetFlag( "tm_helices", tm_pool);
        fold_helper.SetFlag( "symmetry", "C4");
        fold_helper.SetFlag( "protocols", storage::Vector< std::string>::Create( "Membrane", "Multimer"));
        fold_helper.SetFlag( "prefix", prefix);
        fold_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( storage_path, "Overwrite"));

        // check the command line and run
        BCL_ExampleAssert( fold_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fold_helper.RunCommand(), 0);

        // check the output pdb
        const size_t expected_helices( 3);
        assemble::ProteinModel protein_model( Proteins::GetModel( output_pdb));
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().HELIX).GetSize() >= expected_helices, true);

        // clean up
        io::DirectoryEntry output_file( output_pdb);
        output_file.Remove();
      }

      // NMR restraint example
      // -native example/example_files/input/biology/1ubi.pdb
      // -use_native_pool
      // -sspred JUFO PSIPRED
      // -sequence_data example/example_files/input/biology/ 1ubi
      // -protocols Restraint
      // -restraint_types NOE
      // -restraint_prefix example/example_files/input/biology/1ubi
      // -diable_sse_resize
      // -mc_number_iterations 5000 400
      // -prefix 1ubi_restraint_test_
      // -protein_storage example/example_files/output/fold/ Overwrite
      {
        ApplicationExampleHelper fold_helper( app_enum_fold);

        // get paths and filenames
        const std::string native_pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
        const std::string restraint_prefix( AddExampleInputPathToFilename( e_Biology, "1ubi"));
        const std::string storage_path( AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), ""));
        const std::string prefix( "1ubi_restraint_test_");
        const std::string output_pdb
        (
          AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), prefix + "final_0000.pdb")
        );
        const std::string ss_path( AddExampleInputPathToFilename( e_Biology, ""));

        // test using native pool, default scores, and single stage
        fold_helper.SetFlag( "native", native_pdb);
        fold_helper.SetFlag( "use_native_pool");
        fold_helper.SetFlag( "sspred", storage::Vector< std::string>::Create( "JUFO", "PSIPRED"));
        fold_helper.SetFlag( "sequence_data", storage::Vector< std::string>::Create( ss_path, "1ubi"));
        fold_helper.SetFlag( "protocols", "Restraint");
        fold_helper.SetFlag( "restraint_types", "NOE(Star)");
        fold_helper.SetFlag( "restraint_prefix", restraint_prefix);
        fold_helper.SetFlag( "mc_number_iterations", storage::Vector< std::string>::Create( "5000", "400"));
        fold_helper.SetFlag( "prefix", prefix);
        fold_helper.SetFlag( "protein_storage", storage::Vector< std::string>::Create( storage_path, "Overwrite"));

        // check the command line and run
        BCL_ExampleAssert( fold_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fold_helper.RunCommand(), 0);

        // check the output pdb
        const size_t expected_helices( 1);
        const size_t expected_strands( 4);
        assemble::ProteinModel protein_model( Proteins::GetModel( output_pdb));
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().HELIX).GetSize(), expected_helices);
        BCL_ExampleCheck( protein_model.GetSSEs( biol::GetSSTypes().STRAND).GetSize(), expected_strands);

        // clean up
        io::DirectoryEntry output_file( output_pdb);
        output_file.Remove();
      }

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppFold

  const ExampleClass::EnumType ExampleAppFold::s_Instance
  (
    GetExamples().AddEnum( ExampleAppFold())
  );

} // namespace bcl
