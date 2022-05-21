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
//#include "internal/biol/bcl_app_fusion_protein.cpp"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "example_interface.h"
#include "biol/bcl_biol.h"
#include "io/bcl_io_directory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_fusion_protein.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Sep 25, 2012
  //! @remarks status incomplete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppFusionProtein :
    public ExampleInterface
  {
  public:

    ExampleAppFusionProtein *Clone() const
    {
      return new ExampleAppFusionProtein( *this);
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

      const app::ApplicationType app_enum_fold( "protein:FusionProtein");
      BCL_ExampleAssert( app_enum_fold.IsDefined(), true);

      // input folder
      const std::string input_path( AddExampleInputPathToFilename( e_Biology, "fusion_protein/"));

    ////////////////
    // operations //
    ////////////////

      // replace term of 1gl4 with parts of 1b9c
      // FusionProtein -output_prefix 1gl4_1b9c
      // -scaffold 1gl4A.pdb -donor data/1gl4A.pdb
      // -replace_term 1b9cA.pdb locator_sse_1gl4_strand1.bcl 1b9cA_strand1.bcl 1gl4A.pdb locator_sse_1gl4_strand_last.bcl locator_sse_1gl4_strand_last.bcl
      // -quality RMSD 0.5 -assignment_score identity -min_number_aligned_residues 3
      // -cutpoint_optimal_peptide
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1gl4A.pdb");
        const std::string donor_pdb( input_path + "1b9cA.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1gl4_1b9c/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "donor", scaffold_pdb);
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "replace_term", storage::Vector< std::string>::Create( donor_pdb, input_path + "locator_sse_1gl4_strand1.bcl", input_path + "1b9cA_strand1.bcl", scaffold_pdb, input_path + "locator_sse_1gl4_strand_last.bcl", input_path + "locator_sse_1gl4_strand_last.bcl"));
        fusion_helper.SetFlag( "quality", storage::Vector< std::string>::Create( "RMSD", "0.5"));
        fusion_helper.SetFlag( "assignment_score", "identity");
        fusion_helper.SetFlag( "min_number_aligned_residues", "3");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop5 loop6  from 1thf onto 1qo2
      // FusionProtein -output_prefix 1thf_1qo2/
      // -scaffold 1qo2.pdb -scaffold_fragment 1qo2_loop1.locators 1qo2_loop5.locators 1qo2_loop6.locators
      // -donor 1thf.pdb 1thf.pdb 1thf.pdb -donor_sites 1thf_loop1.site 1thf_loop5.site 1thf_loop6.site
      // -cutpoint_optimal_peptide
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1qo2.pdb");
        const std::string donor_pdb( input_path + "1thf.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1thf_1qo2/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1qo2_loop1.locators", input_path + "1qo2_loop5.locators", input_path + "1qo2_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 3, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1thf_loop1.site", input_path + "1thf_loop5.site", input_path + "1thf_loop6.site"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "write_rosetta_res_file");
        fusion_helper.SetFlag( "cross_over_heatmap");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop2 loop3 loop5 loop6  from 1qo2 onto 1thf
      // FusionProtein -output_prefix 1qo2_1thf/
      // -scaffold 1thf.pdb -scaffold_fragment 1thf_loop1.locators 1thf_loop2.locators 1thf_loop3.locators 1thf_loop5.locators 1thf_loop6.locators
      // -donor 1qo2.pdb 1qo2.pdb 1qo2.pdb 1qo2.pdb 1qo2.pdb
      // -donor_sites 1qo2_loop1.site 1qo2_loop2.site 1qo2_loop3.site 1qo2_loop5.site 1qo2_loop6.site
      // -cutpoint_optimal_peptide
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1thf.pdb");
        const std::string donor_pdb( input_path + "1qo2.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1qo2_1thf/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1thf_loop1.locators", input_path + "1thf_loop2.locators", input_path + "1thf_loop3.locators", input_path + "1thf_loop5.locators", input_path + "1thf_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 5, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1qo2_loop1.site", input_path + "1qo2_loop2.site", input_path + "1qo2_loop3.site", input_path + "1qo2_loop5.site", input_path + "1qo2_loop6.site"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // create heatmaps
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1thf.pdb");
        const std::string donor_pdb( input_path + "1qo2.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1qo2_1thf/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1thf_loop1.locators", input_path + "1thf_loop2.locators", input_path + "1thf_loop3.locators", input_path + "1thf_loop5.locators", input_path + "1thf_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 5, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1qo2_loop1.site", input_path + "1qo2_loop2.site", input_path + "1qo2_loop3.site", input_path + "1qo2_loop5.site", input_path + "1qo2_loop6.site"));
        fusion_helper.SetFlag( "cross_over_heatmap");
        fusion_helper.SetFlag( "write_selected_models");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop4 loop6  from 1nsj onto 1thf
      // FusionProtein -output_prefix 1nsj_1thf/
      // -scaffold 1thf.pdb -scaffold_fragment 1thf_loop1.locators 1thf_loop4.locators 1thf_loop6.locators
      // -donor 1nsj.pdb 1nsj.pdb 1nsj.pdb -donor_sites 1nsj_loop1.site 1nsj_loop4.site 1nsj_loop6.site
      // -cutpoint_optimal_peptide
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1thf.pdb");
        const std::string donor_pdb( input_path + "1nsj.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1nsj_1thf/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1thf_loop1.locators", input_path + "1thf_loop4.locators", input_path + "1thf_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 3, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1nsj_loop1.site", input_path + "1nsj_loop4.site", input_path + "1nsj_loop6.site"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop4 loop6  from 1nsj onto 1qo2
      // FusionProtein -output_prefix 1nsj_1qo2/
      // -scaffold 1qo2.pdb -scaffold_fragment 1qo2_loop1.locators 1qo2_loop4.locators 1qo2_loop6.locators
      // -donor 1nsj.pdb 1nsj.pdb 1nsj.pdb -donor_sites 1nsj_loop1.site 1nsj_loop4.site 1nsj_loop6.site
      // -cutpoint_optimal_peptide
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1qo2.pdb");
        const std::string donor_pdb( input_path + "1nsj.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1nsj_1qo2/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1qo2_loop1.locators", input_path + "1qo2_loop4.locators", input_path + "1qo2_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 3, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1nsj_loop1.site", input_path + "1nsj_loop4.site", input_path + "1nsj_loop6.site"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop5 loop6  from 1a53 onto 1thf
      // FusionProtein -output_prefix 1a53_1thf/
      // -scaffold 1thf.pdb -scaffold_fragment 1thf_loop1.locators 1thf_loop5.locators 1thf_loop6.locators
      // -donor 1a53.pdb 1a53.pdb 1a53.pdb -donor_sites 1a53_loop1.site 1a53_loop5.site 1a53_loop6.site
      // -cutpoint_optimal_peptide
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1thf.pdb");
        const std::string donor_pdb( input_path + "1a53.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1a53_1thf/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1thf_loop1.locators", input_path + "1thf_loop5.locators", input_path + "1thf_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 3, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1a53_loop1.site", input_path + "1a53_loop5.site", input_path + "1a53_loop6.site"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop4 loop6  from 1nsj onto 1a53
      // FusionProtein -output_prefix 1nsj_1a53/
      // -scaffold 1a53.pdb -scaffold_fragment 1a53_loop1.locators 1a53_loop4.locators 1a53_loop6.locators
      // -donor 1nsj.pdb 1nsj.pdb 1nsj.pdb -donor_sites 1nsj_loop1.site 1nsj_loop4.site 1nsj_loop6.site @general.flags
      // -cutpoint_optimal_peptide
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1a53.pdb");
        const std::string donor_pdb( input_path + "1nsj.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1nsj_1a53/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1a53_loop1.locators", input_path + "1a53_loop4.locators", input_path + "1a53_loop6.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 3, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1nsj_loop1.site", input_path + "1nsj_loop4.site", input_path + "1nsj_loop6.site"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "cutpoint_optimal_peptide");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

      // transferring: loop1 loop4 loop6  from 1nsj onto 1thf with every possible cutpoint, only writing model 0 and 4
      // FusionProtein -output_prefix 1nsj_1thf/
      // -scaffold 1thf.pdb -scaffold_fragment 1thf_loop1.locators 1thf_loop4.locators 1thf_loop6.locators
      // -donor 1nsj.pdb 1nsj.pdb 1nsj.pdb -donor_sites 1nsj_loop1.site 1nsj_loop4.site 1nsj_loop6.site
      // -write_selected_models 0 4
      // @general.flags
      {
        ApplicationExampleHelper fusion_helper( app_enum_fold);

        // get paths and filenames
        const std::string scaffold_pdb( input_path + "1thf.pdb");
        const std::string donor_pdb( input_path + "1nsj.pdb");
        const std::string prefix( AddExampleOutputPathToFilename( biol::GetNamespaceIdentifier(), "1nsj_1thf_combi/"));
        io::Directory output_prefix( prefix);
        if( !output_prefix.DoesExist())
        {
          BCL_ExampleAssert( output_prefix.Make(), true);
        }

        // test
        fusion_helper.SetFlag( "scaffold", scaffold_pdb);
        fusion_helper.SetFlag( "scaffold_fragment", storage::Vector< std::string>::Create( input_path + "1thf_loop1.locators", input_path + "1thf_loop4.locators"));
        fusion_helper.SetFlag( "donor", storage::Vector< std::string>( 2, donor_pdb));
        fusion_helper.SetFlag( "donor_sites", storage::Vector< std::string>::Create( input_path + "1nsj_loop1.site", input_path + "1nsj_loop4.site"));
        fusion_helper.SetFlag( "write_selected_models", storage::Vector< std::string>::Create( "0", "4"));
        fusion_helper.SetFlag( "replace_combinatorial");
        fusion_helper.SetFlag( "output_prefix", prefix);
        fusion_helper.AddParameter( "@" + input_path + "general.flags");

        // check the command line and run
        BCL_ExampleAssert( fusion_helper.CheckCommandString( true), true);
        BCL_ExampleAssert( fusion_helper.RunCommand(), 0);
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppFusionProtein

  const ExampleClass::EnumType ExampleAppFusionProtein::s_Instance
  (
    GetExamples().AddEnum( ExampleAppFusionProtein())
  );

} // namespace bcl
