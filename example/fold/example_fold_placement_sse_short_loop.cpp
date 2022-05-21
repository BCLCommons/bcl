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
#include "fold/bcl_fold_placement_sse_short_loop.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_placement_sse_short_loop.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPlacementSSEShortLoop :
    public ExampleInterface
  {
  public:

    ExampleFoldPlacementSSEShortLoop *Clone() const
    {
      return new ExampleFoldPlacementSSEShortLoop( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));
      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // create an empty protein model
      BCL_MessageStd( "Initializing chain information for empty model");
      assemble::ProteinModel initial_model( model.GetEmptyChains());

//      // Please do not remove the commented sections below, this is required for detailed testing
//      // Mert
//      {
//        // read in the protein model 1J27
//        const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1J27A.pdb"));
//        storage::Map< biol::SSType, size_t> ssetype_min_size;
//        ssetype_min_size[ biol::GetSSTypes().HELIX] = 7;
//        ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;
//        assemble::ProteinModel
//          protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));
//
//        // create the sse pool
//        BCL_MessageStd( "Creating the SSEPool");
//        util::ShPtr< assemble::SSEPool> sse_pool( new assemble::SSEPool( protein_model.GetSSEs()));
//
//        // locate SSEs
//        util::SiPtr< const assemble::SSE> si_helix_18_37( assemble::LocatorSSE( 'A', 18, 37).Locate( protein_model));
//        util::SiPtr< const assemble::SSE> si_helix_63_81( assemble::LocatorSSE( 'A', 63, 81).Locate( protein_model));
//        util::SiPtr< const assemble::SSE> si_strand_2_13( assemble::LocatorSSE( 'A',  2, 13).Locate( protein_model));
//        util::SiPtr< const assemble::SSE> si_strand_40_45( assemble::LocatorSSE( 'A', 40, 45).Locate( protein_model));
//        util::SiPtr< const assemble::SSE> si_strand_52_61( assemble::LocatorSSE( 'A', 52, 61).Locate( protein_model));
//        util::SiPtr< const assemble::SSE> si_strand_85_97( assemble::LocatorSSE( 'A', 85, 97).Locate( protein_model));
//
//      //////////////////////////////////
//      // construction and destruction //
//      //////////////////////////////////
//
//        // construct placement using default constructor
//        fold::PlacementSSEShortLoop place_default;
//
//        // construct from a locator
//        fold::PlacementSSEShortLoop place_side( 7, 0.0);
//        fold::PlacementSSEShortLoop place_top( 7, 1.0);
//
//        // create an empty model
//        assemble::ProteinModel model;
//        model.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( protein_model.GetChain( 'A')->GetSequence())));
//
//        // make copies of SSEs
//        util::ShPtr< assemble::SSE> sp_helix_18_37( si_helix_18_37->HardCopy());
//        util::ShPtr< assemble::SSE> sp_helix_63_81( si_helix_63_81->HardCopy());
//        util::ShPtr< assemble::SSE> sp_strand_2_13( si_strand_2_13->HardCopy());
//        util::ShPtr< assemble::SSE> sp_strand_40_45( si_strand_40_45->HardCopy());
//        util::ShPtr< assemble::SSE> sp_strand_52_61( si_strand_52_61->HardCopy());
//
//        sp_helix_18_37->SetToIdealConformationInPlace();
//        sp_strand_52_61->SetToIdealConformationInPlace();
//        sp_helix_63_81->SetToIdealConformationInPlace();
//        sp_strand_2_13->SetToIdealConformationInPlace();
//        sp_strand_40_45->SetToIdealConformationInPlace();
//
//        assemble::ProteinModel copy_model;
//
//        util::ShPtr< assemble::SSE> sse;
//        size_t nr_iters( 10);
//        // HELIX_HELIX
//        for( size_t i( 0); i < nr_iters; ++i)
//        {
//          util::GetLogger() << "HELIX_HELIX #" << i <<  std::endl;
//          sse = sp_helix_63_81.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_helix_18_37).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_18_37);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hha1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_helix_18_37.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_helix_63_81).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_63_81);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hha2_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_helix_63_81.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_helix_18_37).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_18_37);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hhb1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_helix_18_37.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_helix_63_81).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_63_81);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hhb2_" + util::Format()( i) + ".pdb"));
//        }
//
//        // HELIX_SHEET
//        for( size_t i( 0); i < nr_iters; ++i)
//        {
//          util::GetLogger() << "HELIX_SHEET #" << i <<  std::endl;
//          sse = sp_helix_18_37.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_strand_2_13).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_2_13);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hsa1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_helix_18_37.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_strand_40_45).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_40_45);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hsa2_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_helix_18_37.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_strand_2_13).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_2_13);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hsb1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_helix_18_37.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_strand_40_45).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_40_45);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_hsb2_" + util::Format()( i) + ".pdb"));
//        }
//
//        // SHEET_HELIX
//        for( size_t i( 0); i < nr_iters; ++i)
//        {
//          util::GetLogger() << "SHEET_HELIX #" << i <<  std::endl;
//          sse = sp_strand_40_45.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_helix_18_37).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_18_37);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_sha1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_strand_40_45.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_helix_63_81).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_40_45);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_sha2_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_strand_40_45.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_helix_18_37).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_18_37);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_shb1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_strand_52_61.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_helix_63_81).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_helix_63_81);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_shb2_" + util::Format()( i) + ".pdb"));
//        }
//
//        // SHEET_SHEET
//        for( size_t i( 0); i < nr_iters; ++i)
//        {
//          util::GetLogger() << "SHEET_SHEET #" << i <<  std::endl;
//          sse = sp_strand_40_45.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_strand_52_61).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_52_61);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_ssa1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_strand_52_61.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_side.Place( *sse, *sp_strand_40_45).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_40_45);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_ssa2_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_strand_40_45.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_strand_52_61).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_52_61);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_ssb1_" + util::Format()( i) + ".pdb"));
//
//          sse = sp_strand_52_61.HardCopy();
//          sse->SetToIdealConformationAtOrigin();
//          sse->Transform( place_top.Place( *sse, *sp_strand_40_45).First());
//          copy_model = model.HardCopyProteinModel();
//          copy_model.Insert( sp_strand_40_45);
//          copy_model.Insert( sse);
//          Proteins::WriteModelToPDB( copy_model, AddExampleOutputPathToFilename( place_side, "test_ssb2_" + util::Format()( i) + ".pdb"));
//        }
//      }

    /////////////////
    // helix_58_79 //
    /////////////////

      // initialize locator sse
      assemble::LocatorSSE sse_locator( 'A', 10, 17);

      // locate helix_58_79
      BCL_MessageStd( "Locating strand_10_17");
      util::SiPtr< const assemble::SSE> strand_10_17( sse_locator.Locate( model));
      BCL_ExampleIndirectAssert( strand_10_17.IsDefined(), true, "sse_locator.Locate( model)");

      // copy helix_58_79 and idealize
      util::ShPtr< assemble::SSE> strand_10_17_copy( strand_10_17->HardCopy());
      strand_10_17_copy->SetToIdealConformationInPlace();

      // insert helix_58_79 into empty_model
      initial_model.Insert( strand_10_17_copy);

      // initialize the placement
      BCL_MessageStd( "Creating the PlacementSSEShortLoop");
      fold::PlacementSSEShortLoop placement( 7, 0.5);

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write placement_sse_short_loop_before.pdb");
      Proteins::WriteModelToPDB
      (
        initial_model, AddExampleOutputPathToFilename( placement, "placement_sse_short_loop_before.pdb")
      );

    ////////////////
    // strand_1_7 //
    ////////////////

      // create mutated_model_a
      BCL_MessageStd( "Creating mutated_model_a")
      assemble::ProteinModel mutated_model_a( initial_model);

      // pick strand_1_7
      sse_locator.SetSSEID( 1, 7);
      BCL_MessageStd( "Locating strand_1_7");
      util::SiPtr< const assemble::SSE> strand_1_7
      (
        sse_locator.Locate( model)
      );
      BCL_Example_Check
      (
        strand_1_7.IsDefined(), "could not locate strand_1_7"
      );

      // call the placement
      BCL_MessageStd( "Calling the PlacementSSEShortLoop for strand_1_7");

      // insert strand_1_7 into mutated_model_a
      storage::Pair< math::TransformationMatrix3D, bool> transformation_1_7
      (
        placement.Place( *strand_1_7, mutated_model_a)
      );

      // assert the placement was correct
      BCL_Example_Check
      (
        transformation_1_7.Second(), "The boolean returned from placement strand_1_7 is false!"
      );

      // copy the selected sse, idealize it, transform it and insert
      util::ShPtr< assemble::SSE> strand_1_7_copy( strand_1_7->HardCopy());
      strand_1_7_copy->SetToIdealConformationAtOrigin();
      strand_1_7_copy->Transform( transformation_1_7.First());
      mutated_model_a.Insert( strand_1_7_copy);

      // write mutated protein model to an example pdb
      BCL_MessageStd( "placement_sse_short_loop_a.pdb");
      Proteins::WriteModelToPDB
      (
        mutated_model_a, AddExampleOutputPathToFilename( placement, "placement_sse_short_loop_a.pdb")
      );

    /////////////////
    // helix_23_34 //
    /////////////////

      // create mutated_model_a
      BCL_MessageStd( "Creating mutated_model_b")
      assemble::ProteinModel mutated_model_b( mutated_model_a);

      // pick helix_23_34
      sse_locator.SetSSEID( 23, 34);
      BCL_MessageStd( "Locating helix_23_34");
      util::SiPtr< const assemble::SSE> helix_23_34( sse_locator.Locate( model));
      BCL_ExampleIndirectAssert( helix_23_34.IsDefined(), true, "could not locate helix_23_34");

      // call the placement
      BCL_MessageStd( "Calling the PlacementSSEShortLoop for helix_23_34");

      storage::Pair< math::TransformationMatrix3D, bool> transformation_23_34
      (
        placement.Place( *helix_23_34, mutated_model_a)
      );

      // assert the placement was correct
      BCL_Example_Check
      (
        transformation_23_34.Second(), "The boolean returned from placement helix_23_34 is false!"
      );

      // copy the selected sse, idealize it, transform it and insert
      util::ShPtr< assemble::SSE> helix_23_34_copy( helix_23_34->HardCopy());
      helix_23_34_copy->SetToIdealConformationAtOrigin();
      helix_23_34_copy->Transform( transformation_23_34.First());
      mutated_model_b.Insert( helix_23_34_copy);

      // write mutated protein model to an example pdb
      BCL_MessageStd( "placement_sse_short_loop_b.pdb");
      Proteins::WriteModelToPDB
      (
        mutated_model_b, AddExampleOutputPathToFilename( placement, "placement_sse_short_loop_b.pdb")
      );

    //////////////////
    // strand_40_45 //
    //////////////////

      // create mutated_model_a
      BCL_MessageStd( "Creating mutated_model_c")
      assemble::ProteinModel mutated_model_c( mutated_model_b);

      // pick strand_40_45
      sse_locator.SetSSEID( 40, 45);
      BCL_MessageStd( "Locating strand_40_45");
      util::SiPtr< const assemble::SSE> strand_40_45
      (
        sse_locator.Locate( model)
      );
      BCL_Example_Check
      (
        strand_40_45.IsDefined(), "could not locate strand_40_45"
      );

      // call the placement
      BCL_MessageStd( "Calling the PlacementSSEShortLoop for strand_40_45");

      // insert strand_1_7 into mutated_model_a
      storage::Pair< math::TransformationMatrix3D, bool> transformation_40_45
      (
        placement.Place( *strand_40_45, mutated_model_c)
      );

      // assert the placement was correct
      BCL_Example_Check
      (
        transformation_1_7.Second(), "The boolean returned from placement strand_1_7 is false!"
      );

      // copy the selected sse, idealize it, transform it and insert
      util::ShPtr< assemble::SSE> strand_40_45_copy( strand_40_45->HardCopy());
      strand_40_45_copy->SetToIdealConformationAtOrigin();
      strand_40_45_copy->Transform( transformation_40_45.First());
      mutated_model_c.Insert( strand_40_45_copy);

      // write mutated protein model to an example pdb
      BCL_MessageStd( "placement_sse_short_loop_c.pdb");
      Proteins::WriteModelToPDB
      (
        mutated_model_c, AddExampleOutputPathToFilename( placement, "placement_sse_short_loop_c.pdb")
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementSSEShortLoop

  const ExampleClass::EnumType ExampleFoldPlacementSSEShortLoop::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPlacementSSEShortLoop())
  );

} // namespace bcl

