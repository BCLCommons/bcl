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
#include "fold/bcl_fold_placement_sse_next_to_sse.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse_random.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "find/bcl_find_locator_criteria_wrapper.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_placement_sse_next_to_sse.cpp
  //!
  //! @author karakam
  //! @date Dec 17, 2009
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPlacementSSENextToSSE :
    public ExampleInterface
  {
  public:

    ExampleFoldPlacementSSENextToSSE *Clone() const
    {
      return new ExampleFoldPlacementSSENextToSSE( *this);
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
//        fold::PlacementSSENextToSSE place_default;
//
//        // initialize the pick criteria for SSERandom
//        find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::ProteinModel, assemble::SSE>
//          locator_wrap( ( assemble::LocatorSSERandom()));
//
//        // construct from a locator
//        fold::PlacementSSENextToSSE placement( locator_wrap);
//
//        {
//          // create an empty model
//          assemble::ProteinModel model;
//          model.Insert( util::ShPtr< assemble::Chain>( new assemble::Chain( protein_model.GetChain( 'A')->GetSequence())));
//
//          // make copies of SSEs
//          util::ShPtr< assemble::SSE> sp_helix_18_37( si_helix_18_37->HardCopy());
//          util::ShPtr< assemble::SSE> sp_helix_63_81( si_helix_63_81->HardCopy());
//          util::ShPtr< assemble::SSE> sp_strand_2_13( si_strand_2_13->HardCopy());
//          util::ShPtr< assemble::SSE> sp_strand_40_45( si_strand_40_45->HardCopy());
//          util::ShPtr< assemble::SSE> sp_strand_52_61( si_strand_52_61->HardCopy());
//
//          sp_helix_18_37->SetToIdealConformationInPlace();
//          sp_strand_52_61->SetToIdealConformationInPlace();
//
//          sp_helix_63_81->SetToIdealConformationAtOrigin();
//          sp_strand_2_13->SetToIdealConformationAtOrigin();
//          sp_strand_40_45->SetToIdealConformationAtOrigin();
//
//          // insert a single helix
//          model.Insert( sp_helix_18_37);
//
//          // make a copy of the model
//          assemble::ProteinModel copy_model;
//
//          size_t nr_iters( 30);
//
//          // HELIX_HELIX
//          for( size_t i( 0); i < nr_iters; ++i)
//          {
//            util::GetLogger() <<"HELIX_HELIX #" << i << '\n';
//            util::ShPtr< assemble::SSE> new_helix( sp_helix_63_81->HardCopy());
//
//            // try to place
//            storage::Pair< math::TransformationMatrix3D, bool> place_pair( placement.Place( *new_helix, model));
//
//            // make sure it was successful
//            BCL_ExampleIndirectCheck( place_pair.Second(), true, "placement failed!!");
//
//            // print and apply the transformation
//            util::GetLogger() << "placement\n" << place_pair.First() << '\n';
//            new_helix->Transform( place_pair.First());
//
//            // make copy of the model and insert SSE
//            copy_model = model.HardCopyProteinModel();
//            copy_model.Insert( new_helix);
//
//            // write pdb
//            Proteins::WriteModelToPDB
//            (
//              copy_model, AddExampleOutputPathToFilename( placement, "place_test_helhel_" + util::Format()( i) + ".pdb")
//            );
//          }
//
//          // HELIX_SHEET
//          for( size_t i( 0); i < nr_iters; ++i)
//          {
//            util::GetLogger() <<"SHEET_HELIX#" << i << '\n';
//            util::ShPtr< assemble::SSE> new_strand( sp_strand_40_45->HardCopy());
//
//            // try to place
//            storage::Pair< math::TransformationMatrix3D, bool> place_pair( placement.Place( *new_strand, model));
//
//            // make sure it was successful
//            BCL_ExampleIndirectCheck( place_pair.Second(), true, "placement failed!!");
//
//            // print and apply the transformation
//            util::GetLogger() << "placement\n" << place_pair.First() << '\n';
//            new_strand->Transform( place_pair.First());
//
//            // make copy of the model and insert SSE
//            copy_model = model.HardCopyProteinModel();
//            copy_model.Insert( new_strand);
//
//            // write pdb
//            Proteins::WriteModelToPDB
//            (
//              copy_model, AddExampleOutputPathToFilename( placement, "place_test_shthel_" + util::Format()( i) + ".pdb")
//            );
//          }
//
//          // remove the helix, add strand
//          model.Remove( *sp_helix_18_37);
//          model.Insert( sp_strand_52_61);
//
//          // SHEET_HELIX
//          for( size_t i( 0); i < nr_iters; ++i)
//          {
//            util::GetLogger() <<"HELIX_SHEET#" << i << '\n';
//            util::ShPtr< assemble::SSE> new_helix( sp_helix_63_81->HardCopy());
//
//            // try to place
//            storage::Pair< math::TransformationMatrix3D, bool> place_pair( placement.Place( *new_helix, model));
//
//            // make sure it was successful
//            BCL_ExampleIndirectCheck( place_pair.Second(), true, "placement failed!!");
//
//            // print and apply the transformation
//            util::GetLogger() << "placement\n" << place_pair.First() << '\n';
//            new_helix->Transform( place_pair.First());
//
//            // make copy of the model and insert SSE
//            copy_model = model.HardCopyProteinModel();
//            copy_model.Insert( new_helix);
//
//            // write pdb
//            Proteins::WriteModelToPDB
//            (
//              copy_model, AddExampleOutputPathToFilename( placement, "place_test_helsht_" + util::Format()( i) + ".pdb")
//            );
//          }
//
//          // SHEET_SHEET
//          for( size_t i( 0); i < nr_iters; ++i)
//          {
//            util::GetLogger() <<"SHEET_SHEET#" << i << '\n';
//            util::ShPtr< assemble::SSE> new_strand( sp_strand_40_45->HardCopy());
//
//            // try to place
//            storage::Pair< math::TransformationMatrix3D, bool> place_pair( placement.Place( *new_strand, model));
//
//            // make sure it was successful
//            BCL_ExampleIndirectCheck( place_pair.Second(), true, "placement failed!!");
//
//            // print and apply the transformation
//            util::GetLogger() << "placement\n" << place_pair.First() << '\n';
//            new_strand->Transform( place_pair.First());
//
//            // make copy of the model and insert SSE
//            copy_model = model.HardCopyProteinModel();
//            copy_model.Insert( new_strand);
//
//            // write pdb
//            Proteins::WriteModelToPDB
//            (
//              copy_model, AddExampleOutputPathToFilename( placement, "place_test_shtsht_" + util::Format()( i) + ".pdb")
//            );
//          }
//
//
//        }
//      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
      //build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      assemble::SSEPool
      sse_pool
      (
        protein_model.GetSSEs()
      );

      // create an empty protein model
      BCL_MessageStd( "Initializing chain information for empty model");
      assemble::ProteinModel empty_model( protein_model.GetEmptyChains());

      // pick a random sse from the pool
      util::SiPtr< const assemble::SSE> sse_from_pool
      (
        assemble::PickSSERandom().Pick( sse_pool.GetNonOverlappingSSEs( empty_model))
      );

      BCL_ExampleIndirectAssert
      (
        sse_from_pool.IsDefined(),
        true,
        "assemble::PickSSERandom().Pick( sse_pool.GetNonOverlappingSSEs( empty_model))"
      );

      // copy the sse from pool
      util::ShPtr< assemble::SSE> first_sse( sse_from_pool->Clone());

      // insert into empty_model
      empty_model.Insert( first_sse);

      // initialize the pick criteria for SSERandom
      find::LocatorCriteriaWrapper< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>
      locator_wrap
      (
        ( assemble::LocatorSSERandom())
      );

      // initialize the placement
      BCL_MessageStd( "Creating the PlacementSSENextToSSE");
      fold::PlacementSSENextToSSE placement( locator_wrap);

      // initialize mutated models
      assemble::ProteinModel mutated_model, mutated_model_2, mutated_model_3;

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write placement_sse_next_to_sse_before.pdb");
      Proteins::WriteModelToPDB( empty_model, AddExampleOutputPathToFilename( placement, "placement_sse_next_to_sse_before.pdb"));

      // mutate the protein model till it is for adding the 2nd sse
      size_t counter( 0);
      while( true)
      {
        BCL_MessageStd( "Calling the Placement");
        mutated_model = empty_model;

        // collect nonoverlapping sses from the pool
        util::SiPtrList< const assemble::SSE> eligible_sses( sse_pool.GetNonOverlappingSSEs( mutated_model));

        // pick an sse
        util::SiPtr< const assemble::SSE> picked_sse( assemble::PickSSERandom().Pick( eligible_sses));

        // assert that picked sse is defined
        BCL_Example_Check
        (
          picked_sse.IsDefined(), "The picked sse is undefined!"
        );

        // do the placement
        storage::Pair< math::TransformationMatrix3D, bool> transformation
        (
          placement.Place( *picked_sse, empty_model)
        );

        // if the placement was succesfull
        if( transformation.Second())
        {
          // copy the sse from pool
          util::ShPtr< assemble::SSE> this_sse( picked_sse->Clone());

          // apply the transformation
          this_sse->Transform( transformation.First());

          // insert into empty_model
          mutated_model.Insert( this_sse);
        }
        else
        {
          BCL_MessageStd( "The boolean returned from placement is false!");
        }

        if( ++counter > 10 || mutated_model.GetSSEs().GetSize() > 1)
        {
          // write mutated protein model to an example pdb
          BCL_MessageStd( "placement_sse_next_to_sse_after_1.pdb");
          Proteins::WriteModelToPDB
          (
            mutated_model, AddExampleOutputPathToFilename( placement, "placement_sse_next_to_sse_after_1.pdb")
          );
          break;
        }
      }

      // mutate the protein model till it adds 2 more sses
      counter = 0;
      while( true)
      {
        BCL_MessageStd( "Calling the Placement");
        mutated_model_2 = mutated_model;

        // collect non-overlapping sses from the pool
        util::SiPtrList< const assemble::SSE> eligible_sses( sse_pool.GetNonOverlappingSSEs( mutated_model_2));

        // pick an sse
        util::SiPtr< const assemble::SSE> picked_sse( assemble::PickSSERandom().Pick( eligible_sses));

        // do the placement
        storage::Pair< math::TransformationMatrix3D, bool> transformation
        (
          placement.Place( *picked_sse, mutated_model)
        );

        // if the placement was successful
        if( transformation.Second())
        {
          // copy the sse from pool
          util::ShPtr< assemble::SSE> this_sse( picked_sse->Clone());

          // apply the transformation
          this_sse->Transform( transformation.First());

          // insert into empty_model
          mutated_model_2.Insert( this_sse);
        }
        else
        {
          BCL_MessageStd( "The boolean returned from placement is false!");
        }

        if( ++counter > 20 || mutated_model_2.GetSSEs().GetSize() > 2)
        {
          // write mutated protein model to an example pdb
          BCL_MessageStd( "placement_sse_next_to_sse_after_2.pdb");
          Proteins::WriteModelToPDB
          (
            mutated_model_2, AddExampleOutputPathToFilename( placement, "placement_sse_next_to_sse_after_2.pdb")
          );
          break;
        }
      }

      // mutate the protein model till it adds 2 more sses
      counter = 0;
      while( true)
      {
        BCL_MessageStd( "Calling the mutate");
        mutated_model_3 = mutated_model_2;

        // collect nonoverlapping sses from the pool
        util::SiPtrList< const assemble::SSE> eligible_sses( sse_pool.GetNonOverlappingSSEs( mutated_model_3));

        // pick an sse
        util::SiPtr< const assemble::SSE> picked_sse( assemble::PickSSERandom().Pick( eligible_sses));

        // do the placement
        storage::Pair< math::TransformationMatrix3D, bool> transformation
        (
          placement.Place( *picked_sse, mutated_model_2)
        );

        // if the placement was succesfull
        if( transformation.Second())
        {
          // copy the sse from pool
          util::ShPtr< assemble::SSE> this_sse( picked_sse->Clone());

          // apply the transformation
          this_sse->Transform( transformation.First());

          // insert into empty_model
          mutated_model_3.Insert( this_sse);
        }
        else
        {
          BCL_MessageStd( "The boolean returned from placement is false!");
        }

        if( ++counter > 20 || mutated_model_3.GetSSEs().GetSize() > 3)
        {
          // write mutated protein model to an example pdb
          BCL_MessageStd( "placement_sse_next_to_sse_after_3.pdb");
          Proteins::WriteModelToPDB
          (
            mutated_model_3, AddExampleOutputPathToFilename( placement, "placement_sse_next_to_sse_after_3.pdb")
          );
          break;
        }
      }

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementSSENextToSSE

  const ExampleClass::EnumType ExampleFoldPlacementSSENextToSSE::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPlacementSSENextToSSE())
  );

} // namespace bcl

