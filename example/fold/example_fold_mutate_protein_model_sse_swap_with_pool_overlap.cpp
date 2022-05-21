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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool_overlap.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_compare_type.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "find/bcl_find_collector_criteria_combined.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_swap_with_pool_overlap.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSESwapWithPoolOverlap :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSESwapWithPoolOverlap *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSESwapWithPoolOverlap( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      // build models from pdb sequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel
        protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // initialize chain
      util::ShPtr< assemble::Chain> this_chain( new assemble::Chain( protein_model.GetChains()( 0)->GetSequence()));

      // locate sses
      util::SiPtr< const assemble::SSE>  helix_23_34( assemble::LocatorSSE( 'A', 23, 34).Locate( protein_model));
      util::SiPtr< const assemble::SSE> strand_40_45( assemble::LocatorSSE( 'A', 40, 45).Locate( protein_model));
      util::SiPtr< const assemble::SSE> strand_48_50( assemble::LocatorSSE( 'A', 48, 50).Locate( protein_model));

      // create new SSE that is 42 to 50
      util::ShPtr< assemble::SSE> strand_42_50
      (
        new assemble::SSE( this_chain->GetSequence()->SubSequence( 41, 9), biol::GetSSTypes().STRAND)
      );
      strand_42_50->SetToIdealConformationAtOrigin();

      // initialize pool and read it from 1ubi.pool file
      BCL_MessageStd( "reading pool from 1ubi.pool");
      util::ShPtr< assemble::SSEPool> sse_pool( new assemble::SSEPool());
      const std::string pool_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pool"));
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, pool_filename);
      sse_pool->ReadSSEPool( read, protein_model, 0, 0);
      io::File::CloseClearFStream( read);
      // create the sse pool
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sse_pool);
      protein_model.SetProteinModelData( sp_model_data);

      // scheme
      const std::string scheme( "test_scheme");

      // initialize collector for use in constructor
      const util::ShPtr< find::CollectorCriteriaCombined< assemble::SSE> > sp_collector
      (
        new find::CollectorCriteriaCombined< assemble::SSE>
        (
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSE, assemble::SSE, bool> >
          (
            new assemble::SSECompareType()
          )
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "test default constructor");
      fold::MutateProteinModelSSESwapWithPoolOverlap swap_default;

      // test constructor
      BCL_MessageStd( "test constructor from pool and scheme");
      fold::MutateProteinModelSSESwapWithPoolOverlap swap( sp_collector, false, false, scheme);
      BCL_Example_Check
      (
        swap.GetScheme() == scheme,
        "The scheme for swap should be " + scheme + " not " + swap.GetScheme()
      );

      // test copy constructor
      BCL_MessageStd( "test Clone()");
      fold::MutateProteinModelSSESwapWithPoolOverlap swap_copy( swap);
      BCL_Example_Check
      (
        swap_copy.GetScheme() == swap.GetScheme(),
        "The scheme for swap_copy should be " + swap.GetScheme() + " not " + swap_copy.GetScheme()
      );

      // test Clone()
      BCL_MessageStd( "test Clone()");
      util::ShPtr< fold::MutateProteinModelSSESwapWithPoolOverlap> sp_swap( swap.Clone());
      BCL_Example_Check
      (
        sp_swap->GetScheme() == swap.GetScheme(),
        "The scheme for sp_swap should be " + swap.GetScheme() + " not " + sp_swap->GetScheme()
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( swap.GetClassIdentifier(), GetStaticClassName( swap));

      // test GetScheme()
      BCL_ExampleCheck( swap.GetScheme(), scheme);

    ////////////////
    // operations //
    ////////////////

      {
        BCL_MessageStd( "Testing operator() with model_a");
        // initialize new_model
        assemble::ProteinModel new_model_a( this_chain);
        new_model_a.SetProteinModelData( protein_model.GetProteinModelData());
        new_model_a.Insert( util::ShPtr< assemble::SSE>( helix_23_34->Clone()));

        // write mutated proteinmodel to an example pdb
        Proteins::WriteModelToPDB
        (
          new_model_a,
          AddExampleOutputPathToFilename( swap, "mutated_swap_with_pool_overlap_a_before.pdb")
        );

        // mutate the protein model
        BCL_MessageStd( "calling mutate on new_model_a");
        assemble::ProteinModel mutated_model_a( *swap( new_model_a).GetArgument());

        // write mutated proteinmodel to an example pdb
        Proteins::WriteModelToPDB
        (
          mutated_model_a,
          AddExampleOutputPathToFilename( swap, "mutated_swap_with_pool_overlap_a_after.pdb")
        );
      }

    ////////////////
    // operations //
    ////////////////

      {
        BCL_MessageStd( "Testing operator() with model_b");

        // initialize new_model_b
        BCL_MessageStd( "creating new_model_b");
        assemble::ProteinModel new_model_b( this_chain);
        new_model_b.SetProteinModelData( protein_model.GetProteinModelData());
        new_model_b.Insert( util::ShPtr< assemble::SSE>( strand_40_45->Clone()));
        new_model_b.Insert( util::ShPtr< assemble::SSE>( strand_48_50->Clone()));

        // write mutated proteinmodel to an example pdb
        Proteins::WriteModelToPDB
        (
          new_model_b,
          AddExampleOutputPathToFilename( swap, "mutated_swap_with_pool_overlap_b_before.pdb")
        );

        // whether it pick strand_40_45 or strand_48_50, the sse from pool will strand_44_50 that overlaps with both
        BCL_MessageStd( "calling mutate on new_model_b");
        assemble::ProteinModel mutated_model_b( *swap( new_model_b).GetArgument());

        // write mutated proteinmodel to an example pdb
        Proteins::WriteModelToPDB
        (
          mutated_model_b,
          AddExampleOutputPathToFilename( swap, "mutated_swap_with_pool_overlap_b_after.pdb")
        );

        // the resultant pdb should have only one strand which is strand_44_50 therefore assert that
        BCL_Example_Check
        (
          mutated_model_b.GetSSEs().GetSize() == 1 &&
          mutated_model_b.GetSSEs()( 0)->GetFirstAA()->GetSeqID() == 44 &&
          mutated_model_b.GetSSEs()( 0)->GetLastAA()->GetSeqID() == 50,
          "The replacement did not work correctly, there should be only one sse left which is strand 44 to 50"
        );
      }

    ////////////////
    // operations //
    ////////////////

      {
        BCL_MessageStd( "Testing operator() with model_c");

        // initialize new_model_b
        BCL_MessageStd( "creating new_model_c");
        assemble::ProteinModel new_model_c( this_chain);
        new_model_c.SetProteinModelData( protein_model.GetProteinModelData());
        new_model_c.Insert( util::ShPtr< assemble::SSE>( strand_42_50));

        // write mutated proteinmodel to an example pdb
        Proteins::WriteModelToPDB
        (
          new_model_c,
          AddExampleOutputPathToFilename( swap, "mutated_swap_with_pool_overlap_c_before.pdb")
        );

        // it should remove 42_50 and insert both 40_45 and 48_50 at the same tim
        BCL_MessageStd( "calling mutate on new_model_c");
        assemble::ProteinModel mutated_model_c( *swap( new_model_c).GetArgument());

        // write mutated proteinmodel to an example pdb
        Proteins::WriteModelToPDB
        (
          mutated_model_c,
          AddExampleOutputPathToFilename( swap, "mutated_swap_with_pool_overlap_c_after.pdb")
        );

        // the resultant pdb should have two strands 40 to 45 and 48 50
        BCL_Example_Check
        (
          mutated_model_c.GetSSEs().GetSize() == 2 &&
          mutated_model_c.GetSSEs()( 0)->GetFirstAA()->GetSeqID() == 40 &&
          mutated_model_c.GetSSEs()( 0)->GetLastAA()->GetSeqID() == 45 &&
          mutated_model_c.GetSSEs()( 1)->GetFirstAA()->GetSeqID() == 48 &&
          mutated_model_c.GetSSEs()( 1)->GetLastAA()->GetSeqID() == 50,
          "The replacement did not work correctly, there should be two SSEs, 40-45 and 48_50"
        );

        // create mutate w/ bending
        fold::MutateProteinModelSSESwapWithPoolOverlap swap_bend( sp_collector, false, true, scheme);

        // mutate the protein model (random sse will be moved)
        assemble::ProteinModel bent_sse_model( *swap_bend( new_model_c).GetArgument());
        Proteins::WriteModelToPDB
        (
          bent_sse_model,
          AddExampleOutputPathToFilename( swap_bend, "mutated_swap_with_pool_overlap_bend.pdb")
        );
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSESwapWithPoolOverlap

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSESwapWithPoolOverlap::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSESwapWithPoolOverlap())
  );

} // namespace bcl
