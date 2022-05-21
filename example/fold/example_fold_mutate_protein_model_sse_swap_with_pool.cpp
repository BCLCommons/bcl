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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_with_pool.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_compare_type.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "find/bcl_find_collector_criteria_combined.h"
#include "find/bcl_find_locator_criteria.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_swap_with_pool.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSESwapWithPool :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSESwapWithPool *Clone() const
    { return new ExampleFoldMutateProteinModelSSESwapWithPool( *this);}

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
      // build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // initialize pool and read it from 1ubi.pool file
      io::IFStream read;
      BCL_MessageStd( "reading pool from 1ubi.pool");
      util::ShPtr< assemble::SSEPool> sse_pool( new assemble::SSEPool());
      const std::string pool_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pool"));
      BCL_ExampleMustOpenInputFile( read, pool_filename);
      sse_pool->ReadSSEPool( read, protein_model, 0, 0);
      io::File::CloseClearFStream( read);

      // create the sse pool
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sse_pool);
      protein_model.SetProteinModelData( sp_model_data);

      // initialize collector for use in locator
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

      // initialize random picker
      const util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
      > sp_picker
      (
        new find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::SSE
        >
        (
          // set up picker with random picker
          assemble::PickSSERandom()
        )
      );

      // create locator to locate sses of same type in pool
      const find::LocatorCriteria
      <
        util::SiPtr< const assemble::SSE>,
        util::SiPtrList< const assemble::SSE>,
        assemble::SSE,
        util::SiPtrList< const assemble::SSE>
      > locator_ss_type( sp_collector, sp_picker);

      // create mutate object "mutate_random" from max translation and rotation and SSELocatorRandom
      BCL_MessageStd( "creating MutateProteinModelSSESwapWithPool");
      fold::MutateProteinModelSSESwapWithPool mutate_swap_sses_with_pool( locator_ss_type);

      // initialize new_model
      BCL_MessageStd( "creating new_model_a");
      util::ShPtr< assemble::Chain> this_chain( new assemble::Chain( protein_model.GetChains()( 0)->GetSequence()));

      assemble::ProteinModel new_model( this_chain);
      new_model.SetProteinModelData( protein_model.GetProteinModelData());

      // initialize sse locator
      BCL_MessageStd( "creating LocatorSSE");
      assemble::LocatorSSE sse_locator( protein_model.GetChains()( 0)->GetChainID(), 23, 34);

      // initialize located sse
      BCL_MessageStd( "Locating helix_23_34");
      util::SiPtr< const assemble::SSE> helix_23_34( sse_locator.Locate( protein_model));

      // assert the located sse is defined
      BCL_Example_Check
      (
        helix_23_34.IsDefined(), "helix_23_34 was not correctly located!"
      );

      // insert helix_23_34
      BCL_MessageStd( "Inserting helix_23_34");
      new_model.Insert( util::ShPtr< assemble::SSE>( helix_23_34->Clone()));

      // mutate the protein model (random sse will be moved
      BCL_MessageStd( "calling mutate on new_model");
      assemble::ProteinModel mutated_model( *mutate_swap_sses_with_pool( new_model).GetArgument());

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write mutated_swap_with_pool_before.pdb");
      Proteins::WriteModelToPDB
      (
        new_model, AddExampleOutputPathToFilename( mutate_swap_sses_with_pool, "mutated_swap_with_pool_before.pdb")
      );

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write mutated_swap_with_pool_after.pdb");
      Proteins::WriteModelToPDB
      (
        mutated_model, AddExampleOutputPathToFilename( mutate_swap_sses_with_pool, "mutated_swap_with_pool_after.pdb")
      );

      // the resultant pdb should have only one SSE
      BCL_Example_Check
      (
        mutated_model.GetSSEs().GetSize() == 1,
        "The replacement did not work correctly, there should be only one sse in the protein model"
      );

      BCL_Example_Check
      (
        !sse_locator.Locate( mutated_model).IsDefined(),
        "The replacement did not work correctly, helix_23_34 should have been replaced with another sse from pool"
      );

      // create mutate w/ bending
      fold::MutateProteinModelSSESwapWithPool mutate_swap_pool_bend( locator_ss_type, true);

      // mutate the protein model (random sse will be moved)
      assemble::ProteinModel bent_sse_model( *mutate_swap_pool_bend( new_model).GetArgument());
      Proteins::WriteModelToPDB
      (
        bent_sse_model,
        AddExampleOutputPathToFilename( mutate_swap_pool_bend, "mutated_swap_pool_bend.pdb")
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSESwapWithPool

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSESwapWithPool::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSESwapWithPool())
  );

} // namespace bcl
