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
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_locator_sse_furthest.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "find/bcl_find_collector_criteria_wrapper.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "fold/bcl_fold_placement_sse_next_to_sse.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_add.cpp
  //!
  //! @author karakam
  //! @date Dec 17, 2009
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSEAdd :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSEAdd *Clone() const
    {
      return new ExampleFoldMutateProteinModelSSEAdd( *this);}

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

      // read in the protein model 1J27
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1J27A.pdb"));
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 7;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;
      assemble::ProteinModel
        protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      // create the sse pool
      BCL_MessageStd( "Creating the SSEPool");
      const util::ShPtr< assemble::SSEPool> sp_sse_pool( new assemble::SSEPool( protein_model.GetSSEs()));
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData());
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sp_sse_pool);
      protein_model.SetProteinModelData( sp_model_data);

      // initialize the sse_pool_picker
      BCL_MessageStd( "Creating the PickSSERandom");
      util::ShPtr
      <
        find::PickCriteriaInterface
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        >
      >
      sse_pool_picker
      (
        new find::PickCriteriaWrapper
        <
          util::SiPtr< const assemble::SSE>, util::SiPtrList< const assemble::SSE>, assemble::DomainInterface
        >
        (
          assemble::PickSSERandom()
        )
      );

      // initialize collector criteria wrapper for collector sse
      util::ShPtr
      <
        find::CollectorCriteriaInterface
        <
          util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::SSE
        >
      >
      collector_sse_wrap
      (
        new find::CollectorCriteriaWrapper
        <
          util::SiPtrList< const assemble::SSE>, assemble::DomainInterface, assemble::SSE
        >
        (
          assemble::CollectorSSE()
        )
      );

      // initialize the locator to locate a sse within the protein model to which next to a second sse is placed
      BCL_MessageStd( "Creating the SSEFurthestLocator");

      util::ShPtr
      <
        find::LocatorCriteriaInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface, assemble::SSE>
      >
      sse_furthest
      (
        new assemble::LocatorSSEFurthest< assemble::SSE>
        (
          collector_sse_wrap
        )
      );

      // initialize the placement
      BCL_MessageStd( "Creating the PlacementSSENextToSSE");
      util::ShPtr< fold::PlacementInterface< assemble::SSE, assemble::ProteinModel> >
        placement( new fold::PlacementSSENextToSSE( sse_furthest));

      // create mutate object that places sse picked from pool into a given protein model
      BCL_MessageStd( "Creating the MutateProteinModelSSEAdd");
      fold::MutateProteinModelSSEAdd mutate( sse_pool_picker, placement);

      // create an empty protein model
      assemble::ProteinModel empty_model;
      empty_model.SetProteinModelData( protein_model.GetProteinModelData());

      // Initialize chain information into the empty model
      BCL_MessageStd( "Initializing chain information for empty model");
      empty_model.Insert
      (
        util::ShPtr< assemble::Chain>
        (
          new assemble::Chain
          (
            protein_model.GetChains()( 0)->GetSequence()
          )
        )
      );

      {
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator sse_itr
        (
          protein_model.GetChains()( 0)->GetData().Begin()
        );
        empty_model.Insert( *( sse_itr++));
        empty_model.Insert( *( sse_itr++));
        empty_model.Insert( *( sse_itr));
      }

      // write empty protein model to an example pdb
      BCL_MessageStd( "write empty_model.pdb");
      Proteins::WriteModelToPDB( empty_model, AddExampleOutputPathToFilename( mutate, "empty_model.pdb"));

      assemble::ProteinModel mutated_model;
      // mutate the protein model till it is successful
      for( size_t counter( 0); counter < 10; ++counter)
      {
        BCL_MessageStd( "Calling the mutate");
        math::MutateResult< assemble::ProteinModel> mutate_result( mutate( empty_model));
        // if the mutate failed continue
        if( !mutate_result.GetArgument().IsDefined())
        {
          continue;
        }
        mutated_model = *mutate_result.GetArgument();
        if( mutated_model.GetSSEs().GetSize() > 3)
        {
          break;
        }
      }

      // write mutated protein model to an example pdb
      BCL_MessageStd( "write add_mutated_model.pdb");
      Proteins::WriteModelToPDB( mutated_model, AddExampleOutputPathToFilename( mutate, "add_mutated_model.pdb"));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSEAdd

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSEAdd::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSEAdd())
  );

} // namespace bcl
