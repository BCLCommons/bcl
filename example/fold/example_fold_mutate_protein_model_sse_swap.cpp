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
#include "fold/bcl_fold_mutate_protein_model_sse_swap.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_sse_compare_type.h"
#include "find/bcl_find_collector_criteria_combined.h"
#include "find/bcl_find_locator_criteria.h"
#include "find/bcl_find_pick_criteria_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_swap.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSESwap :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSESwap *Clone() const
    { return new ExampleFoldMutateProteinModelSSESwap( *this);}

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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //build models from pdbsequences
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 5;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

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
      fold::MutateProteinModelSSESwap mutate_swap_sses( locator_ss_type);

      // mutate the protein model (random sse will be moved)
      assemble::ProteinModel mutated_model( *mutate_swap_sses( protein_model).GetArgument());

      const double rmsd_start( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, protein_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      const double rmsd_mutate( assemble::Quality::Calculate( quality::GetMeasures().e_RMSD, mutated_model, biol::GetAtomTypes().GetBackBoneAtomTypes()));
      BCL_MessageStd( "rmsd before swapping sses " + util::Format()( rmsd_start));
      BCL_MessageStd( "rmsd after swapping sses " + util::Format()( rmsd_mutate));
      BCL_ExampleCheck( protein_model.GetSSEs().GetSize(), mutated_model.GetSSEs().GetSize());

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write mutated_swap.pdb");
      Proteins::WriteModelToPDB( mutated_model, AddExampleOutputPathToFilename( mutate_swap_sses, "mutated_swap.pdb"));

      // create mutate w/ bending
      fold::MutateProteinModelSSESwap mutate_swap_bend( locator_ss_type, true);

      // mutate the protein model (random sse will be moved)
      assemble::ProteinModel bent_sse_model( *mutate_swap_bend( protein_model).GetArgument());
      Proteins::WriteModelToPDB
      (
        bent_sse_model,
        AddExampleOutputPathToFilename( mutate_swap_bend, "mutated_swap_bend.pdb")
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldMutateProteinModelSSESwap

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSESwap::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSESwap())
  );

} // namespace bcl
