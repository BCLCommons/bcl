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
#include "fold/bcl_fold_mutate_protein_model_sse_swap_body.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_pick_body_random.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "math/bcl_math_mutate_result.h"
#include "restraint/bcl_restraint_contains_body_origin.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_protein_model_sse_swap_body.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldMutateProteinModelSSESwapBody :
    public ExampleInterface
  {
  public:

    ExampleFoldMutateProteinModelSSESwapBody *Clone() const
    { return new ExampleFoldMutateProteinModelSSESwapBody( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));

      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));
      BCL_MessageStd( "protein model built");

      // get the strand secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> strand_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().STRAND)
      );

      // create a ShPtrVector of bodies which will be created from "strand_sses"
      util::ShPtrVector< assemble::SSEGeometryInterface> bodies;

      // add one of the coord::Bodies of "strand_sses" to "bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( strand_sses.Begin()), sse_itr_end( sse_itr + 1);
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> temp_sse( ( *sse_itr)->Clone());
        bodies.PushBack( temp_sse);
      }

      // create ShPtr to a assemble::PickCriteriaInterface "pick" which is used to determine how a restraint body is
      // selected
      util::ShPtr
      <
        find::PickCriteriaInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface>
      > pick
      (
        // initialize with a PickCriteriaWrapper initialized with a PickBody Random to randomly select bodies
        // to place SSEs in
        new find::PickCriteriaWrapper< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface>
        (
          util::ShPtr< find::PickInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface> > >
          (
            new find::PickBodyRandom()
          )
        )
      );

      // create ShPtr to a restraint::Body "body_restraint"
      util::ShPtr< restraint::Body> body_restraint
      (
        new restraint::Body
        (
          util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> >( bodies.Clone()), //< initialize Body with clone of "bodies"
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          (
            new restraint::ContainsBodyOrigin() //<initialize body with ContainsBodyOrigin occupancy check object
          )
        )
      );

      // create ShPtr to a LocatorInterface which will define how the sse to be swapped is located
      // initialize with
      util::ShPtr< find::LocatorInterface< util::SiPtr< const assemble::SSE>, assemble::DomainInterface> > sse_locator
      (
        new assemble::LocatorSSE( 'A', 23, 34)
      );

      // create PlacementSSEIntoBody "placement" and initialize with "body_restraint", "pick", and "mutate_tm_3d_null"
      fold::MutateProteinModelSSESwapBody mutate_swap
      (
        body_restraint,
        pick,
        sse_locator
      );

      // get the helix secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> helix_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().HELIX)
      );

      // create copy of "protein_model" named "swap_model" which is the model where a helix will be swapped to a new
      // body
      assemble::ProteinModel swap_model( protein_model);

      // remove all the strands from "swap_model" so that the restraint bodies aren't all occupied and an SSE can be
      // swapped to one of the randomly selected restraint_bodies
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator itr( strand_sses.Begin()), itr_end( strand_sses.End());
        itr != itr_end;
        ++itr
      )
      {
        // remove the strand currently denoted by "itr"
        swap_model.Remove( **itr);
      }

      // create MutateResult "mutate_result" and initialize with result of "mutate_swap"
      math::MutateResult< assemble::ProteinModel> mutate_result( mutate_swap( swap_model));

      // make sure the swap worked
      const linal::Vector3D correct_coord( 27.4796, 32.4955, 9.6245);
      const linal::Vector3D swapped_coord( sse_locator->Locate( *mutate_result.GetArgument())->GetCenter());
      BCL_ExampleIndirectCheckWithinTolerance( correct_coord, swapped_coord, 0.001, "SSE Swap");

      // write mutated proteinmodel to an example pdb
      const std::string
        out_filename( AddExampleOutputPathToFilename( mutate_swap, "mutate_protein_model_sse_swap_body.pdb"));
      BCL_MessageStd( "write mutate_protein_model_sse_swap_body.pdb");
      Proteins::WriteModelToPDB( *mutate_result.GetArgument(), out_filename);

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementSSEIntoBody

  const ExampleClass::EnumType ExampleFoldMutateProteinModelSSESwapBody::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldMutateProteinModelSSESwapBody())
  );

} // namespace bcl
