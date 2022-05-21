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
#include "fold/bcl_fold_placement_sse_into_body.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_pick_body_random.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "restraint/bcl_restraint_body.h"
#include "restraint/bcl_restraint_contains_body_origin.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_placement_sse_into_body.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFoldPlacementSSEIntoBody :
    public ExampleInterface
  {
  public:

    ExampleFoldPlacementSSEIntoBody *Clone() const
    { return new ExampleFoldPlacementSSEIntoBody( *this);}

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
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));
      BCL_MessageStd( "protein model built");

      // get the strand secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> strand_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().STRAND)
      );

      // create a ShPtrVector of bodies which will be created from "strand_sses"
      util::ShPtrVector< assemble::SSEGeometryInterface> bodies;

      // add the coord::Bodies of "strand_sses" to "bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( strand_sses.Begin()), sse_itr_end( strand_sses.End());
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
            new restraint::ContainsBodyOrigin() //initialize body with ContainsBodyOrigin occupance check object
          )
        )
      );

      // create ShPtr to a math::MutateInterface< math::TransformationMatrix3D> "mutate_tm_3d_null"
      util::ShPtr< math::MutateInterface< math::TransformationMatrix3D> > mutate_tm_3d_null
      (
        // initialize with MutateTransformationMatrix3DNull to place SSE into body with same orientation as body
        new restraint::MutateTransformationMatrix3DNull()
      );

      // create PlacementSSEIntoBody "placement" and initialize with "body_restraint", "pick", and "mutate_tm_3d_null"
      fold::PlacementSSEIntoBody placement
      (
        body_restraint,
        pick,
        mutate_tm_3d_null
      );

      // get the helix secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> helix_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().HELIX)
      );

      // create copy of "protein_model" named "place_model" which is the model the helix will be placed in
      assemble::ProteinModel place_model( protein_model);

      // remove all the strands from "place_model" so that the restraint bodies aren't all occupied and an sse can be
      // added to one of the randomly selected restraint_bodies
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator itr( strand_sses.Begin()), itr_end( strand_sses.End());
        itr != itr_end;
        ++itr
      )
      {
        // remove the strand currently denoted by "itr"
        place_model.Remove( **itr);
      }

      // make sure that "place_model" only has two SSEs left (the helices)
      BCL_Example_Check
      (
        place_model.GetSSEs().GetSize() == 2, "place model should only have 2 sses but has "
        + util::Format()( place_model.GetSSEs().GetSize())
      );

      // create storage::Pair "transformation" and initialize with the transformation matrix which will be applied
      // to a helix of "helix_sses"
      storage::Pair< math::TransformationMatrix3D, bool> transformation
      (
        placement.Place( **helix_sses.Begin(), place_model)
      );

      // make sure "transformation" indicates the placement was successful
      BCL_Example_Check
      (
        transformation.Second(), "placement returned false"
      );

      // make a copy of the protein model "protein_model" into which an SSE will be placed
      assemble::ProteinModel new_model( protein_model);

      BCL_MessageStd( "new_model protein model created");

      // create a copy of the selected sse
      util::ShPtr< assemble::SSE> copy_sse( ( **helix_sses.Begin()).Clone());

      // get the transformation matrix of "copy_sse" at the origin
      math::TransformationMatrix3D origin( math::Inverse( copy_sse->GetOrientation()));

      // apply the "transformation" to "origin" so that "origin" has the rotational and translational
      // information of "transformation"
      origin( transformation.First());

      // move and rotate "copy_sse" according to "origin"
      copy_sse->Transform( origin);

      // replace the sse with "copy_sse"
      new_model.Replace( copy_sse);

      // write mutated proteinmodel to an example pdb
      BCL_MessageStd( "write placement_sse_into_body.pdb");
      Proteins::WriteModelToPDB( new_model, AddExampleOutputPathToFilename( placement, "placement_sse_into_body.pdb"));

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFoldPlacementSSEIntoBody

  const ExampleClass::EnumType ExampleFoldPlacementSSEIntoBody::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldPlacementSSEIntoBody())
  );

} // namespace bcl
