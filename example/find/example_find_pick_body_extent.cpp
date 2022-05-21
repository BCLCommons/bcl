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
#include "find/bcl_find_pick_body_extent.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_find_pick_body_extent.cpp
  //!
  //! @author linders
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFindPickBodyExtent :
    public ExampleInterface
  {
  public:

    ExampleFindPickBodyExtent *Clone() const
    { return new ExampleFindPickBodyExtent( *this);}

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

    /////////////////
    // preparation //
    /////////////////

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"));
      // create ProteinModel "protein_model" from "pdb"
      BCL_MessageStd( "building models from pdb chains and sse information");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 0;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 0;
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      // get the helix secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> helix_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().HELIX)
      );

      // get all the secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> all_sses
      (
        protein_model.GetChains()( 0)->GetSSEs()
      );

      // create a ShPtrVector of helix_bodies which will be created from "helix_sses"
      util::ShPtrVector< assemble::SSEGeometryInterface> helix_bodies;

      // add the coord::Bodies of "helix_sses" to "helix_bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( helix_sses.Begin()), sse_itr_end( helix_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> temp_sse( ( *sse_itr)->Clone());
        helix_bodies.PushBack( temp_sse);
      }
      BCL_ExampleCheck( helix_bodies.GetSize(), 2);

      // create a ShPtrVector of all_bodies which will be created from "all_sses"
      util::ShPtrVector< assemble::SSEGeometryInterface> all_bodies;

      // add the coord::Bodies of "all_sses" to "all_bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( all_sses.Begin()), sse_itr_end( all_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> temp_sse( ( *sse_itr)->Clone());
        all_bodies.PushBack( temp_sse);
      }
      BCL_ExampleCheck( all_bodies.GetSize(), 7);

      // define extent_specification
      linal::Vector3D extent_specification( 0.5, 0.5, 5.0);

      // create assemble::PickBodyExtent
      find::PickBodyExtent body_pick( extent_specification, extent_specification);

      // use "body_pick" to get a random body that meets the extent specification from the first body in helix_sses
      // this will always return the first helix in helix_sses (all the other ones differ too much in length)
      util::ShPtr< assemble::SSEGeometryInterface> extent_body_helix( body_pick.Pick( helix_bodies, **helix_sses.Begin()));

      BCL_ExampleIndirectAssert( extent_body_helix.IsDefined(), true, "Pick");

      BCL_MessageStd
      (
        "the sse body extents are: \n"
        + util::Format()( ( *helix_sses.Begin())->GetExtent( coord::GetAxes().e_Z))
        + " and the picked body extents are: \n" + util::Format()( extent_body_helix->GetExtent( coord::GetAxes().e_Z))
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          ( *helix_sses.Begin())->GetExtent( coord::GetAxes().e_Z),
          extent_body_helix->GetExtent( coord::GetAxes().e_Z)
        ),
        "the sse body extents are: \n"
        + util::Format()( ( *helix_sses.Begin())->GetExtent( coord::GetAxes().e_Z))
        + " and the picked body extents are: \n" + util::Format()( extent_body_helix->GetExtent( coord::GetAxes().e_Z))
      );

      // use "body_pick" to get a random body that meets the extent specification from the first body in helix_sses
      // this should return any of the first, second or third strand in all_sses (all the other ones differ too much in length)
      util::ShPtr< coord::GeometryInterface> extent_body( body_pick.Pick( all_bodies, **all_sses.Begin()));
      BCL_ExampleIndirectAssert( extent_body.IsDefined(), true, "Pick");

      BCL_MessageStd
      (
        "the sse body extents are: \n"
        + util::Format()( ( *all_sses.Begin())->GetExtent( coord::GetAxes().e_Z))
        + " and the picked body extents are: \n" + util::Format()( extent_body->GetExtent( coord::GetAxes().e_Z))
      );

      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance
        (
          ( *all_sses.Begin())->GetExtent( coord::GetAxes().e_Z),
          extent_body->GetExtent( coord::GetAxes().e_Z),
          extent_specification.Z()
        ),
        "the sse body extents are: \n"
        + util::Format()( ( *all_sses.Begin())->GetExtent( coord::GetAxes().e_Z))
        + " and the picked body extents are: \n" + util::Format()( extent_body->GetExtent( coord::GetAxes().e_Z))
      );

      // initialize counter that counts the number of bodies in all_bodies for which ExtentsWithinTolerance returns true
      size_t counter( 0);

      // calculate extents for the first sse in all_sses
      linal::Vector3D first_sse_extent( ( *all_sses.Begin())->GetExtents());

      //explicitly test the ExtentsWithinTolerance function
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator body_itr( all_bodies.Begin()),
          body_itr_end( all_bodies.End());
        body_itr != body_itr_end; ++body_itr
      )
      {
        if( body_pick.ExtentsWithinTolerance( first_sse_extent, ( *body_itr)->GetExtents()))
        {
          ++counter;
        }
      }

      BCL_ExampleIndirectCheck( counter, 4, "ExtentsWithinTolerance");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFindPickBodyExtent

  const ExampleClass::EnumType ExampleFindPickBodyExtent::s_Instance
  (
    GetExamples().AddEnum( ExampleFindPickBodyExtent())
  );

} // namespace bcl
