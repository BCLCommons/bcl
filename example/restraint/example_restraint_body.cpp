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
#include "restraint/bcl_restraint_body.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_contains_body_origin.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_body.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintBody :
    public ExampleInterface
  {
  public:

    ExampleRestraintBody *Clone() const
    { return new ExampleRestraintBody( *this);}

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
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // get the helix secondary structure elements of "protein_model"
      util::SiPtrVector< const assemble::SSE> helix_sses
      (
        protein_model.GetChains()( 0)->GetSSEs( biol::GetSSTypes().HELIX)
      );
      BCL_MessageStd( "here 1 helix_sses size is " + util::Format()( helix_sses.GetSize()));
      // create a ShPtrVector of bodies which will be created from "helix_sses" and will be used to creat the body
      // restraint
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > bodies( new util::ShPtrVector< assemble::SSEGeometryInterface>());

      // add the coord::Bodies of "helix_sses" to "bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( helix_sses.Begin()), sse_itr_end( helix_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        util::ShPtr< assemble::SSE> temp_sse( ( *sse_itr)->Clone());
        bodies->PushBack( temp_sse);
      }
      BCL_MessageStd
      (
        "center of first restraint body "
        + util::Format()( bodies->operator()( 0)->GetCenter()) + "\n "
        + util::Format()( bodies->operator()( 0)->GetExtent( coord::GetAxes().e_Z))
        + " center of second restraint body "
        + util::Format()( bodies->operator()( 1)->GetCenter()) + "\n "
        + util::Format()( bodies->operator()( 1)->GetExtent( coord::GetAxes().e_Z))
      );

      // create restraint::Body out of "bodies" and way to check the occupancies
      restraint::Body body_restraint
      (
        // initialize with "bodies"
        bodies,
        // initialize with ShPtr to a restraint::ContainsBodyOrigin to determine occupancies
        util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
        (
          new restraint::ContainsBodyOrigin()
        )
      );

      // create SiPtrVector of assemble::SSEGeometryInterface which will hold the bodies of "protein_model"
      util::SiPtrVector< const assemble::SSEGeometryInterface> protein_model_bodies;

      // add the coord::Bodies of "helix_sses" to "protein_model_bodies"
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator sse_itr( helix_sses.Begin()), sse_itr_end( helix_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        protein_model_bodies.PushBack( **sse_itr);
      }

      // print the body information of the two SSEs from "protein_model
      BCL_MessageStd
      (
        "center of first protein body "
        + util::Format()( protein_model_bodies( 0)->GetCenter()) + "\n "
        + util::Format()( protein_model_bodies( 0)->GetExtent( coord::GetAxes().e_Z))
        + " center of second protein body "
        + util::Format()( protein_model_bodies( 1)->GetCenter()) + "\n "
        + util::Format()( protein_model_bodies( 1)->GetExtent( coord::GetAxes().e_Z))
      );

      // get unoccupied restraints: should be zero; all should be occupied
      BCL_Example_Check
      (
        !( body_restraint.GetUnoccupied( protein_model_bodies).GetSize()),
        " all restraints should be occupied but "
        + util::Format()( body_restraint.GetUnoccupied( protein_model_bodies).GetSize())
        + " are not occupied"
      );

      // change the origin of first body
//      bodies->operator()( 0)->SetOrigin( linal::Vector3D( 5.0, 6.0, 7.0));
      // now get unoccupied restraints: should be one
//      BCL_Assert
//      (
//        body_restraint.GetUnoccupied( protein_model_bodies).GetSize() ==  1,
//        " one restraints should be unoccupied but "
//        + util::Format()( body_restraint.GetUnoccupied( protein_model_bodies).GetSize())
//        + " are unoccupied"
//      );

//      // reset the origin of the first body to its original origin
//      bodies->operator()( 0)->SetOrigin( linal::Vector3D(  37.2137, 27.377, 13.1844));

      //! create Assignment "body_assignment" and initialize with the assignments of "body_restraint"
      //! with "protein_model_bodies"
//      restraint::SSEAssignment body_assignment
//      (
//        body_restraint.GenerateAssignment( protein_model_bodies)
//      );
//
//      // make sure that "body_assignment" has the correct number of members (should be 2)
//      BCL_Assert
//      (
//        body_assignment.GetGroupCollection().TotalDepth() == 2,
//        "should be two bodies from \"protein_model_bodies\" but instead there are "
//        + util::Format()( body_assignment.GetGroupCollection().TotalDepth())
//      );
//
//      // create size_t "index" and initialize with the key of the first element in "body_assignment" (should be 0)
//      size_t index( body_assignment.GetGroupCollection().Begin()->first);
//
//      // make sure that "body_assignment" is correct
//      BCL_Assert
//      (
//
////            body_assignment.GetRestraint()->operator()( index)->GetMatrix() ==
////         ( *body_assignment.GetGroupCollection().Begin()->second.Begin())->GetMatrix()
////
////      &&
//          math::EqualWithinTolerance
//          (
//            body_assignment.GetRestraint()->operator()( index)->GetExtents(),
//            ( *body_assignment.GetGroupCollection().Begin()->second.Begin())->GetExtents()
//          )
//        && math::EqualWithinTolerance
//          (
//            body_assignment.GetRestraint()->operator()( index)->GetOrigin(),
//            ( *body_assignment.GetGroupCollection().Begin()->second.Begin())->GetOrigin()
//          )
//        && index == 0,
//        " the key should be zero but is " + util::Format()( index)
//        + " and this key should give the body of \"body_restraint\" which is occupied by the body accessed by key zero"
//        + " \nthe origin of the restraint body " + util::Format()( index)
//        + " is \n" + util::Format()( body_assignment.GetRestraint()->operator()( index)->GetOrigin())
//        + " \nand the extents are\n "
//        + util::Format()( body_assignment.GetRestraint()->operator()( index)->GetExtents())
//        + "\n and the matrix is \n"
//        + util::Format()( body_assignment.GetRestraint()->operator()( index)->GetMatrix())
//        + " \nThe origin of the body accessed by key " + util::Format()( index)
//        + " \n is " + util::Format()( ( *body_assignment.GetGroupCollection().Begin()->second.Begin())->GetOrigin())
//        + " \n and the extents are \n "
//        + util::Format()( ( *body_assignment.GetGroupCollection().Begin()->second.Begin())->GetExtents())
//        + " \n and the matrix is \n"
//        + util::Format()( ( *body_assignment.GetGroupCollection().Begin()->second.Begin())->GetMatrix())
//      );
//
//      // change the origin of the first body in "bodies" so that it will be unoccupied
//      bodies->operator()( 0)->SetOrigin( linal::Vector3D( 5.0, 6.0, 7.0));
//
//      //! create Assignment "body_assignment_b" and initialize with the assignments of "body_restraint"
//      //! with "protein_model_bodies"
//      restraint::SSEAssignment body_assignment_b
//      (
//        body_restraint.GenerateAssignment( protein_model_bodies)
//      );
//
//      // make sure that "body_assignment_b" has the correct number of members (should be 1)
//      BCL_Assert
//      (
//        body_assignment_b.GetGroupCollection().TotalDepth() == 1,
//        "should be 1 bodies from \"protein_model_bodies\" but instead there are "
//        + util::Format()( body_assignment.GetGroupCollection().TotalDepth())
//      );
//
//      // create size_t "index_b" and initialize with the key of the first element in "body_assignment" (should be 1)
//      size_t index_b( body_assignment_b.GetGroupCollection().Begin()->first);
//      BCL_Assert
//      (
//
//            body_assignment.GetRestraint()->operator()( index_b)->GetMatrix() ==
//            ( *body_assignment_b.GetGroupCollection().Begin()->second.Begin())->GetMatrix()
//
//        && math::EqualWithinTolerance
//          (
//            body_assignment_b.GetRestraint()->operator()( index_b)->GetExtents(),
//            ( *body_assignment_b.GetGroupCollection().Begin()->second.Begin())->GetExtents()
//          )
//        && math::EqualWithinTolerance
//          (
//            body_assignment_b.GetRestraint()->operator()( index_b)->GetOrigin(),
//            ( *body_assignment_b.GetGroupCollection().Begin()->second.Begin())->GetOrigin()
//          )
//        && index_b == 1,
//        " the key should be 1 but is " + util::Format()( index_b)
//        + " and this key should give the body of \"body_restraint\" which is occupied by the body accessed by key 1"
//        + " \nthe origin of the restraint body " + util::Format()( index_b)
//        + " is \n" + util::Format()( body_assignment_b.GetRestraint()->operator()( index_b)->GetOrigin())
//        + " \nand the extents are\n "
//        + util::Format()( body_assignment_b.GetRestraint()->operator()( index_b)->GetExtents())
//        + "\n and the matrix is \n"
//        + util::Format()( body_assignment_b.GetRestraint()->operator()( index_b)->GetMatrix())
//        + " \nThe origin of the body accessed by key " + util::Format()( index_b)
//        + " \n is " + util::Format()( ( *body_assignment_b.GetGroupCollection().Begin()->second.Begin())->GetOrigin())
//        + " \n and the extents are \n "
//        + util::Format()( ( *body_assignment_b.GetGroupCollection().Begin()->second.Begin())->GetExtents())
//        + " \n and the matrix is \n"
//        + util::Format()( ( *body_assignment_b.GetGroupCollection().Begin()->second.Begin())->GetMatrix())
//      );

      // write mutated proteinmodel to an example pdb
      //BCL_MessageStd( "write mutated_swap.pdb");
      //write.open( AddExamplePathToFilename( "mutated_swap.pdb").c_str());
      //factory.WriteModelToPDB( mutated_model, write);
      //io::File::CloseClearFStream( write);

       // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintBody

  const ExampleClass::EnumType ExampleRestraintBody::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintBody())
  );

} // namespace bcl
