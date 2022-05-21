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
#include "score/bcl_score_restraint_body_protein_model.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "restraint/bcl_restraint_assignment.h"
#include "restraint/bcl_restraint_body.h"
#include "restraint/bcl_restraint_contains_body_origin.h"
#include "score/bcl_score_body_assignment.h"
#include "score/bcl_score_body_extent_position_agreement.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_body_protein_model.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintBodyProteinModel :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintBodyProteinModel *Clone() const
    { return new ExampleScoreRestraintBodyProteinModel( *this);}

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
      // print static classname
      BCL_MessageStd
      (
        "static classname is " + GetStaticClassName< score::RestraintBodyProteinModel>()
      );

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

      // create a ShPtrVector of bodies which will be created from "helix_sses" and will be used to creat the
      // assignment used for testing purposes
      util::ShPtr< util::ShPtrVector< assemble::SSEGeometry> > bodies( new util::ShPtrVector< assemble::SSEGeometry>());

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

      // create GroupCollection "group_collection" to hold elements of "helix_sses" assigned with "bodies"
      restraint::GroupCollection< size_t, assemble::SSEGeometry> group_collection;

      // insert the first helix of "helix_sses" assigned with the second body in "bodies"
      group_collection.Insert( 1, restraint::Group< assemble::SSEGeometry>( 1, assemble::SSEGeometry( **helix_sses.Begin())));

      // insert the second helix of "helix_sses" assigned with the first body in "bodies"
      group_collection.Insert( 0, restraint::Group< assemble::SSEGeometry>( 1, assemble::SSEGeometry( **++helix_sses.Begin())));

      // create Assignment "assignment" for scoring test purposes
      restraint::Assignment< util::ShPtrVector< assemble::SSEGeometry>, size_t, assemble::SSEGeometry> assignment
      (
        bodies, group_collection
      );

      // create restraint::Body out of "bodies" and way to check the occupancies
      util::ShPtr< restraint::Body> body_restraint
      (
        new restraint::Body
        (
          // initialize with "bodies"
          bodies,
          // initialize with ShPtr to a restraint::ContainsBodyOrigin to determine occupancies
          util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
          (
            new restraint::ContainsBodyOrigin()
          )
        )
      );

      // create RestraintBodyProteinModel "score"
      score::RestraintBodyProteinModel score
      (
        util::ShPtr< util::ShPtrVector< restraint::Body> >
        (
          new util::ShPtrVector< restraint::Body>( 1, body_restraint)
        ),
        score::BodyAssignment( score::BodyExtentPositionAgreement())
      );

      // create double "calculated_score" and initialize with score calculated by "score" for "protein_model"
//      double calculated_score( score( protein_model));
//
//      // create double "expected_score" and initialize with zero
//      double expected_score( 0);
//
//      // make sure that "calculated_score" is correct
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        math::EqualWithinTolerance( expected_score, calculated_score),
//        " expected_score is " + util::Format()( expected_score)
//        + " calculated score is " + util::Format()( calculated_score)
//      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRestraintBodyProteinModel

  const ExampleClass::EnumType ExampleScoreRestraintBodyProteinModel::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintBodyProteinModel())
  );

} // namespace bcl
