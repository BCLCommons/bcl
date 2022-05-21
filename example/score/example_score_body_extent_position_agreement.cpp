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
#include "score/bcl_score_body_extent_position_agreement.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_body_extent_position_agreement.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreBodyExtentPositionAgreement :
    public ExampleInterface
  {
  public:

    ExampleScoreBodyExtentPositionAgreement *Clone() const
    { return new ExampleScoreBodyExtentPositionAgreement( *this);}

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
        "static classname is " + GetStaticClassName< score::BodyExtentPositionAgreement>()
      );

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

      // print the origin and extents of the two helices in "helix_sses"
      BCL_MessageStd
      (
        "The origin of the first helix in \"helix_sses\" is : \n "
        + util::Format()( ( *helix_sses.Begin())->GetCenter())
        + "\n and the extents of the first helix in \"helix_sses\" is : \n"
        + util::Format()( ( *helix_sses.Begin())->GetExtent( coord::GetAxes().e_Z))
        +"\n\nThe origin of the secpmd helix in \"helix_sses\" is : \n "
        + util::Format()( ( *( ++helix_sses.Begin()))->GetCenter())
        + "\n and the extents of the second helix in \"helix_sses\" is : \n"
        + util::Format()( ( *( ++helix_sses.Begin()))->GetExtent( coord::GetAxes().e_Z))
      );

      // create ScoreBodyExtentPositionAgreement "score_agreement"
      score::BodyExtentPositionAgreement score_agreement;

      // test GetExtentPosition
//      double position( score_agreement.GetExtentPosition( **helix_sses.Begin(), coord::GetAxes().e_X));
//      double expected_position( 39.3337);
//      BCL_Example_Check
//      (
//      ExampleClass::ExampleResult::e_Trivial,
//        math::EqualWithinTolerance( expected_position, position),
//        "expected position is " + util::Format()( expected_position) + " .\n Position is "
//        +util::Format()( position)
//      );

      // test GetAbsoluteExtentDeviation
//      double deviation
//      (
//        score_agreement.GetAbsoluteExtentDeviation( **helix_sses.Begin(), **( ++helix_sses.Begin()),  coord::GetAxes().e_Y)
//      );
//      double expected_deviation( 7.944);
//      BCL_Example_Check
//      (
//      ExampleClass::ExampleResult::e_Trivial,
//        math::EqualWithinTolerance( expected_deviation, deviation),
//        "expected deviation is " + util::Format()( expected_deviation) + " .\n deviation is "
//        +util::Format()( deviation)
//      );

      // create double "agreement" which is the agreement of the first and second helices of "helix_sses"
//      double agreement
//      (
//        score_agreement( storage::VectorND< 2, assemble::SSEGeometry>( **helix_sses.Begin(), **( ++helix_sses.Begin())))
//      );

      // make sure that the correct agreement was calculated
//      double expected_agreement( 22.69247);
//      BCL_Example_Check
//      (
//      ExampleClass::ExampleResult::e_Trivial,
//        math::EqualWithinTolerance( expected_agreement, agreement),
//        "agreement should be " +util::Format()( expected_agreement)
//        + " but instead is " + util::Format()( agreement)
//      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreBodyExtentPositionAgreement

  const ExampleClass::EnumType ExampleScoreBodyExtentPositionAgreement::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreBodyExtentPositionAgreement())
  );

} // namespace bcl
