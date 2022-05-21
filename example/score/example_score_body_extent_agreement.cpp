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
#include "score/bcl_score_body_extent_agreement.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_body_extent_agreement.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreBodyExtentAgreement :
    public ExampleInterface
  {
  public:

    ExampleScoreBodyExtentAgreement *Clone() const
    { return new ExampleScoreBodyExtentAgreement( *this);}

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
        "static classname is " + GetStaticClassName< score::BodyExtentAgreement>()
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
      const util::SiPtrVector< const assemble::SSE> helix_sses
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
        +"\n\nThe origin of the second helix in \"helix_sses\" is : \n "
        + util::Format()( ( *( ++helix_sses.Begin()))->GetCenter())
        + "\n and the extents of the second helix in \"helix_sses\" is : \n"
        + util::Format()( ( *( ++helix_sses.Begin()))->GetExtent( coord::GetAxes().e_Z))
      );

      // create ScoreBodyExtentPositionAgreement "score_agreement"
      score::BodyExtentAgreement score_agreement( 1.5, 1.0, 0.5, 1.5, -1.0, coord::GetAxes().e_Z);

      // create ShPtr to helix_sses.Begin() and ++helix_sses.Begin()
      const assemble::SSE &first_helix( **helix_sses.Begin());
      const assemble::SSE &second_helix( **( ++helix_sses.Begin()));

      // test GetExtentDeviation
      const double deviation( score_agreement.GetExtentDeviation( first_helix, second_helix)
      );
      const double expected_deviation( 12.0198);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_deviation, deviation),
        "expected deviation is " + util::Format()( expected_deviation) + " .\n deviation is "
        + util::Format()( deviation)
      );

      // create double "agreement" to test calculating agreement when outside of the upper bounds of the restraint
      const double agreement( score_agreement( first_helix, second_helix));

      // make sure that the correct agreement was calculated
      const double expected_agreement( 10.0198);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_agreement, agreement),
        "agreement should be " + util::Format()( expected_agreement)
        + " but instead is " + util::Format()( agreement)
      );

      // create double "agreement_upper_tolerance" to test calculating agreement
      // when outside of the upper bounds of the restraint (again)
      score::BodyExtentAgreement score_agreement_upper_tolerance( 1.5, 1.0, 2.5, 1.5, -1.0, coord::GetAxes().e_Z);
      const double agreement_upper_tolerance( score_agreement_upper_tolerance( first_helix, second_helix));

      // make sure that the correct agreement was calculated
      const double expected_agreement_upper_tolerance( 8.0198);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_agreement_upper_tolerance, agreement_upper_tolerance),
        "agreement should be " + util::Format()( expected_agreement_upper_tolerance)
        + " but instead is " + util::Format()( agreement_upper_tolerance)
      );

      // create ScoreBodyExtentPositionAgreement "score_agreement" to check scoring in total agreement region
      score::BodyExtentAgreement score_agreement_total( 1.5, 1.0, 14.0, 1.5, -1.0, coord::GetAxes().e_Z);
      const double agreement_total( score_agreement_total( first_helix, second_helix));

      // make sure that the correct agreement was calculated
      const double expected_agreement_total( -1.0);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_agreement_total, agreement_total),
        "agreement should be " + util::Format()( expected_agreement_total)
        + " but instead is " + util::Format()( agreement_total)
      );

      // test transition region for upper bound
      for( double upper_tolerance( 10.0198); upper_tolerance <= 12.0198; upper_tolerance += 0.2)
      {
        score::BodyExtentAgreement score_agreement_transition_upper( 1.5, 1.0, upper_tolerance, 2.0, -1.0, coord::GetAxes().e_Z);
        const double agreement_transition_upper( score_agreement_transition_upper( first_helix, second_helix));

        BCL_Example_Check
        (
          agreement_transition_upper < 0.0 && agreement_transition_upper >= -1.0,
          "agreement should be between 0.0 and -1.0, but instead is " +
          util::Format()( agreement_transition_upper)
        );
        BCL_MessageStd
        (
          "upper tolerance: " + util::Format()( upper_tolerance) + ", agreement: " + util::Format()( agreement_transition_upper)
        );
      }

      // test transition region for lower bound
      for( double lower_tolerance( 10.0198); lower_tolerance <= 12.0198; lower_tolerance += 0.2)
      {
        score::BodyExtentAgreement score_agreement_transition_lower( lower_tolerance, 2.0, 2.5, 2.0, -1.0, coord::GetAxes().e_Z);
        const double agreement_transition_lower( score_agreement_transition_lower( second_helix, first_helix));

        BCL_Example_Check
        (
          agreement_transition_lower < 0.0 && agreement_transition_lower >= -1.0,
          "agreement should be between 0.0 and -1.0, but instead is " +
          util::Format()( agreement_transition_lower)
        );
        BCL_MessageStd
        (
          "lower tolerance: " + util::Format()( lower_tolerance) + ", agreement: " + util::Format()( agreement_transition_lower)
        );
      }

      // create double "agreement_lower_tolerance" to test calculating agreement
      // when outside of the lower bounds of the restraint
      score::BodyExtentAgreement score_agreement_lower_tolerance( 2.5, 1.5, 3.5, 4.5, -1.0, coord::GetAxes().e_Z);
      const double agreement_lower_tolerance( score_agreement_lower_tolerance( second_helix, first_helix));

      // make sure that the correct agreement was calculated
      const double expected_agreement_lower_tolerance( 8.0198);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_agreement_lower_tolerance, agreement_lower_tolerance),
        "agreement should be " + util::Format()( expected_agreement_lower_tolerance)
        + " but instead is " + util::Format()( agreement_lower_tolerance)
      );

      // create ScoreBodyExtentPositionAgreement "score_agreement"
      score::BodyExtentAgreement score_agreement_lower( 14.0, 1.5, 0.5, 1.5, -1.0, coord::GetAxes().e_Z);
      // create double "agreement_lower" to test calculating agreement when inside the lower bounds of the restraint
      const double agreement_lower( score_agreement_lower( second_helix, first_helix));

      // make sure that the correct agreement was calculated
      const double expected_agreement_lower( -1.0);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_agreement_lower, agreement_lower),
        "agreement should be " + util::Format()( expected_agreement_lower)
        + " but instead is " + util::Format()( agreement_lower)
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreBodyExtentAgreement

  const ExampleClass::EnumType ExampleScoreBodyExtentAgreement::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreBodyExtentAgreement())
  );

} // namespace bcl
