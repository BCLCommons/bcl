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
#include "score/bcl_score_restraint_distance_spin_label.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "score/bcl_score_restraint_atom_attraction.h"
#include "score/bcl_score_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_distance_spin_label.cpp
  //!
  //! @author alexanns
  //! @date Feb 09, 2009
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintDistanceSpinLabel :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintDistanceSpinLabel *Clone() const
    { return new ExampleScoreRestraintDistanceSpinLabel( *this);}

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
      // test Default constructor
      BCL_MessageStd( "test Default constructor");
      score::RestraintDistanceSpinLabel def_constr;

      // test ScoreDistance function
      BCL_MessageStd( "test ScoreDistance function");
      // SL-CB too positive
      BCL_MessageStd( "test ScoreDistance function SL-CB too positive");
      double score( def_constr.ScoreDistance( 30.0, 60.0));
      double correct_score( 0);
      BCL_MessageDbg
      (
        "calculated score is " + util::Format()( score) +
        " but correct score is " + util::Format()( correct_score)
      );
      BCL_ExampleCheckWithinAbsTolerance( correct_score, score, 1.0e-15);
      // SL-CB too negative
      BCL_MessageStd( "test ScoreDistance function SL-CB too negative");
      score = def_constr.ScoreDistance( 60.0, 30.0);
      BCL_MessageDbg
      (
        "calculated score is " + util::Format()( score) +
        " but correct score is " + util::Format()( correct_score)
      );
      BCL_ExampleCheckWithinTolerance( correct_score, score, 0.001);

      // TODO Fix the bug that follows this, which occurs because of an assertion in restraint::Distance
      BCL_MessageStd( "test ScoreDistance function SL-CB just right");
      BCL_ExampleCheckWithinAbsTolerance( def_constr.ScoreDistance( 18, 24.0), -0.935243, 0.000001);
      BCL_ExampleCheck( def_constr.ScoreDistance( 18, 24.0) < def_constr.ScoreDistance( 24, 18.0), true);

      // print out dsl-dcb and the score values for plotting if message level is debug
      const double dsl( 50.0);
      for( double dcb( 0.0); dcb <= 100; dcb += 0.1)
      {
        BCL_MessageDbg
        (
          "Dsl - Dcb example value : " + util::Format()( dsl - dcb) + " score : " +
          util::Format()( def_constr.ScoreDistance( dcb, dsl))
        );
      }

      // GetEPRDistanceScores
      const util::ShPtrVector< score::ProteinModel> scores
      (
        restraint::EPRDistanceData().GetScores()
      );
      BCL_ExampleAssert( scores.GetSize(), 3);
      // test knowledge based potential part of score
      {
        const util::ShPtr< score::RestraintAtomDistance> kb_potential_model_score( scores( 0));
        BCL_ExampleAssert( kb_potential_model_score.IsDefined(), true);
        const util::ShPtr< score::RestraintDistanceSpinLabel> kb_potential
        (
          kb_potential_model_score->GetScore()->Clone()
        );
        BCL_ExampleAssert( kb_potential.IsDefined(), true);
        BCL_ExampleCheckWithinAbsTolerance( kb_potential->ScoreDistance( 12, 24.5), 0.0, 0.0000001);
        BCL_ExampleCheckWithinAbsTolerance( kb_potential->ScoreDistance( 24.5, 12), 0.0, 0.0000001);
        BCL_ExampleCheck( kb_potential->ScoreDistance( 18.5, 12) > kb_potential->ScoreDistance( 14.5, 12), true);
        for( double x = -30.0; x < 30; x += 0.5)
        {
          const double current_score( kb_potential->ScoreDistance( 0, x));
          BCL_MessageDbg( "testKB: x " + util::Format()( x) + " sc " + util::Format()( current_score));
        }
      }
      const double width( score::RestraintAtomAttraction::GetDefaultTransitionWidth());
      // test penalty on left side of knowledge-based potential
      {
        const util::ShPtr< score::RestraintAtomDistance> left_penalty_model_score( scores( 1));
        BCL_ExampleAssert( left_penalty_model_score.IsDefined(), true);
        const util::ShPtr< score::RestraintAtomAttraction> left_penalty
        (
          left_penalty_model_score->GetScore()->Clone()
        );
        BCL_ExampleAssert( left_penalty.IsDefined(), true);
        BCL_ExampleCheck( left_penalty->GetFunction()( -8.0), -1);
        BCL_ExampleCheck( left_penalty->GetFunction()( -8.0 - 2 * width), 0);
        BCL_ExampleCheck
        (
          left_penalty->GetFunction()
          (
            -8.0 - width / 4.0) > left_penalty->GetFunction()( -8.0 - width / 6.0
          ),
          true
        );
        for( double x = -30.0; x < 30; x += 0.5)
        {
          const double current_score( left_penalty->GetFunction()( x));
          BCL_MessageDbg( "testLS: x " + util::Format()( x) + " sc " + util::Format()( current_score));
        }
      }
      // test penalty on right side of knowledge-based potential
      {
        BCL_MessageDbg( "right boundary score example test");
        const util::ShPtr< score::RestraintAtomDistance> right_penalty_model_score( scores( 2));
        BCL_ExampleAssert( right_penalty_model_score.IsDefined(), true);
        const util::ShPtr< score::RestraintAtomAttraction> right_penalty
        (
          right_penalty_model_score->GetScore()->Clone()
        );
        BCL_ExampleAssert( right_penalty.IsDefined(), true);
        BCL_ExampleCheck( right_penalty->GetFunction()( 8.0), -1);
        BCL_ExampleCheck( right_penalty->GetFunction()( 8.0 + 2 * width), 0);
        BCL_ExampleCheck
        (
          right_penalty->GetFunction()
          (
            8.0 + width / 4.0) > right_penalty->GetFunction()( 8.0 + width / 6.0
          ), true
        );
        for( double x = -30.0; x < 30; x += 0.5)
        {
          const double current_score( right_penalty->GetFunction()( x));
          BCL_MessageDbg( "testRS: x " + util::Format()( x) + " sc " + util::Format()( current_score));
        }
      }

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRestraintDistanceSpinLabel

  const ExampleClass::EnumType ExampleScoreRestraintDistanceSpinLabel::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintDistanceSpinLabel())
  );

} // namespace bcl
