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
#include "score/bcl_score_aa_assignment_pam.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_pam.cpp
  //!
  //! @author heinzes1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentPam :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentPam *Clone() const
    { return new ExampleScoreAAAssignmentPam( *this);}

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
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_pam_100( new score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_100));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_pam_120( new score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_120));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_pam_160( new score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_pam_250( new score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_250));

      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      BCL_MessageStd
      (
        "pam 100 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentPAM::Probability( score::AAAssignmentPAM::e_PAM_100, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );
      BCL_MessageStd
      (
        "pam 120 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentPAM::Probability( score::AAAssignmentPAM::e_PAM_120, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );
      BCL_MessageStd
      (
        "pam 160 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentPAM::Probability( score::AAAssignmentPAM::e_PAM_160, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );
      BCL_MessageStd
      (
        "pam 250 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()( score::AAAssignmentPAM::Probability( score::AAAssignmentPAM::e_PAM_250, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()))
      );

      const double assign_score_100( score_assignment_pam_100->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_120( score_assignment_pam_120->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_160( score_assignment_pam_160->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_250( score_assignment_pam_250->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double expected_assign_score_100( 0.0);
      const double expected_assign_score_120( 0.1);
      const double expected_assign_score_160( 0.1);
      const double expected_assign_score_250( 0.1);
      BCL_MessageStd( "assignment score 100 is " + util::Format()( assign_score_100));
      BCL_MessageStd( "assignment score 120 is " + util::Format()( assign_score_120));
      BCL_MessageStd( "assignment score 160 is " + util::Format()( assign_score_160));
      BCL_MessageStd( "assignment score 250 is " + util::Format()( assign_score_250));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_100, assign_score_100),
        "Assignment score 100 does not return expected score of " + util::Format()( assign_score_100)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_120, assign_score_120),
        "Assignment score 120 does not return expected score of " + util::Format()( assign_score_120)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_160, assign_score_160),
        "Assignment score 160 does not return expected score of " + util::Format()( assign_score_160)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_250, assign_score_250),
        "Assignment score 250 does not return expected score of " + util::Format()( assign_score_250)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( *score_assignment_pam_100);
      // create instance of class "AAAssignmentPAM" and read from file
      score::AAAssignmentPAM read_score;
      ReadBCLObject( read_score);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentPam

  const ExampleClass::EnumType ExampleScoreAAAssignmentPam::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentPam())
  );

} // namespace bcl
