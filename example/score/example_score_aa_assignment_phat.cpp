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
#include "score/bcl_score_aa_assignment_phat.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_phat.cpp
  //!
  //! @author dongen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentPhat :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentPhat *Clone() const
    {
      return new ExampleScoreAAAssignmentPhat( *this);
    }

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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor from a phat_table
      score::AAAssignmentPHAT AAA_assign_phat1( score::AAAssignmentPHAT::e_PHAT_85);
      // check that "AAA_assign_phat1" was constructed properly
      BCL_Example_Check
      (
        AAA_assign_phat1.GetTableType() == score::AAAssignmentPHAT::e_PHAT_85,
        "GetTableType gives" +
        score::AAAssignmentPHAT::GetTableTypeString( AAA_assign_phat1.GetTableType()) +
        "but should give " +
        score::AAAssignmentPHAT::GetTableTypeString( score::AAAssignmentPHAT::e_PHAT_85)
      );

      // check copy constructor
      score::AAAssignmentPHAT AAA_assign_phat_copy_constr( AAA_assign_phat1);
      // check that "AAA_assign_phat1" was constructed properly
      BCL_Example_Check
      (
        AAA_assign_phat_copy_constr.GetTableType() == AAA_assign_phat1.GetTableType(),
        "GetTableType gives" +
        score::AAAssignmentPHAT::GetTableTypeString( AAA_assign_phat_copy_constr.GetTableType())
        + "but should give " +
        score::AAAssignmentPHAT::GetTableTypeString( AAA_assign_phat1.GetTableType())
      );

      // clone
      util::ShPtr< score::AAAssignmentPHAT> sp_AAA_assign_phat2
      (
        AAA_assign_phat1.Clone()
      );

      BCL_Example_Check
      (
        sp_AAA_assign_phat2->GetTableType() == AAA_assign_phat1.GetTableType(),
        "GetTableType gives" +
        score::AAAssignmentPHAT::GetTableTypeString( sp_AAA_assign_phat2->GetTableType()) +
        "but should give " +
        score::AAAssignmentPHAT::GetTableTypeString( AAA_assign_phat1.GetTableType())
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::score::AAAssignmentPHAT");
      BCL_Example_Check
      (
        GetStaticClassName< score::AAAssignmentPHAT>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< score::AAAssignmentPHAT>() + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< score::AAAssignmentPHAT>() == sp_AAA_assign_phat2->GetClassIdentifier(),
        "GetClassIdentifier gives " + sp_AAA_assign_phat2->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test Probability, but first we must set up the parameters.

      // create a scoring function with each of the PHAT tables
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_phat_85( new score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_85));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_phat_80( new score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_80));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_phat_75( new score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_75));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_phat_70( new score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_70));

      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // Output sample probabilites for each PHAT table
      BCL_MessageStd
      (
        "phat 85 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()
        (
          score::AAAssignmentPHAT::Probability
          (
            score::AAAssignmentPHAT::e_PHAT_85,
            seq1.GetFirstAA()->GetType(),
            seq2.GetFirstAA()->GetType()
          )
        )
      );
      BCL_MessageStd
      (
        "phat 80 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()
        (
          score::AAAssignmentPHAT::Probability
          (
            score::AAAssignmentPHAT::e_PHAT_80,
            seq1.GetFirstAA()->GetType(),
            seq2.GetFirstAA()->GetType()
          )
        )
      );
      BCL_MessageStd
      (
        "phat 75 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()
        (
          score::AAAssignmentPHAT::Probability
          (
            score::AAAssignmentPHAT::e_PHAT_75,
            seq1.GetFirstAA()->GetType(),
            seq2.GetFirstAA()->GetType()
          )
        )
      );
      BCL_MessageStd
      (
        "phat 70 matrix entry of first aa of 1fms_ and first aa of 1f5mA " +
        util::Format()
        (
          score::AAAssignmentPHAT::Probability
          (
            score::AAAssignmentPHAT::e_PHAT_70,
            seq1.GetFirstAA()->GetType(),
            seq2.GetFirstAA()->GetType()
          )
        )
      );

      //check the probability of PHAT_70
      const double probability_phat70
      (
        score::AAAssignmentPHAT::Probability
        (
          score::AAAssignmentPHAT::e_PHAT_70, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()
        )
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( probability_phat70, 0.0792447),
        "Probability gives " + util::Format()( probability_phat70) +
        " but should give " + util::Format()( 0.0792447)
      );

      //check the probability of PHAT_75
      const double probability_phat75
      (
        score::AAAssignmentPHAT::Probability
        (
          score::AAAssignmentPHAT::e_PHAT_75, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType()
        )
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( probability_phat75, 0.0629463),
        "Probability gives " + util::Format()( probability_phat75) +
        " but should give " + util::Format()( 0.0629463)
      );

      //check the probability of PHAT_80
      const double probability_phat80
      (
        score::AAAssignmentPHAT::Probability( score::AAAssignmentPHAT::e_PHAT_80, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType())
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( probability_phat80, 0.0792447),
        "Probability gives " + util::Format()( probability_phat80) +
        " but should give " + util::Format()( 0.0792447)
      );

      //check the probability of PHAT_85
      const double probability_phat85
      (
        score::AAAssignmentPHAT::Probability( score::AAAssignmentPHAT::e_PHAT_85, seq1.GetFirstAA()->GetType(), seq2.GetFirstAA()->GetType())
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( probability_phat85, 0.0629463),
        "Probability gives " + util::Format()( probability_phat85) +
        " but should give " + util::Format()( 0.0629463)
      );

      //test operator() for each table
      const double assign_score_85( score_assignment_phat_85->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_80( score_assignment_phat_80->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_75( score_assignment_phat_75->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));
      const double assign_score_70( score_assignment_phat_70->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));

      // create expected scores for each table
      const double expected_assign_score_85( 0.1);
      const double expected_assign_score_80( 0.2);
      const double expected_assign_score_75( 0.1);
      const double expected_assign_score_70( 0.2);

      // output the scores for each table
      BCL_MessageStd( "assignment score 85 is " + util::Format()( assign_score_85));
      BCL_MessageStd( "assignment score 80 is " + util::Format()( assign_score_80));
      BCL_MessageStd( "assignment score 75 is " + util::Format()( assign_score_75));
      BCL_MessageStd( "assignment score 70 is " + util::Format()( assign_score_70));

      // make sure that the scores are what is expected
      BCL_Example_Check
        (
          math::EqualWithinTolerance( expected_assign_score_85, assign_score_85),
          "Assignment score 85 does not return expected score of " + util::Format()( assign_score_85)
        );
      BCL_Example_Check
        (
          math::EqualWithinTolerance( expected_assign_score_80, assign_score_80),
          "Assignment score 80 does not return expected score of " + util::Format()( assign_score_80)
        );
      BCL_Example_Check
        (
          math::EqualWithinTolerance( expected_assign_score_75, assign_score_75),
          "Assignment score 75 does not return expected score of " + util::Format()( assign_score_75)
        );
      BCL_Example_Check
        (
          math::EqualWithinTolerance( expected_assign_score_70, assign_score_70),
          "Assignment score 70 does not return expected score of " + util::Format()( assign_score_70)
        );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( AAA_assign_phat1);
      // create instance of class "AAAssignmentPHAT" and read from file
      score::AAAssignmentPHAT read_score;
      ReadBCLObject( read_score);

      BCL_MessageStd( "read in AAAssignmentPHAT : " + util::Format()( read_score));

      // make sure the "read_score" is correct
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
      read_score_clone( read_score.Clone());

      // create double "read_calculated_score" to hold the score calculated by "read_score"
      const double read_calculated_score( read_score_clone->operator()( *seq1.GetFirstAA(), *seq2.GetFirstAA()));

      // make sure that the score calculated by "read_calculated_score" is what is expected
      BCL_ExampleIndirectCheck( read_calculated_score, expected_assign_score_85, "I/O");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentPhat

  const ExampleClass::EnumType ExampleScoreAAAssignmentPhat::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentPhat())
  );

} // namespace bcl
