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
#include "score/bcl_score_aa_assignment_ss_prediction.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_assignment_ss_prediction.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAAssignmentSSPrediction :
    public ExampleInterface
  {
  public:

    ExampleScoreAAAssignmentSSPrediction *Clone() const
    { return new ExampleScoreAAAssignmentSSPrediction( *this);}

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
      // declare sspred::MethodHandler
      sspred::MethodHandler method_handler;

      // declare scoring functions
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_PSIPRED( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_PSIPRED));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_JUFO( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_JUFO));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_SAM( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_SAM));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_TMHMM( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_TMHMM));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_TMMOD( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_TMMOD));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_B2TMPRED( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_B2TMPRED));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_PROFTMB( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_PROFTMB));
      util::ShPtr< function::BinaryInterface< const biol::AABase, const biol::AABase, double> >
        score_assignment_ss_prediction_CONPRED( new score::AAAssignmentSSPrediction( sspred::GetMethods().e_CONPRED));

      io::IFStream read;

      // read seq from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.fasta"));
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_MessageStd( "reading PSIPRED");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.psipred_ss2"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read jufo ss prediction
      BCL_MessageStd( "reading JUFO");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.jufo"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_JUFO);
      io::File::CloseClearFStream( read);

      // read sam ss prediction
      BCL_MessageStd( "reading SAM");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.rdb6"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_SAM);
      io::File::CloseClearFStream( read);

      // read tmhmm ss prediction
      BCL_MessageStd( "reading TMHMM");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_TMHMM.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_TMHMM);
      io::File::CloseClearFStream( read);

      // read tmmod ss prediction
      BCL_MessageStd( "reading TMMOD");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_TMMOD.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_TMMOD);
      io::File::CloseClearFStream( read);

      // read B2TMPRED ss prediction
      BCL_MessageStd( "reading B2TMPRED");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_B2TMPRED.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_B2TMPRED);
      io::File::CloseClearFStream( read);

//      // read proftmb ss prediction
//      BCL_MessageStd( "reading PROFTMB");
//      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_PROFTMB.txt"));
//      seq.ReadSSPrediction( read, sspred::GetMethods().e_PROFTMB);
//      io::File::CloseClearFStream( read);

      // read conpred ss prediction
      BCL_MessageStd( "reading CONPRED");
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_CONPRED.txt"));
      method_handler.ReadPredictionsForAASequence( read, seq, sspred::GetMethods().e_CONPRED);
      io::File::CloseClearFStream( read);

      BCL_MessageStd
      (
        "aatype of first and second aa of 1eco_ " +
        seq.GetData()( 0)->GetType()->GetThreeLetterCode() + " " + seq.GetData()( 1)->GetType()->GetThreeLetterCode()
      );

      const double assign_score_PSIPRED( score_assignment_ss_prediction_PSIPRED->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
      const double assign_score_JUFO( score_assignment_ss_prediction_JUFO->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
      const double assign_score_SAM( score_assignment_ss_prediction_SAM->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
      const double assign_score_TMHMM( score_assignment_ss_prediction_TMHMM->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
      const double assign_score_TMMOD( score_assignment_ss_prediction_TMMOD->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
      const double assign_score_B2TMPRED( score_assignment_ss_prediction_B2TMPRED->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
//      const double assign_score_PROFTMB( score_assignment_ss_prediction_PROFTMB->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));
      const double assign_score_CONPRED( score_assignment_ss_prediction_CONPRED->operator()( *seq.GetData()( 0), *seq.GetData()( 1)));

      const double expected_assign_score_PSIPRED( 0.443059);
      const double expected_assign_score_JUFO( 0.41961);
      const double expected_assign_score_SAM( 0.351174);
      const double expected_assign_score_TMHMM( 0.477121);
      const double expected_assign_score_TMMOD( 0.477121);
      const double expected_assign_score_B2TMPRED( 0.477121);
//      const double expected_assign_score_PROFTMB(0);
      const double expected_assign_score_CONPRED( 0.477121);

      BCL_MessageStd( "assignment score PSIPRED is " + util::Format()( assign_score_PSIPRED));
      BCL_MessageStd( "assignment score JUFO is " + util::Format()( assign_score_JUFO));
      BCL_MessageStd( "assignment score SAM is " + util::Format()( assign_score_SAM));
      BCL_MessageStd( "assignment score TMHMM is " + util::Format()( assign_score_TMHMM));
      BCL_MessageStd( "assignment score TMMOD is " + util::Format()( assign_score_TMMOD));
      BCL_MessageStd( "assignment score B2TMPRED is " + util::Format()( assign_score_B2TMPRED));
//      BCL_MessageStd( "assignment score PROFTMB is " + util::Format()( assign_score_PROFTMB));
      BCL_MessageStd( "assignment score CONPRED is " + util::Format()( assign_score_CONPRED));

      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_PSIPRED, assign_score_PSIPRED),
        "Assignment score PSIPRED does not return expected score of " + util::Format()( expected_assign_score_PSIPRED)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_JUFO, assign_score_JUFO),
        "Assignment score JUFO does not return expected score of " + util::Format()( expected_assign_score_JUFO)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_SAM, assign_score_SAM),
        "Assignment score SAM does not return expected score of " + util::Format()( expected_assign_score_SAM)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_TMHMM, assign_score_TMHMM),
        "Assignment score TMHMM does not return expected score of " + util::Format()( expected_assign_score_TMHMM)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_TMMOD, assign_score_TMMOD),
        "Assignment score TMMOD does not return expected score of " + util::Format()( expected_assign_score_TMMOD)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_B2TMPRED, assign_score_B2TMPRED),
        "Assignment score B2TMPRED does not return expected score of " + util::Format()( expected_assign_score_B2TMPRED)
      );
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        math::EqualWithinTolerance( expected_assign_score_PROFTMB, assign_score_PROFTMB),
//        "Assignment score PROFTMB does not return expected score of " + util::Format()( expected_assign_score_PROFTMB)
//      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_assign_score_CONPRED, assign_score_CONPRED),
        "Assignment score CONPRED does not return expected score of " + util::Format()( expected_assign_score_CONPRED)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAAssignmentSSPrediction

  const ExampleClass::EnumType ExampleScoreAAAssignmentSSPrediction::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAAssignmentSSPrediction())
  );

} // namespace bcl
