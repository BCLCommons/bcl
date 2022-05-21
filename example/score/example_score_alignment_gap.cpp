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
#include "score/bcl_score_alignment_gap.h"

// includes from bcl - sorted alphabetically
//#include "align/bcl_align_factory.h"
#include "align/bcl_align_handler_pir.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_alignment_gap.cpp
  //!
  //! @author heinzes1
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAlignmentGap :
    public ExampleInterface
  {
  public:

    ExampleScoreAlignmentGap *Clone() const
    {
      return new ExampleScoreAlignmentGap( *this);
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
      // prepare everything: read in sequences to test
      io::IFStream read;

      // read seq1 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      biol::AASequence seq1( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read seq2 from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      biol::AASequence seq2( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // write sequences
      seq1.WriteFasta( util::GetLogger());
      seq2.WriteFasta( util::GetLogger());

      // create align::Factory which will handle the amino acid sequences and create simple test alignment
//      align::Factory< biol::AABase> align_factory;
//      align::AlignmentSimple< biol::AABase> alignment1( align_factory.CopySequenceToAlignment( seq1));
//      align::AlignmentSimple< biol::AABase> alignment2( align_factory.CopySequenceToAlignment( seq2));
//      align::AlignmentSimple< biol::AABase> alignment( align::Factory< biol::AABase>::MergeAlignments( alignment1, alignment2));
//
//      // create an AlignmentHandlerPIR and output the Alignment in pir format
//      align::HandlerPIR< biol::AABase> alignment_handler_pir( 50);
//      alignment_handler_pir.WriteAlignment( util::GetLogger(), alignment);
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // create score::Assignment and score::AlignmentAssignment
//      util::ShPtr< math::FunctionInterfaceSerializable< align::Assignment< biol::AABase>, double> > sp_assign_score
//      (
//        new score::Assignment< biol::AABase>( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160))
//      );
//      util::ShPtr< score::AlignmentAssignment< biol::AABase> >
//        sp_align_score( new score::AlignmentAssignment< biol::AABase>( sp_assign_score));
//
//      // create score::AlignmentGap
////      util::ShPtr< score::AlignmentGap< biol::AABase> >
////        sp_align_gap_score( new score::AlignmentGap< biol::AABase>( -0.1));
//
//      // create sum of all alignment scores
//      util::ShPtr< math::SumFunctionMixin< align::AlignmentInterface< biol::AABase>, double> >
//        sp_align_sum_score( new math::SumFunctionMixin< align::AlignmentInterface< biol::AABase>, double>( sp_align_score, 1.0, 0.0));
////      sp_align_sum_score->NewOperand( sp_align_gap_score);
//
//      // create old assignment score for comparison
//      util::ShPtr< score::AssignmentWithGap< biol::AABase> > sp_assign_gap_score_old
//      (
//        new score::AssignmentWithGap< biol::AABase>
//        (
//          double( 0.0) * score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160),
//          -0.1, -0.1, -0.1, -0.1
//        )
//      );
//
//      // create old assignment score for comparison
//      util::ShPtr< score::AssignmentWithGap< biol::AABase> > sp_assign_sum_score_old
//      (
//        new score::AssignmentWithGap< biol::AABase>
//        (
//          score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160),
//          -0.1, -0.1, -0.1, -0.1
//        )
//      );
//
//    /////////////////
//    // data access //
//    /////////////////
//
//    ///////////////
//    // operators //
//    ///////////////
//
//      // get new gap score
////      double score( sp_align_gap_score->operator()( alignment));
//      double score;
//      // calculate old gap score by using the operator on each assignment
//      double score_old( 0.0);
////      for
////      (
////        align::AlignmentSimple< biol::AABase>::const_iterator itr( alignment.GetAssignments().Begin()), itr_end( alignment.GetAssignments().End());
////        itr != itr_end;
////        ++itr
////      )
////      {
////        score_old += sp_assign_gap_score_old->operator()( **itr);
////      }
////
////      BCL_MessageStd( "Score=" + util::Format()( score) + "|" + util::Format()( score_old));
////      BCL_Example_Check
////      (
////        ExampleClass::ExampleResult::e_Trivial,
////        math::EqualWithinTolerance( score, score_old),
////        "Alignment gap score differs: " + util::Format()( score) + "!=" + util::Format()( score_old) + "\n"
////      );
//
//      // get new score
//      score = sp_align_sum_score->operator()( alignment);
//      // calculate old score by using the operator on each assignment
//      score_old = 0.0;
//      for
//      (
//        align::AlignmentSimple< biol::AABase>::const_iterator itr( alignment.GetAssignments().Begin()), itr_end( alignment.GetAssignments().End());
//        itr != itr_end;
//        ++itr
//      )
//      {
//        // calculate score using the new and old method, and output assignment and its score
//        score_old += sp_assign_sum_score_old->operator()( **itr);
//      }
//
//      BCL_MessageStd( "Score=" + util::Format()( score) + "|" + util::Format()( score_old));
////      BCL_Example_Check
////      (
////        ExampleClass::ExampleResult::e_Trivial,
////        math::EqualWithinTolerance( score, score_old),
////        "Alignment sum score differs: " + util::Format()( score) + "!=" + util::Format()( score_old) + "\n"
////      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAlignmentGap

  const ExampleClass::EnumType ExampleScoreAlignmentGap::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAlignmentGap())
  );

} // namespace bcl
