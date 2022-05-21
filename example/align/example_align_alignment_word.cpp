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
#include "align/bcl_align_alignment_word.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_handler_pir.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_blosum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_alignment_word.cpp
  //!
  //! @author heinzes1
  //! @date Mar 26, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignmentWord :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignmentWord *Clone() const
    {
      return new ExampleAlignAlignmentWord( *this);
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
      // prepare everything: read in a sequence to test
      io::IFStream read;

      // read seq from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence>
        seq( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // create alignment
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment( new align::AlignmentLeaf< biol::AABase>( seq));

      // create score: all scores for single sequence alignments have to be 0.0
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "2: AlignmentWord()")
      align::AlignmentWord< biol::AABase> default_word_alignment;
      bool exp_empty( true);
      size_t exp_size( 0);
      double exp_score( 0.0), exp_self_score( 0.0);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( default_word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "2: " + ConditionString( default_word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test constructor taking alignment
      BCL_MessageStd( "3: AlignmentWord(alignment)")
      align::AlignmentWord< biol::AABase> empty_word_alignment( alignment);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( empty_word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "3: " + ConditionString( empty_word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test constructor taking alignment and itrs
      BCL_MessageStd( "4: AlignmentWord(alignment, itr_begin, itr_end)")
      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
        itr_begin( alignment->GetAssignments().Begin()),
        itr_end( itr_begin);
      storage::AdvanceIterator( itr_begin, alignment->GetAssignments().End(), 4);
      storage::AdvanceIterator( itr_end, alignment->GetAssignments().End(), 8);
      align::AlignmentWord< biol::AABase> word_alignment( alignment, itr_begin, itr_end);
      exp_empty = false;
      exp_size = 4;
      exp_self_score = 2.2;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "4: " + ConditionString( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test ScoreWith
      BCL_MessageStd( "5: ScoreWith()")
      biol::AASequence score_sequence( biol::AASequenceFactory::BuildSequenceFromFASTAString( "AAAA"));
      double exp_score_with( -0.2);
      BCL_ExampleIndirectCheck
      (
        word_alignment.ScoreWith( score_sequence, assign_score),
        exp_score_with,
        "5: ScoreWith==" + util::Format()( word_alignment.ScoreWith( score_sequence, assign_score))
          + "!=" + util::Format()( exp_score_with)
      );

      // test Prepend
      BCL_MessageStd( "6: Prepend(assignment)")
      word_alignment.Prepend( *--itr_begin);
      exp_size = 5;
      exp_self_score = 2.8;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "6: " + ConditionString( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test Prepend list
      BCL_MessageStd( "7: Prepend(assignment_list)")
      util::ShPtrList< align::Assignment< biol::AABase> > assignment_list;
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator itr( alignment->GetAssignments().Begin());
        itr != itr_begin;
        ++itr
      )
      {
        assignment_list.PushBack( *itr);
      }
      word_alignment.Prepend( assignment_list);
      exp_size = 8;
      exp_self_score = 4.9;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "7: " + ConditionString( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test Append
      BCL_MessageStd( "8: Append(assignment)")
      word_alignment.Append( *itr_end++);
      exp_size = 9;
      exp_self_score = 5.6;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "8: " + ConditionString( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test Append list
      BCL_MessageStd( "9: Append(assignment_list)")
      assignment_list.Reset();
      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator new_end_itr( itr_end);
      storage::AdvanceIterator( new_end_itr, alignment->GetAssignments().End(), 4);
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator itr( itr_end);
        itr != new_end_itr;
        ++itr
      )
      {
        assignment_list.PushBack( *itr);
      }
      word_alignment.Append( assignment_list);
      exp_size = 13;
      exp_self_score = 8.0;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "9: " + ConditionString( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test ResetAssignments
      BCL_MessageStd( "10: ResetAssignments()")
      word_alignment.ResetAssignments();
      exp_empty = true;
      exp_size = 0;
      exp_self_score = 0.0;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score),
        true,
        "10: " + ConditionString( word_alignment, assign_score, exp_empty, exp_size, exp_score, exp_self_score)
      );

      // test IsSubAlignment with itself
      BCL_MessageStd( "11: IsSubAlignment( self)")
      BCL_ExampleIndirectCheck
      (
        alignment->IsSubAlignment( *alignment),
        true,
        "11: alignment is not a subalignment of itself, but it should"
      );

      // test IsSubAlignment with empty subalignment
      BCL_MessageStd( "12: IsSubAlignment( empty_alignment)")
      word_alignment = align::AlignmentWord< biol::AABase>( alignment); // set to defined empty state
      BCL_ExampleIndirectCheck
      (
        alignment->IsSubAlignment( word_alignment),
        false,
        "12: empty word_alignment is a subalignment of itself, but it should not"
      );

      // test IsSubAlignment with non-empty subalignment
      BCL_MessageStd( "13: IsSubAlignment( sub_alignment)")
      word_alignment = align::AlignmentWord< biol::AABase>( alignment, itr_begin, itr_end); // set to defined non-empty state
      BCL_ExampleIndirectCheck
      (
        alignment->IsSubAlignment( word_alignment),
        true,
        "13: non-empty word_alignment is not a subalignment of itself, but it should"
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    //! @brief prints some information about the given alignment word
    //! @param ALIGNMENT the alignment word to print
    //! @param ASSIGN_SCORE the assignment score to be used for scoring
    void PrintAlignmentData
    (
      const align::AlignmentWord< biol::AABase> &ALIGNMENT,
      const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE
    ) const
    {
      align::HandlerPIR< biol::AABase> handler;
      handler.WriteAlignment( util::GetLogger(), ALIGNMENT);
      BCL_MessageStd
      (
        "Size=" + util::Format()( ALIGNMENT.GetSize())
          + " | Empty=" + util::Format()( ALIGNMENT.IsEmpty())
          + " | Score=" + util::Format()( ALIGNMENT.Score( ASSIGN_SCORE))
          + " | SelfScore=" + util::Format()( ALIGNMENT.ScoreSelf( ASSIGN_SCORE))
      );
    }

    //! @brief checks if the alignment word and assignment score agree with the given expected values
    //! @param ALIGNMENT the alignment word to be used
    //! @param ASSIGN_SCORE the assignment score
    //! @param EXP_EMPTY the expected value if the alignment word should be empty
    //! @param EXP_SIZE the expected size of the alignment word
    //! @param EXP_SCORE the expected score of the alignment word using the assignment score
    //! @param EXP_SELF_SCORE the expected self score of the alignment word
    //! @return if all values agree with their expected values
    bool ConditionCheck
    (
      const align::AlignmentWord< biol::AABase> &ALIGNMENT,
      const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE,
      const bool EXP_EMPTY,
      const size_t EXP_SIZE,
      const double EXP_SCORE,
      const double EXP_SELF_SCORE
    ) const
    {
      return ALIGNMENT.IsEmpty() == EXP_EMPTY
        && ALIGNMENT.GetSize() == EXP_SIZE
        && math::EqualWithinTolerance( EXP_SCORE, ALIGNMENT.Score( ASSIGN_SCORE))
        && math::EqualWithinTolerance( EXP_SELF_SCORE, ALIGNMENT.ScoreSelf( ASSIGN_SCORE));
    }

    //! @brief generates a string for disagreeing values of the alignment word with the given expected values
    //! @param ALIGNMENT the alignment word to be used
    //! @param ASSIGN_SCORE the assignment score
    //! @param EXP_EMPTY the expected value if the alignment word should be empty
    //! @param EXP_SIZE the expected size of the alignment word
    //! @param EXP_SCORE the expected score of the alignment word using the assignment score
    //! @param EXP_SELF_SCORE the expected self score of the alignment word
    //! @return a string describing the disagreement
    const std::string ConditionString
    (
      const align::AlignmentWord< biol::AABase> &ALIGNMENT,
      const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE,
      const bool EXP_EMPTY,
      const size_t EXP_SIZE,
      const double EXP_SCORE,
      const double EXP_SELF_SCORE
    ) const
    {
      std::string result;
      if( ALIGNMENT.IsEmpty() != EXP_EMPTY)
      {
        result.append( "Empty==" + util::Format()( ALIGNMENT.IsEmpty()) + "!=" + util::Format()( EXP_EMPTY) + "; ");
      }
      if( ALIGNMENT.GetSize() != EXP_SIZE)
      {
        result.append( "Size==" + util::Format()( ALIGNMENT.GetSize()) + "!=" + util::Format()( EXP_SIZE) + "; ");
      }
      if( !math::EqualWithinTolerance( EXP_SCORE, ALIGNMENT.Score( ASSIGN_SCORE)))
      {
        result.append
        (
          "Score==" + util::Format()( ALIGNMENT.Score( ASSIGN_SCORE))
          + "!=" + util::Format()( EXP_SCORE) + "; "
        );
      }
      if( !math::EqualWithinTolerance( EXP_SELF_SCORE, ALIGNMENT.ScoreSelf( ASSIGN_SCORE)))
      {
        result.append
        (
          "SelfScore==" + util::Format()( ALIGNMENT.ScoreSelf( ASSIGN_SCORE))
          + "!=" + util::Format()( EXP_SELF_SCORE) + "; "
        );
      }
      return result;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignAlignmentWord

  const ExampleClass::EnumType ExampleAlignAlignmentWord::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignAlignmentWord())
  );

} // namespace bcl
