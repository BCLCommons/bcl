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
#include "align/bcl_align_alignment_hit.h"

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
  //! @example example_align_alignment_hit.cpp
  //!
  //! @author heinzes1
  //! @date Apr 12, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignmentHit :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignmentHit *Clone() const
    {
      return new ExampleAlignAlignmentHit( *this);
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

      // create word alignment
      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
        itr_begin( alignment->GetAssignments().Begin()),
        itr_end( itr_begin);
      storage::AdvanceIterator( itr_begin, alignment->GetAssignments().End(), 4);
      storage::AdvanceIterator( itr_end, alignment->GetAssignments().End(), 8);
      util::ShPtr< align::AlignmentWord< biol::AABase> >
        alignment_word( new align::AlignmentWord< biol::AABase>( alignment, itr_begin, itr_end));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "2: AlignmentHit()")
      align::AlignmentHit< biol::AABase> default_alignment_hit;
      bool exp_empty( true);
      size_t exp_size( 0);
      double exp_score( 0.0);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( default_alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "2: " + ConditionString( default_alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // test constructor taking two word alignments
      BCL_MessageStd( "3: AlignmentHit( alignment, alignment)")
      align::AlignmentHit< biol::AABase> alignment_hit( alignment_word, alignment_word);
      exp_empty = false;
      exp_size = alignment_word->GetSize();
      exp_score = 2.2;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "3: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // Clone() and Empty() are tested below

    /////////////////
    // data access //
    /////////////////

      // test GetSequences
      BCL_MessageStd( "4: GetSequences()")
      const util::ShPtrList< align::SequenceInterface< biol::AABase> > sequences_list( alignment_hit.GetSequences());
      const size_t exp_list_size( 2);
      BCL_ExampleIndirectCheck
      (
        sequences_list.GetSize(),
        exp_list_size,
        "4: sequence_list.Size=" + util::Format()( sequences_list.GetSize()) + "!=" + util::Format()( exp_list_size)
      );

      // test GetChildAlignments
      BCL_MessageStd( "5: GetChildAlignments()")
      const util::ShPtrList< align::AlignmentInterface< biol::AABase> >
        child_alignment_list( alignment_hit.GetChildAlignments());
      BCL_ExampleIndirectCheck
      (
        child_alignment_list.GetSize(),
        exp_list_size,
        "5: child_alignment_list.Size=" + util::Format()( sequences_list.GetSize())
          + "!=" + util::Format()( exp_list_size)
      );

      // test GetDepth
      BCL_MessageStd( "6: GetDepth()")
      const size_t depth( alignment_hit.GetDepth());
      BCL_ExampleIndirectCheck
      (
        depth,
        alignment_word->GetDepth() * 2,
        "6: Depth=" + util::Format()( depth) + "!=" + util::Format()( alignment_word->GetDepth() * 2)
      );

      // GetSize and indirectly GetAssignments are tested in ConditionCheck() below

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // IsEmpty is tested in ConditionCheck below

      // test Prepend assignment
      BCL_MessageStd( "7: Prepend(assignment)")
      alignment_hit.Prepend( alignment_hit.GetPrependAssignment()); // prepend
      --itr_begin; // keep itr updated
      exp_size++;
      exp_score = 2.8;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "7: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // test Prepend assignment list
      BCL_MessageStd( "8: Prepend(assignment_list)")
      util::ShPtrList< align::Assignment< biol::AABase> > assignment_list;
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          itr( alignment->GetAssignments().Begin());
        itr != itr_begin;
        ++itr
      )
      {
        util::ShPtr< align::Assignment< biol::AABase> > assignment( new align::Assignment< biol::AABase>());
        assignment->Append( ( **itr).GetMembers()); // add members for each child alignment
        assignment->Append( ( **itr).GetMembers());
        assignment_list.Append( assignment);
      }
      alignment_hit.Prepend( assignment_list);
      exp_size += 3;
      exp_score = 4.9;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "8: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // test PrependNextAssignment
      BCL_MessageStd( "9: PrependNextAssignment()");
      alignment_hit.PrependNextAssignment(); // must not change anything, hit starts at beginning of complete alignment
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "9: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // test Append assignment
      BCL_MessageStd( "10: Append(assignment)")
      alignment_hit.Append( alignment_hit.GetAppendAssignment()); // append
      exp_size++;
      exp_score = 5.6;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "10: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // test Append assignment list
      BCL_MessageStd( "11: Append(assignment_list)")
      assignment_list.Reset(); // remove all previous assignments
      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator itr_new_end( itr_end);
      // add 2 assignments, itr_new_end points behind the last assignment
      storage::AdvanceIterator( itr_new_end, alignment->GetAssignments().End(), 3);
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator itr( itr_end);
        itr != itr_new_end;
        ++itr
      )
      {
        util::ShPtr< align::Assignment< biol::AABase> > assignment( new align::Assignment< biol::AABase>());
        assignment->Append( ( **itr).GetMembers()); // add members for each child alignment
        assignment->Append( ( **itr).GetMembers());
        assignment_list.Append( assignment);
      }
      alignment_hit.Append( assignment_list);
      exp_size += 2;
      exp_score = 7.0;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "11: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // test AppendNextAssignment
      BCL_MessageStd( "12: AppendNextAssignment()");
      alignment_hit.AppendNextAssignment();
      exp_size++;
      exp_score = 7.5;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "12: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
      );

      // Score is tested in ConditionCheck() below

      // test Clone()
      BCL_MessageStd( "13: Clone()")
      util::ShPtr< align::AlignmentHit< biol::AABase> > alignment_hit_clone( alignment_hit.Clone());
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( *alignment_hit_clone, assign_score, exp_empty, exp_size, exp_score),
        true,
        "13: " + ConditionString( *alignment_hit_clone, assign_score, exp_empty, exp_size, exp_score)
      );

      // test Empty()
      BCL_MessageStd( "14: Empty()")
      util::ShPtr< align::AlignmentHit< biol::AABase> > alignment_hit_empty( alignment_hit.Empty());
      exp_size = alignment_word->GetSize();
      exp_score = 2.2;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( *alignment_hit_empty, assign_score, exp_empty, exp_size, exp_score),
        true,
        "14: " + ConditionString( *alignment_hit_empty, assign_score, exp_empty, exp_size, exp_score)
      );

      // test ResetAssignments
      BCL_MessageStd( "15: ResetAssignments()")
      alignment_hit.ResetAssignments();
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_hit, assign_score, exp_empty, exp_size, exp_score),
        true,
        "15: " + ConditionString( alignment_hit, assign_score, exp_empty, exp_size, exp_score)
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
    void PrintAlignmentData
    (
      const align::AlignmentHit< biol::AABase> &ALIGNMENT,
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
      );
    }

    //! @brief tests if the alignment hit properties agree with their expected values
    //! @param ALIGNMENT the alignment hit to be used
    //! @param ASSIGN_SCORE the assignment score
    //! @param EXP_EMPTY the expected value if the alignment hit should be empty
    //! @param EXP_SIZE the expected size of the alignment hit
    //! @param EXP_SCORE the expected score of the alignment hit using the assignment score
    //! @return true if all expected values agree with their the actual values
    bool ConditionCheck
    (
      const align::AlignmentHit< biol::AABase> &ALIGNMENT,
      const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE,
      const bool EXP_EMPTY,
      const size_t EXP_SIZE,
      const double EXP_SCORE
    ) const
    {
      return ALIGNMENT.IsEmpty() == EXP_EMPTY
        && ALIGNMENT.GetSize() == EXP_SIZE
        && math::EqualWithinTolerance( EXP_SCORE, ALIGNMENT.Score( ASSIGN_SCORE));
    }

    //! @brief generates a string for disagreeing values of the alignment word with the given expected values
    //! @param ALIGNMENT the alignment hit to be used
    //! @param ASSIGN_SCORE the assignment score
    //! @param EXP_EMPTY the expected value if the alignment hit should be empty
    //! @param EXP_SIZE the expected size of the alignment hit
    //! @param EXP_SCORE the expected score of the alignment hit using the assignment score
    //! @return a string describing the disagreement
    const std::string ConditionString
    (
      const align::AlignmentHit< biol::AABase> &ALIGNMENT,
      const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE,
      const bool EXP_EMPTY,
      const size_t EXP_SIZE,
      const double EXP_SCORE
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
      return result;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignAlignmentHit

  const ExampleClass::EnumType ExampleAlignAlignmentHit::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignAlignmentHit())
  );

} // namespace bcl
