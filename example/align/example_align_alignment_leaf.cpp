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
#include "align/bcl_align_alignment_leaf.h"

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
  //! @example example_align_alignment_leaf.cpp
  //!
  //! @author heinzes1
  //! @date 2010/10/29
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignmentLeaf :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignmentLeaf *Clone() const
    {
      return new ExampleAlignAlignmentLeaf( *this);
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
      // create fasta_sequence and in_file_stream and read stream in fasta_sequence
      io::IFStream in_file_stream;
      BCL_ExampleMustOpenInputFile( in_file_stream, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence> fasta_sequence
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( in_file_stream))
      );
      io::File::CloseClearFStream( in_file_stream);

      // create an incorrect assignment to test, with incorrect number aa in the assignment (assignment.size=0)
      util::ShPtr< align::Assignment< biol::AABase> > assignment_incorrect( new align::Assignment< biol::AABase>());
      // create two correct assignments to test (only one aa "assignments", since this is for a leaf alignment)
      util::ShPtr< align::Assignment< biol::AABase> > assignment_a( new align::Assignment< biol::AABase>());
      assignment_a->Append( fasta_sequence->GetFirstAA());
      util::ShPtr< align::Assignment< biol::AABase> > assignment_b( new align::Assignment< biol::AABase>());
      assignment_b->Append( fasta_sequence->GetAA( fasta_sequence->GetSize() / 2));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test AlignmentLeaf( SequenceInterface) constructor
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_a( new align::AlignmentLeaf< biol::AABase>( fasta_sequence));
      // expected values: parent_alignment_size=0, depth=1, size=0, empty=true
      size_t exp_child_alignment_size( 0);
      size_t exp_depth( 1);
      size_t exp_size( 140);
      bool exp_empty( false);
      BCL_MessageStd( "1: AlignmentLeaf( SequenceInterface)");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "1: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Clone
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_b( alignment_a->Clone());
      BCL_MessageStd( "2: AlignmentLeaf.Clone()");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_b, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "2: " + ConditionString( alignment_b, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentLeaf( list of SequenceInterface) constructor
      util::ShPtrList< align::SequenceInterface< biol::AABase> > sequence_list( 1, fasta_sequence);
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_c( new align::AlignmentLeaf< biol::AABase>( sequence_list));
      // expected values: parent_alignment_size=0, depth=1, size=0, empty=true
      exp_depth = 1;
      exp_size = 0;
      exp_empty = true;
      BCL_MessageStd( "3: AlignmentLeaf( list of SequenceInterface)");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "3: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentLeaf( list of SequenceInterface) constructor
      sequence_list.Append( fasta_sequence);
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_d( new align::AlignmentLeaf< biol::AABase>( sequence_list));
      // expected values: parent_alignment_size=0, depth=1, size=0, empty=true
      exp_depth = 2;
      BCL_MessageStd( "4: AlignmentLeaf( list of SequenceInterface)");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_d, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "4: " + ConditionString( alignment_d, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

    /////////////////
    // data access //
    /////////////////

      // tests for GetParentAlignments(), GetDepth(), GetSequences(), GetSize(), GetAssignments(), IsEmpty() are in
      // ConditionCheck() and ConditionString() functions below

    ////////////////
    // operations //
    ////////////////

      // test ResetAssignments
      alignment_a->ResetAssignments();
      std::string success_str( "5: ResetAssignments()");
      exp_depth = 1;
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "5: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Prepend with incorrect assignment
      bool success( alignment_a->Prepend( assignment_incorrect));
      success_str = "6: Prepend.success==" + util::Format()( success);
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, false, success_str + "!=false");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "6: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Append with incorrect assignment, should not change the size either
      success = alignment_a->Append( assignment_incorrect);
      success_str = "7: PushBack.success==" + util::Format()( success);
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, false, success_str + "!=false");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "7: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Prepend with correct assignment
      success = alignment_a->Prepend( assignment_a);
      success_str = "8: Prepend.success==" + util::Format()( success);
      exp_size++;
      exp_empty = false;
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, true, success_str + "!=true");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "8: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Append with correct assignment
      success = alignment_a->Append( assignment_b);
      success_str = "9: Append.success==" + util::Format()( success);
      exp_size++;
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, true, success_str + "!=true");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "9: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test with GetAssignments if assignments are in correct order
      // first aa in alignment should be same as first aa in sequence
      biol::AAType aatype_alignment( ( **alignment_a->GetAssignments().Begin()).GetMembers().FirstElement()->GetType());
      biol::AAType aatype_sequence( fasta_sequence->GetFirstAA()->GetType());
      BCL_MessageStd( "10: alignment.first.aa==" + util::Format()( aatype_alignment));
      BCL_MessageStd( "10: sequence.first.aa==" + util::Format()( aatype_sequence));
      BCL_ExampleIndirectCheck
      (
        aatype_alignment == aatype_sequence,
        true,
        "10: alignment.first.aa==" + util::Format()( aatype_alignment)
          + "!=" + util::Format()( aatype_sequence) + "==sequence.first.aa"
      );
      // second aa in alignment should be same as center aa in sequence
      aatype_alignment = ( **( ++alignment_a->GetAssignments().Begin())).GetMembers().FirstElement()->GetType();
      aatype_sequence = fasta_sequence->GetAA( fasta_sequence->GetSize() / 2)->GetType();
      BCL_MessageStd( "10: alignment.second.aa==" + util::Format()( aatype_alignment));
      BCL_MessageStd( "10: sequence.center.aa==" + util::Format()( aatype_sequence));
      BCL_ExampleIndirectCheck
      (
        aatype_alignment == aatype_sequence,
        true,
        "10: alignment.second.aa==" + util::Format()( aatype_alignment)
          + "!=" + util::Format()( aatype_sequence) + "==sequence.center.aa"
      );

      // test Prepend( assignment_list)
      success = alignment_c->Prepend( alignment_a->GetAssignments());
      success_str = "11: Prepend_list.success==" + util::Format()( success);
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, true, success_str + "!=true");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "11: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test with GetAssignments if assignments are in correct order
      // first aa in alignment should be same as first aa in sequence
      aatype_alignment = ( **alignment_c->GetAssignments().Begin()).GetMembers().FirstElement()->GetType();
      aatype_sequence = fasta_sequence->GetFirstAA()->GetType();
      BCL_MessageStd( "12: alignment.first.aa==" + util::Format()( aatype_alignment));
      BCL_MessageStd( "12: sequence.first.aa==" + util::Format()( aatype_sequence));
      BCL_ExampleIndirectCheck
      (
        aatype_alignment == aatype_sequence,
        true,
        "12: alignment.first.aa==" + util::Format()( aatype_alignment)
          + "!=" + util::Format()( aatype_sequence) + "==sequence.first.aa"
      );
      // second aa in alignment should be same as center aa in sequence
      aatype_alignment = ( **( ++alignment_c->GetAssignments().Begin())).GetMembers().FirstElement()->GetType();
      aatype_sequence = fasta_sequence->GetAA( fasta_sequence->GetSize() / 2)->GetType();
      BCL_MessageStd( "12: alignment.second.aa==" + util::Format()( aatype_alignment));
      BCL_MessageStd( "12: sequence.center.aa==" + util::Format()( aatype_sequence));
      BCL_ExampleIndirectCheck
      (
        aatype_alignment == aatype_sequence,
        true,
        "12: alignment.second.aa==" + util::Format()( aatype_alignment)
          + "!=" + util::Format()( aatype_sequence) + "==sequence.center.aa"
      );

      // test Append( assignment_list)
      success = alignment_a->Append( alignment_c->GetAssignments());
      success_str = "13: Append_list.success==" + util::Format()( success);
      exp_size = 4;
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, true, success_str + "!=true");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "13: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test constructor from AlignmentInterface with AlignmentLeaf
      align::AlignmentInterface< biol::AABase> &alignment_e( *alignment_b); // create reference to an AlignmentInterface
      storage::Pair< align::AlignmentLeaf< biol::AABase>, int> alignment_pair( alignment_e, 42); // call constructor
      util::ShPtr< align::AlignmentLeaf< biol::AABase> > alignment_f( alignment_pair.First().Clone());
      exp_size = 140;
      BCL_MessageStd( "14: AlignmentLeaf( AlignmentInterface) constructor");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_f, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "14: " + ConditionString( alignment_f, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test constructor from AlignmentInterface with AlignmentNode in ExampleAlignmentNode

      // create alignment with two sequences for Score() tests
      align::AlignerMerge< biol::AABase> aligner;
      align::AlignmentNode< biol::AABase> alignment_g( aligner.AlignPair( alignment_a, alignment_a).First());

      // test Score()
      double score( alignment_g.Score( score::AssignmentWithGap< biol::AABase>()));
      double expected_score( 0.0);
      BCL_MessageStd( "15: Score(default_score) = " + util::Format()( score));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( expected_score, score),
        true,
        "15: expected Score(default_score) =" + util::Format()( expected_score)
      );

      // create assignment score
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );

      // test Score with assign_score
      score = alignment_g.Score( assign_score);
      expected_score = 3.2;
      BCL_MessageStd( "16: Score(assign_score) = " + util::Format()( score));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( expected_score, score),
        true,
        "16: expected Score(default_score) =" + util::Format()( expected_score)
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    //! @brief function to check if all sub-conditions are met
    //! @param ALIGNMENT the alignment to check
    //! @param EXP_CHILD_ALIGNMENT_SIZE expected number of parent alignments
    //! @param EXP_DEPTH expected depth of the alignment
    //! @param EXP_SIZE expected size, i.e. number of assignments
    //! @param EXP_EMPTY expected is-empty-value
    //! @return if all sub-conditions are true
    bool ConditionCheck
    (
      const util::ShPtr< align::AlignmentInterface< biol::AABase> > &ALIGNMENT,
      const size_t EXP_CHILD_ALIGNMENT_SIZE,
      const size_t EXP_DEPTH,
      const size_t EXP_SIZE,
      const bool EXP_EMPTY
    ) const
    {
      return ALIGNMENT->GetChildAlignments().GetSize() == EXP_CHILD_ALIGNMENT_SIZE
        && ALIGNMENT->GetDepth() == EXP_DEPTH
        && ALIGNMENT->GetDepth() == ALIGNMENT->GetSequences().GetSize()
        && ALIGNMENT->GetSize() == EXP_SIZE
        && ALIGNMENT->GetSize() == ALIGNMENT->GetAssignments().GetSize()
        && ALIGNMENT->IsEmpty() == EXP_EMPTY;
    }

    //! @brief function assembles an error string based on the sub-condition not fulfilled
    //! @param ALIGNMENT the alignment to check
    //! @param EXP_CHILD_ALIGNMENT_SIZE expected number of parent alignments
    //! @param EXP_DEPTH expected depth of the alignment
    //! @param EXP_SIZE expected size, i.e. number of assignments
    //! @param EXP_EMPTY expected is-empty-value
    //! @return the error string
    const std::string ConditionString
    (
      const util::ShPtr< align::AlignmentInterface< biol::AABase> > &ALIGNMENT,
      const size_t EXP_CHILD_ALIGNMENT_SIZE,
      const size_t EXP_DEPTH,
      const size_t EXP_SIZE,
      const bool EXP_EMPTY
    ) const
    {
      std::string result;
      if( ALIGNMENT->GetChildAlignments().GetSize() != EXP_CHILD_ALIGNMENT_SIZE)
      {
        result.append
        (
          "child_alignments.Size==" + util::Format()( ALIGNMENT->GetChildAlignments().GetSize())
          + "!=" + util::Format()( EXP_CHILD_ALIGNMENT_SIZE) + "; "
        );
      }
      if( ALIGNMENT->GetDepth() != EXP_DEPTH)
      {
        result.append( "Depth==" + util::Format()( ALIGNMENT->GetDepth()) + "!=" + util::Format()( EXP_DEPTH) + "; ");
      }
      if( ALIGNMENT->GetDepth() != ALIGNMENT->GetSequences().GetSize())
      {
        result.append
        (
          "Depth==" + util::Format()( ALIGNMENT->GetDepth())
          + "!=" + util::Format()( ALIGNMENT->GetSequences().GetSize()) + "==sequences.Size; "
        );
      }
      if( ALIGNMENT->GetSize() != EXP_SIZE)
      {
        result.append( "Size==" + util::Format()( ALIGNMENT->GetSize()) + "!=" + util::Format()( EXP_SIZE) + "; ");
      }
      if( ALIGNMENT->GetSize() != ALIGNMENT->GetAssignments().GetSize())
      {
        result.append
        (
          "Size==" + util::Format()( ALIGNMENT->GetSize()) + "!="
          + util::Format()( ALIGNMENT->GetAssignments().GetSize()) + "==GetAssignments.Size; "
        );
      }
      if( ALIGNMENT->IsEmpty() != EXP_EMPTY)
      {
        result.append( "IsEmpty==" + util::Format()( ALIGNMENT->IsEmpty()) + "!=" + util::Format()( EXP_EMPTY) + "; ");
      }
      return result;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignmentLeaf

  const ExampleClass::EnumType ExampleAlignAlignmentLeaf::s_Instance( GetExamples().AddEnum( ExampleAlignAlignmentLeaf()));

} // namespace bcl
