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
#include "align/bcl_align_alignment_node.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_leaf.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_alignment_node.cpp
  //!
  //! @author heinzes1
  //! @date 2010/10/29
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignmentNode :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignmentNode *Clone() const
    {
      return new ExampleAlignAlignmentNode( *this);
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

      // create AlignmentLeaf using constructor from SequenceInterface
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_a( new align::AlignmentLeaf< biol::AABase>( fasta_sequence));
      // expected values: parent_alignment_size=0, depth=1, size=0, empty=true
      size_t exp_child_alignment_size( 0);
      size_t exp_depth( 1);
      size_t exp_size( 140);
      bool exp_empty( false);
      BCL_MessageStd( "1: AlignmentLeaf( SequenceInterface) constructor");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "1: " + ConditionString( alignment_a, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // create AlignmentLeaf from sequence_list
      util::ShPtrList< align::SequenceInterface< biol::AABase> > sequence_list( 2, fasta_sequence);
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_b( new align::AlignmentLeaf< biol::AABase>( sequence_list));
      // expected values: parent_alignment_size=0, depth=1, size=0, empty=true
      exp_depth = 2;
      exp_size = 0;
      exp_empty = true;
      BCL_MessageStd( "2: AlignmentLeaf( list of SequenceInterface) constructor");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_b, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "2: " + ConditionString( alignment_b, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // create an incorrect assignment to test, with incorrect number aa in the assignment (assignment.size=0)
      util::ShPtr< align::Assignment< biol::AABase> > assignment_incorrect( new align::Assignment< biol::AABase>());
      // create a correct assignment to test (with two aa "assignments", since this is for a node alignment)
      util::ShPtr< align::Assignment< biol::AABase> > assignment_correct( new align::Assignment< biol::AABase>());
      assignment_correct->Append( fasta_sequence->GetFirstAA());
      assignment_correct->Append( fasta_sequence->GetAA( fasta_sequence->GetSize() / 2));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test AlignmentNode constructor from two Alignments, both being AlignmentLeaf
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_c( new align::AlignmentNode< biol::AABase>( alignment_a, alignment_a));
      BCL_MessageStd( "3: AlignmentNode( AlignmentInterface=Leaf1, AlignmentInterface=Leaf1)");
      exp_child_alignment_size = exp_depth = 2;
      exp_size = 0;
      exp_empty = true;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "3: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from two Alignments with one being AlignmentNode
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_d( new align::AlignmentNode< biol::AABase>( alignment_a, alignment_c));
      BCL_MessageStd( "4: AlignmentNode( AlignmentInterface=Leaf1, AlignmentInterface=Node)");
      exp_depth = 3;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_d, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "4: " + ConditionString( alignment_d, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from two Alignments with one being AlignmentLeaf with multiple sequences
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_e( new align::AlignmentNode< biol::AABase>( alignment_b, alignment_c));
      BCL_MessageStd( "5: AlignmentNode( AlignmentInterface=Leaf2, AlignmentInterface=Node)");
      exp_depth = 4;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_e, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "5: " + ConditionString( alignment_e, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from list of parent alignments with 4 alignments
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > parent_alignment_list;
      parent_alignment_list.Append( alignment_c->GetChildAlignments()); // alignment_c has two alignment_a
      parent_alignment_list.Append( alignment_c->GetChildAlignments());
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_f( new align::AlignmentNode< biol::AABase>( parent_alignment_list));
      BCL_MessageStd( "6: AlignmentNode( ShPtrList< AlignmentInterface> correct)");
      exp_child_alignment_size = 4;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_f, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "6: " + ConditionString( alignment_f, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from list of parent alignments with one alignment, i.e. incorrect list
      parent_alignment_list.Reset();
      parent_alignment_list.PushBack( alignment_a);
      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_g( new align::AlignmentNode< biol::AABase>( parent_alignment_list));
      BCL_MessageStd( "7: AlignmentNode( ShPtrList< AlignmentInterface> incorrect)");
      exp_child_alignment_size = exp_depth = 0;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_g, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "7: " + ConditionString( alignment_g, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode.Clone
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_h( alignment_c->Clone());
      BCL_MessageStd( "8: AlignmentNode.Clone()");
      exp_child_alignment_size = exp_depth = 2;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_h, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "8: " + ConditionString( alignment_h, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // test Prepend with incorrect assignment
      bool success( alignment_c->Prepend( assignment_incorrect));
      std::string success_str( "9: Prepend.success==" + util::Format()( success));
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, false, success_str + "!=false");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "9: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Append with incorrect assignment
      success = alignment_c->Append( assignment_incorrect);
      success_str = "10: Append.success==" + util::Format()( success);
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, false, success_str + "!=false");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "10: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Prepend with correct assignment
      success = alignment_c->Prepend( assignment_correct);
      success_str = "11: Prepend.success==" + util::Format()( success);
      exp_size++;
      exp_empty = false;
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck( success, true, success_str + "!=true");
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "11: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Append with correct assignment
      exp_size++;
      BCL_ExampleCheck( alignment_c->Append( assignment_correct), true);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "12: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Prepend( assignment_list)
      BCL_ExampleCheck( alignment_h->Prepend( alignment_c->GetAssignments()), true);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_h, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "13: " + ConditionString( alignment_h, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test ResetAssignments
      alignment_c->ResetAssignments();
      success_str = "14: ResetAssignments()";
      exp_size = 0;
      exp_empty = true;
      BCL_MessageStd( success_str);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "14: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test Append( assignment_list)
      exp_size = 2;
      exp_empty = false;
      BCL_MessageStd( success_str);
      BCL_ExampleCheck( alignment_c->Append( alignment_h->GetAssignments()), true);
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "15: " + ConditionString( alignment_c, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentLeaf constructor from AlignmentInterface with AlignmentLeaf is in ExampleAlignmentLeaf
      // test AlignmentLeaf constructor from AlignmentInterface with AlignmentNode
      storage::Pair< align::AlignmentLeaf< biol::AABase>, int> alignment_pair( *alignment_c, 42); // call constructor
      util::ShPtr< align::AlignmentLeaf< biol::AABase> > alignment_i( alignment_pair.First().Clone());
      BCL_MessageStd( "16: AlignmentLeaf( AlignmentInterface=Node) ");
      exp_child_alignment_size = 0; // collapsed leaf from a node
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_i, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "16: " + ConditionString( alignment_i, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from AlignmentInterface with AlignmentLeaf; NOT FULLY IMPLEMENTED
      storage::Pair< align::AlignmentNode< biol::AABase>, int> alignment_pair_a( *alignment_a, 33); // call constructor
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment_j( alignment_pair_a.First().Clone());
      BCL_MessageStd( "17: AlignmentNode( AlignmentInterface=Leaf1) ");
      exp_depth = exp_size = 0;
      exp_empty = true;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_j, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "17: " + ConditionString( alignment_j, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from AlignmentInterface with AlignmentLeaf; NOT FULLY IMPLEMENTED
      storage::Pair< align::AlignmentNode< biol::AABase>, int> alignment_pair_b( *alignment_b, 33); // call constructor
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment_k( alignment_pair_b.First().Clone());
      BCL_MessageStd( "18: AlignmentNode( AlignmentInterface=Leaf2) ");
      exp_depth = exp_size = 0;
      exp_empty = true;
      BCL_ExampleIndirectCheck
      (
        ConditionCheck( alignment_k, exp_child_alignment_size, exp_depth, exp_size, exp_empty),
        true,
        "18: " + ConditionString( alignment_k, exp_child_alignment_size, exp_depth, exp_size, exp_empty)
      );

      // test AlignmentNode constructor from AlignmentInterface with AlignmentNode

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

  }; //end ExampleAlignmentNode

  const ExampleClass::EnumType ExampleAlignAlignmentNode::s_Instance( GetExamples().AddEnum( ExampleAlignAlignmentNode()));

} // namespace bcl
