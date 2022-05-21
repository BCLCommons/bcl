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
#include "score/bcl_score_assignment.h"

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
  //! @example example_score_assignment.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAssignment :
    public ExampleInterface
  {
  public:

    ExampleScoreAssignment *Clone() const
    {
      return new ExampleScoreAssignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
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
//      align::AlignmentSimple< biol::AABase> alignment( align_factory.MergeAlignments( alignment1, alignment2));
//
//      // create an AlignmentHandlerPIR and output the Alignment in pir format
//      align::HandlerPIR< biol::AABase> alignment_handler_pir( 50);
//      alignment_handler_pir.WriteAlignment( util::GetLogger(), alignment);
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // create assignment score: test default constructor
//      util::ShPtr< math::FunctionInterfaceSerializable< align::Assignment< biol::AABase>, double> > sp_assign_score
//      (
//        new score::Assignment< biol::AABase>( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160))
//      );
//
//      // create old assignment score for comparison
//      util::ShPtr< score::AssignmentWithGap< biol::AABase> > sp_assign_score_old
//      (
//        new score::AssignmentWithGap< biol::AABase>
//        (
//          score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160),
//          -0.5, -0.1, -0.2, -0.1
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
//      // result string: BCL_Example_Check cannot be used in a loop, multiple checks with the same name are not supported
//      std::string result_string;
//
//      // test operator on each assignment and print out score
//      for
//      (
//        align::AlignmentSimple< biol::AABase>::const_iterator itr( alignment.GetAssignments().Begin()), itr_end( alignment.GetAssignments().End());
//        itr != itr_end;
//        ++itr
//      )
//      {
//        // calculate score using the new and old method, and output assignment and its score
//        double score( sp_assign_score->operator()( **itr));
//        double score_old( sp_assign_score_old->operator()( **itr));
//        util::GetLogger() << "Assignment=" << WriteShortAssignmentInfo( **itr) << " Score=" << score << "|"
//          << score_old << '\n';
//
//        // try to find gaps in current assignment
//        size_t gap_count( 0);
//        for
//        (
//          util::SiPtrList< const biol::AABase>::const_iterator
//            itr_members( ( *itr)->GetMembers().Begin()), itr_members_end( ( *itr)->GetMembers().End());
//          itr_members != itr_members_end;
//          ++itr_members
//        )
//        {
//          if( !( *itr_members).IsDefined()) // if a pointer is not defined, it's a gap
//          {
//            gap_count++;
//            break;
//          }
//        }
//
//        // ensure new and old score are equal, only if assignment has no gaps; gap scores are differently implemented
//        if( score != score_old && gap_count == 0)
//        {
//          result_string.append
//          (
//            "Assignment " + WriteShortAssignmentInfo( **itr) + " score differs: " + util::Format()( score) + "!="
//              + util::Format()( score_old) + "\n"
//          );
//        }
//      }
//
//      // check result
//      BCL_Example_Check( result_string.length() == 0, result_string);

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

  private:

    // write one letter codes of AABases in given Assignment
    std::string WriteShortAssignmentInfo( const align::Assignment< biol::AABase> &ASSIGNMENT) const
    {
      std::string result;
      util::SiPtrList< const biol::AABase> members( ASSIGNMENT.GetMembers());

      // loop over members of the assignment (aabases) and print each one letter code
      for
      (
        util::SiPtrList< const biol::AABase>::const_iterator itr( members.Begin()), itr_end( members.End());
        itr != itr_end;
        ++itr
      )
      {
        // append gap char if SiPtr is empty (!IsDefined)
        result.append( 1, ( *itr).IsDefined() ? ( *itr)->GetType()->GetOneLetterCode() : '_');
      }

      return result;
    }

  public:

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAssignment

  const ExampleClass::EnumType ExampleScoreAssignment::s_Instance( GetExamples().AddEnum( ExampleScoreAssignment()));

} // namespace bcl
