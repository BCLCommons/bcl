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
#include "score/bcl_score_alignment_assignment.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_handler_pir.h"
//#include "align/bcl_align_factory.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_alignment_assignment.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAlignmentAssignment :
    public ExampleInterface
  {
  public:

      ExampleScoreAlignmentAssignment *Clone() const
    {
      return new ExampleScoreAlignmentAssignment( *this);
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
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // create align::Factory which will handle the amino acid sequences and create simple test alignment
//      align::Factory< biol::AABase> align_factory;
//      align::AlignmentSimple< biol::AABase> alignment1( align_factory.CopySequenceToAlignment( seq));
//      align::AlignmentSimple< biol::AABase> alignment2( align_factory.CopySequenceToAlignment( seq));
//      align::AlignmentSimple< biol::AABase> alignment( align_factory.MergeAlignments( alignment1, alignment2));
//
//      // create an AlignmentHandlerPIR and output the Alignment in pir format
//      align::HandlerPIR< biol::AABase> alignment_handler_pir( 50);
//      alignment_handler_pir.WriteAlignment( util::GetLogger(), alignment);
//
//      // create assignment score
//      util::ShPtr< math::FunctionInterfaceSerializable< align::Assignment< biol::AABase>, double> > sp_assign_score
//      (
//        new score::Assignment< biol::AABase>( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_100))
//      );
//
//    //////////////////////////////////
//    // construction and destruction //
//    //////////////////////////////////
//
//      // create alignment score to test constructor
//      util::ShPtr< score::AlignmentAssignment< biol::AABase> >
//        sp_align_score( new score::AlignmentAssignment< biol::AABase>( sp_assign_score));
//
//    /////////////////
//    // data access //
//    /////////////////
//
//    ///////////////
//    // operators //
//    ///////////////
//
//      // score alignment and print score
//      double expected_score( 78.3);
//      double score( sp_align_score->operator()( alignment));
//      BCL_MessageStd( "Score=" + util::Format()( score));
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        score == expected_score,
//        "Score=" + util::Format()( score) + " Expected: Score=" + util::Format()( expected_score)
//      );

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

  }; //end ExampleScoreAlignmentVertically

  const ExampleClass::EnumType ExampleScoreAlignmentAssignment::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAlignmentAssignment())
  );

} // namespace bcl
