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
#include "align/bcl_align_aligner_progressive.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_blosum.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_aligner_progressive.cpp
  //!
  //! @author heinzes1
  //! @date Jan 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignerProgressive :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignerProgressive *Clone() const
    {
      return new ExampleAlignAlignerProgressive( *this);
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

      // read seq_a and seq_b from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test1.fasta"));
      util::ShPtr< biol::AASequence> seq_a
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read))
      );
      io::File::CloseClearFStream( read);
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test2.fasta"));
      util::ShPtr< biol::AASequence> seq_b
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read))
      );
      io::File::CloseClearFStream( read);

      // create alignment
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_a( new align::AlignmentLeaf< biol::AABase>( seq_a));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_b( new align::AlignmentLeaf< biol::AABase>( seq_b));
      // alignment_list
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignment_list;
      alignment_list.PushBack( alignment_a);
      alignment_list.PushBack( alignment_b);
      alignment_list.PushBack( alignment_a);

      // create score
      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );

      // construct a pairwise aligner
      util::ShPtr< align::PairwiseAlignerInterface< biol::AABase> >
        sp_pairwise_aligner( new align::AlignerDynamicProgramming< biol::AABase>());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      align::AlignerProgressive< biol::AABase> default_aligner;

      // test constructor taking pairwise aligner
      align::AlignerProgressive< biol::AABase> aligner( sp_pairwise_aligner, assign_score);

      // test Clone()
      util::ShPtr< align::MultipleAlignerInterface< biol::AABase> > aligner_clone( default_aligner.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test AlignMultiple()
      BCL_MessageStd( "1: aligner.AlignMultiple()");
      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_pair( aligner.AlignMultiple( alignment_list));
      align::AlignmentNode< biol::AABase> result_alignment( alignment_pair.First());
      double result_score( alignment_pair.Second());
      BCL_ExampleIndirectCheck
      (
        result_alignment.GetDepth(),
        3,
        "1a: result_alignment.GetDepth==" + util::Format()( result_alignment.GetDepth()) + "!=3"
      );
      ;
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 3.8, result_score),
        true,
        "1b: result_score==" + util::Format()( result_score) + "!=3.8"
      );

      // test AlignMultiple() with cloned default constructed aligner
      BCL_MessageStd( "2: aligner.AlignMultiple() with cloned aligner");
      alignment_pair = aligner_clone->AlignMultiple( alignment_list);
      result_alignment = alignment_pair.First();
      result_score = alignment_pair.Second();
      BCL_ExampleIndirectCheck
      (
        result_alignment.GetDepth(),
        3,
        "2a: result_alignment.GetDepth==" + util::Format()( result_alignment.GetDepth()) + "!=2"
      );
      ;
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( 0.0, result_score),
        true,
        "2b: result_score==" + util::Format()( result_score) + "!=0.0"
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignAlignerProgressive

  const ExampleClass::EnumType ExampleAlignAlignerProgressive::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignAlignerProgressive())
  );
  
} // namespace bcl
