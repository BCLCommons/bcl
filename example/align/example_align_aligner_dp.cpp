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
#include "align/bcl_align_aligner_dp.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_progressive.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_blast_profile.h"
#include "score/bcl_score_aa_assignment_blosum.h"
#include "score/bcl_score_aa_assignment_pam.h"
#include "score/bcl_score_aa_assignment_ss_prediction.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_aligner_dp.cpp
  //!
  //! @author heinzes1
  //! @date Jan 20, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignerDP :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignerDP *Clone() const
    {
      return new ExampleAlignAlignerDP( *this);
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

      align::AlignerDP< biol::AABase> aligner_default;

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // test aligning short sequences
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test3.fasta"));
      util::ShPtr< biol::AASequence>
        seq_a( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      util::ShPtr< align::AlignmentInterface< biol::AABase> >
        alignment_a( new align::AlignmentLeaf< biol::AABase>( seq_a));

      score::AssignmentWithGap< biol::AABase> assign_score
      (
        util::CloneToShPtr( double( 1.0) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -1.0, -1.0, -1.0, -1.0
      );

      aligner_default.SetScoringFunction( assign_score);
      storage::Pair< align::AlignmentNode< biol::AABase>, double>
        alignment( aligner_default.AlignPair( alignment_a, alignment_a));

      align::HandlerPIR< biol::AABase> handler;
      handler.WriteAlignment( util::GetLogger(), alignment.First());

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // do MSA tests and check the resulting score is still correct

      // define scoring function
      score::AssignmentWithGap< biol::AABase> test_assign_score
      (
        util::CloneToShPtr( double( 0.6) * score::AAAssignmentBlastProfile()
          + double( 0.2) * score::AAAssignmentSSPrediction( sspred::GetMethods().e_PSIPRED)
          + double( 0.1) * score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160)
          + double( 0.1) * score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45)),
        -0.5,
        -0.1,
        -0.2,
        -0.1
      );

      // create progressive and dynamic programming aligner which will perform the alignment
      util::ShPtr< align::PairwiseAlignerInterface< biol::AABase> >
        aligner_dynamic_programming( new align::AlignerDP< biol::AABase>());
      align::AlignerProgressive< biol::AABase> aligner_progressive( aligner_dynamic_programming, test_assign_score);

      // first test: seq_c and seq_d

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.fasta"));
      util::ShPtr< biol::AASequence> seq_c( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.psipred_ss2"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *seq_c, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *seq_c);
      io::File::CloseClearFStream( read);

      // read from fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.fasta"));
      util::ShPtr< biol::AASequence> seq_d( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.psipred_ss2"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *seq_d, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1f5mA.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *seq_d);
      io::File::CloseClearFStream( read);

      // create single-sequence alignments and alignment list
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_c( new align::AlignmentLeaf< biol::AABase>( seq_c));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_d( new align::AlignmentLeaf< biol::AABase>( seq_d));
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignment_list;
      alignment_list.PushBack( alignment_c);
      alignment_list.PushBack( alignment_d);

      // align sequences
      storage::Pair< align::AlignmentNode< biol::AABase>, double>
        alignment_and_score( aligner_progressive.AlignMultiple( alignment_list));

      // check for correct score
      double expected_score( 4.7275);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( expected_score, alignment_and_score.Second()),
        true,
        "alignment score is " + util::Format()( alignment_and_score.Second()) + ", should be "
          + util::Format()( expected_score)
      );

      // second test: seq_e and seq_f

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.fasta"));
      util::ShPtr< biol::AASequence> seq_e( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.psipred_ss2"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *seq_e, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1eco_.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *seq_e);
      io::File::CloseClearFStream( read);

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1IE9.fasta"));
      util::ShPtr< biol::AASequence> seq_f( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1IE9.psipred_ss2"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *seq_f, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1IE9.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *seq_f);
      io::File::CloseClearFStream( read);

      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_e( new align::AlignmentLeaf< biol::AABase>( seq_e));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_f( new align::AlignmentLeaf< biol::AABase>( seq_f));
      alignment_list.Reset();
      alignment_list.PushBack( alignment_e);
      alignment_list.PushBack( alignment_f);

      // align sequences
      alignment_and_score = aligner_progressive.AlignMultiple( alignment_list);

      // check for correct score
      expected_score = -1.73179;
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( expected_score, alignment_and_score.Second()),
        true,
        "alignment score is " + util::Format()( alignment_and_score.Second()) + ", should be "
          + util::Format()( expected_score)
      );

      // third test: seq_g and seq_h

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1HN0.fasta"));
      util::ShPtr< biol::AASequence> seq_g( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1HN0.psipred_ss2"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *seq_g, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1HN0.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *seq_g);
      io::File::CloseClearFStream( read);

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1EPW.fasta"));
      util::ShPtr< biol::AASequence> seq_h( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read psipred ss prediction
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1EPW.psipred_ss2"));
      sspred::MethodHandler::ReadPredictionsForAASequence( read, *seq_h, sspred::GetMethods().e_PSIPRED);
      io::File::CloseClearFStream( read);

      // read blast profile
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1EPW.ascii"));
      biol::BlastProfileHandler::ReadProfileForAASequence( read, *seq_h);
      io::File::CloseClearFStream( read);

      // create single-sequence alignments and alignment list
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_g( new align::AlignmentLeaf< biol::AABase>( seq_g));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_h( new align::AlignmentLeaf< biol::AABase>( seq_h));
      alignment_list.Reset();
      alignment_list.PushBack( alignment_g);
      alignment_list.PushBack( alignment_h);

      // check for correct score
      BCL_ExampleCheckWithinTolerance( 2.55857, aligner_progressive.AlignMultiple( alignment_list).Second(), 0.0001);

      // fourth test: use three sequences

      // set simple scoring function
      aligner_progressive.SetScoringFunction( assign_score);

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test1.fasta"));
      util::ShPtr< biol::AASequence> seq_i( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // read fasta
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "1fms_test2.fasta"));
      util::ShPtr< biol::AASequence> seq_k( new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read)));
      io::File::CloseClearFStream( read);

      // create single-sequence alignments and alignment list
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_i( new align::AlignmentLeaf< biol::AABase>( seq_i));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > alignment_k( new align::AlignmentLeaf< biol::AABase>( seq_k));
      alignment_list.Reset();
      alignment_list.PushBack( alignment_i);
      alignment_list.PushBack( alignment_i); // use this alignment twice
      alignment_list.PushBack( alignment_k);

      // check for correct score
      BCL_ExampleCheckWithinTolerance( 3.8, aligner_progressive.AlignMultiple( alignment_list).Second(), 0.0001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignAlignerDP

  const ExampleClass::EnumType ExampleAlignAlignerDP::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignAlignerDP())
  );

} // namespace bcl
