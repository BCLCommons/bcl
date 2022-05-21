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
#include "align/bcl_align_aligner_shift.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "quality/bcl_quality_lcs.h"
#include "score/bcl_score_alignment_quality.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_aligner_shift.cpp
  //! @brief test the aligner that creates alignments by simply shifting the assignments systematically
  //!
  //! @author woetzen
  //! @date Mar 22, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignAlignerShift :
    public ExampleInterface
  {
  public:

    ExampleAlignAlignerShift *Clone() const
    {
      return new ExampleAlignAlignerShift( *this);
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

    /////////////////
    // preparation //
    /////////////////

      util::ShPtr< quality::MeasureInterface> sp_lcs( new quality::LCS( 3.0, 3));

      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // get model and sse
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));
      util::ShPtr< assemble::SSE> sp_sse_helix_23_34( Proteins::GetSSE( pdb_filename, 'A', 23, 34));
      util::ShPtr< assemble::SSE> sp_sse_strand_40_45( Proteins::GetSSE( pdb_filename, 'A', 40, 45));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      align::AlignerShift< biol::AABase> aligner_shift;

    /////////////////
    // data access //
    /////////////////

      // for bid alignments of the form:
      // ABC--
      // --abc
      // set the minimum number of assignments without a gap
      aligner_shift.SetMinNumberAssignmentsWithoutGap( 3);

      // forbid alignments of the form:
      // ABCD--
      // --ABCD
      // set the maximum number of gaps (from one side)
      aligner_shift.SetMaxNumberGaps( 100, 100);

      // set the scoring function
      aligner_shift.SetScoreComparisonFunction( util::ToSiPtr( sp_lcs->GetComparisonFunction()));

      // set the scoring function
      util::ShPtr< score::AlignmentQuality> sp_score( new score::AlignmentQuality());
      sp_score->SetMeasure( sp_lcs);
      sp_score->SetAtoms( storage::Set< biol::AtomType>( biol::GetAtomTypes().CA));
      aligner_shift.SetScoringFunction( sp_score);

    ////////////////
    // operations //
    ////////////////

      // align
      util::ShPtr< align::AlignmentInterface< biol::AABase> > leaf_seq( new align::AlignmentLeaf< biol::AABase>( model.GetChain( 'A')->GetSequence()));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > leaf_sse_helix( new align::AlignmentLeaf< biol::AABase>( sp_sse_helix_23_34));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > leaf_sse_strand( new align::AlignmentLeaf< biol::AABase>( sp_sse_strand_40_45));

      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_all
      (
        aligner_shift.AlignPair( leaf_sse_helix, leaf_seq)
      );

      BCL_ExampleIndirectCheck( alignment_all.Second(), 12.0, "expected score is 12, alignment score is: " + util::Format()( alignment_all.Second()));

      // align pairpair
      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_pairpair
      (
        aligner_shift.AlignPairPair( leaf_sse_helix, leaf_sse_helix, leaf_sse_strand, leaf_sse_strand, 3, 3, true)
      );
      BCL_ExampleIndirectCheck( alignment_pairpair.Second(), 18, "expected score is 17, alignment score is: " + util::Format()( alignment_pairpair.Second()));

      aligner_shift.SetMaxNumberGaps( 3, 7);
      storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment_3_7
      (
        aligner_shift.AlignPair( leaf_sse_helix, leaf_seq)
      );
      BCL_ExampleIndirectCheck( alignment_3_7.First().IsEmpty(), true, "no alignment should be possible with the max number of gaps");

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( aligner_shift);

      // read the object back in
      align::AlignerShift< biol::AABase> aligner_shift_read;
      ReadBCLObject( aligner_shift_read);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAlignAlignerShift

  const ExampleClass::EnumType ExampleAlignAlignerShift::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignAlignerShift())
  );

} // namespace bcl
