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
#include "score/bcl_score_aa_pair_atom_clash.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_pair_atom_clash.cpp
  //!
  //! @author karakam, alexanns
  //! @date Jun 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAPairAtomClash :
    public ExampleInterface
  {
  public:

    ExampleScoreAAPairAtomClash *Clone() const
    {
      return new ExampleScoreAAPairAtomClash( *this);
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

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_internal_clash.pdb"));
      // get the protein model
      util::ShPtr< assemble::SSE> sp_loop( Proteins::GetSSE( pdb_filename, 'A', 46, 63));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::AAPairAtomClash score_def;

      // construct from variables
      score::AAPairAtomClash score( 1.0, 1, "test");

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( score_def.GetSigmoidWidth(), 1.0);
      BCL_ExampleCheck( score_def.GetMinimalSequenceSeparation(), 0);
      BCL_ExampleCheck( score_def.GetDistanceCutoff(), 12.0);
      BCL_ExampleCheck( score_def.GetConsiderDifferentChain(), true);

      BCL_ExampleCheck( score.GetSigmoidWidth(), 1.0);
      BCL_ExampleCheck( score.GetMinimalSequenceSeparation(), 1);
      BCL_ExampleCheck( score.GetDistanceCutoff(), 12.0);
      BCL_ExampleCheck( score.GetConsiderDifferentChain(), true);

    ///////////////
    // operators //
    ///////////////

      BCL_MessageStd( "Scoring loop " + sp_loop->GetIdentification());

      // score
      double sum_score( 0.0);

      // iterate over all residues on
      for
      (
        biol::AASequence::const_iterator aa_itr_a( sp_loop->Begin()), aa_itr_a_end( sp_loop->End() - 2);
        aa_itr_a != aa_itr_a_end; ++aa_itr_a
      )
      {
        for
        (
          biol::AASequence::const_iterator aa_itr_b( aa_itr_a + 2), aa_itr_b_end( sp_loop->End());
          aa_itr_b != aa_itr_b_end; ++aa_itr_b
        )
        {
          // score this pair
          const double this_score( score( **aa_itr_a, **aa_itr_b));
          sum_score += this_score;

          BCL_MessageDbg
          (
            ( *aa_itr_a)->GetIdentification() + " -- " + ( *aa_itr_b)->GetIdentification() +
            " " + util::Format()( this_score)
          );
        }
      }

      // test with expected score
      BCL_MessageStd( "sum_score: " + util::Format()( sum_score));
      const double expected_sum_score( 2.09239);
      BCL_ExampleCheckWithinTolerance( expected_sum_score, sum_score, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // write and read object back in
      WriteBCLObject( score);
      score::AAPairAtomClash score_read;
      ReadBCLObject( score_read);

      // test that reading was correct
      BCL_ExampleCheck( score.GetSigmoidWidth(), score_read.GetSigmoidWidth());
      BCL_ExampleCheck( score.GetMinimalSequenceSeparation(), score_read.GetMinimalSequenceSeparation());
      BCL_ExampleCheck( score.GetDistanceCutoff(), score_read.GetDistanceCutoff());
      BCL_ExampleCheck( score.GetConsiderDifferentChain(), score_read.GetConsiderDifferentChain());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAPairAtomClash

  const ExampleClass::EnumType ExampleScoreAAPairAtomClash::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAPairAtomClash())
  );

} // namespace bcl
