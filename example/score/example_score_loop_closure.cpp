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
#include "score/bcl_score_loop_closure.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "score/bcl_score_protein_model_sse_neighbors.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_loop_closure.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreLoopClosure :
    public ExampleInterface
  {
  public:

    ExampleScoreLoopClosure *Clone() const
    {
      return new ExampleScoreLoopClosure( *this);
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
      // construct min_sse_sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 0;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 0;
      min_sse_sizes[ biol::GetSSTypes().COIL] = 999;

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_bad_loops.pdb"));

      // get the protein model
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // locate helices
      util::SiPtr< const assemble::SSE> strand_a( assemble::LocatorSSE( 'A', 1, 7).Locate( model));
      util::SiPtr< const assemble::SSE> strand_b( assemble::LocatorSSE( 'A', 10, 17).Locate( model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test constructor from nr_excluded and width");
      score::LoopClosure loop( 1, 1.0, 1.0);

      BCL_MessageStd( "test constructor from nr_excluded and width");
      score::LoopClosure loop2( 2, 1.0, 1.0);

      BCL_MessageStd( "test copy constructor");
      score::LoopClosure loop_copy( loop);

      BCL_MessageStd( "test clone");
      util::ShPtr< score::LoopClosure> sp_loop( loop.Clone());
      BCL_Example_Check
      (
        sp_loop->GetScheme() == loop.GetScheme(),
        "Cloned object should be \n" + util::Format()( loop) + "\nnot\n" + util::Format()( *sp_loop)
      );

      score::ProteinModelSSENeighbors score_loop_model( sp_loop, false);

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( loop.GetClassIdentifier(), "bcl::score::LoopClosure");

      // test GetNrExcludedResidues
      BCL_ExampleCheck( loop2.GetNrExcludedResidues(), 2);

      // test GetSigmoidWidth
      BCL_MessageStd( "test GetSigmoidWidth()");
      BCL_Example_Check
      (
        loop.GetSigmoidWidth() == 1.0,
        "GetSigmoidWidth() should return 1.0 not " + util::Format()( loop.GetSigmoidWidth())
      );

    ////////////////
    // operations //
    ////////////////

      // initialize two loop sequence and euclidean distance pairs
      const storage::Pair< size_t, double> seq_euc_a( 4, 20);
      const storage::Pair< size_t, double> seq_euc_b( 4, 8);

      BCL_MessageStd( "testing ScoreLoop()");

      // test with an unclosable loop
      BCL_Example_Check
      (
        loop.ScoreLoop( seq_euc_a) == 1.0,
        "ScoreLoop() should have returned 1.0 but instead it returned false 4 residues and 20 Angstroms!"
      );

      // test with a closable loop
      BCL_Example_Check
      (
        loop.ScoreLoop( seq_euc_b) == 0.0,
        "ScoreLoop() should have returned 0.0 but instead it returned true 4 residues and 8 Angstroms!"
      );

      BCL_MessageStd( "testing SequenceAndEuclideanDistanceWithExclusion()");
      // initialize expected variables
      storage::Pair< size_t, double> expected_seq_euc_excl_0( 2, 33.4751);
      storage::Pair< size_t, double> expected_seq_euc_excl_1( 4, 27.3127);
      storage::Pair< size_t, double> expected_seq_euc_excl_2( 6, 21.8749);

      BCL_MessageStd( "testing with exclusion 0");
      // calculate distances
      storage::Pair< size_t, double> seq_euc_excl_0
      (
        loop.SequenceAndEuclideanDistanceWithExclusion( *strand_a, *strand_b, 0)
      );
      // check the values
      BCL_Example_Check
      (
        seq_euc_excl_0.First() == expected_seq_euc_excl_0.First() &&
        math::EqualWithinTolerance( seq_euc_excl_0.Second(), expected_seq_euc_excl_0.Second()),
        "The distances with exclusion 0 should be\n" + util::Format()( expected_seq_euc_excl_0) + "\nnot\n" +
        util::Format()( seq_euc_excl_0)
      );

      BCL_MessageStd( "testing with exclusion 1");
      // calculate distances
      storage::Pair< size_t, double> seq_euc_excl_1
      (
        loop.SequenceAndEuclideanDistanceWithExclusion( *strand_a, *strand_b, 1)
      );
      // check the values
      BCL_Example_Check
      (
        seq_euc_excl_1.First() == expected_seq_euc_excl_1.First() &&
        math::EqualWithinTolerance( seq_euc_excl_1.Second(), expected_seq_euc_excl_1.Second()),
        "The distances with exclusion 1 should be\n" + util::Format()( expected_seq_euc_excl_1) + "\nnot\n" +
        util::Format()( seq_euc_excl_1)
      );

      BCL_MessageStd( "testing with exclusion 2");
      // calculate distances
      storage::Pair< size_t, double> seq_euc_excl_2
      (
        loop.SequenceAndEuclideanDistanceWithExclusion( *strand_a, *strand_b, 2)
      );
      // check the values
      BCL_Example_Check
      (
        seq_euc_excl_2.First() == expected_seq_euc_excl_2.First() &&
        math::EqualWithinTolerance( seq_euc_excl_2.Second(), expected_seq_euc_excl_2.Second()),
        "The distances with exclusion 2 should be\n" + util::Format()( expected_seq_euc_excl_2) + "\nnot\n" +
        util::Format()( seq_euc_excl_2)
      );

    ///////////////
    // operators //
    ///////////////

      // test the operator()
      BCL_MessageStd( "testing operator()");

      // calculate the score for the model
      const double expected_score( 2);
      const double score( score_loop_model( model));

      // output the score
      BCL_MessageStd( "The loop score for the model is : " + util::Format()( score));

      // compare the score with the expected score
      BCL_Example_Check
      (
        expected_score == score,
        "The score should have been " + util::Format()( expected_score)
      );

      storage::Pair< size_t, double> seq_euc_c( 3, 9.5);
      // do iterative translation of "ca_cb_aa_b" in the z-direction and score its interaction with "ca_cb_aa_a"
      for( size_t i( 0); i < 15; ++i)
      {
        // increment the loop distance
        seq_euc_c.Second() += 0.1;

        const double loop_score( loop.ScoreLoop( seq_euc_c));

        // message score
        BCL_MessageStd
        (
          "dist:\t" + util::Format()( seq_euc_c.Second()) + "\tscore:\t" + util::Format()( loop_score)
        );

        if( i == 6)
        {
          // check the operator results
          const double expected_value( 0.351029);
          BCL_Example_Check
          (
            math::EqualWithinTolerance
            (
              loop_score,
              expected_value
            ),
            "The LoopClosure score for this loop should be " + util::Format()( expected_value) + " but is " +
            util::Format()( loop_score)
          );
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "testing read write")
      WriteBCLObject( loop);
      score::LoopClosure loop_read( util::GetUndefined< size_t>(), util::GetUndefined< double>(), util::GetUndefined< double>());
      ReadBCLObject( loop_read);
      BCL_ExampleCheck( loop_read.GetLabel(), loop.GetLabel());

      // test WriteDetailedSchemeAndValues()
      BCL_MessageStd( "testing WriteDetailedSchemeAndValues()")
      score_loop_model.WriteDetailedSchemeAndValues( model, util::GetLogger());

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreLoopClosure

  const ExampleClass::EnumType ExampleScoreLoopClosure::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreLoopClosure())
  );

} // namespace bcl
