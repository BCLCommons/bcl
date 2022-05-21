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
#include "score/bcl_score_protein_model_sse_linear_loop_proximity.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "score/bcl_score_protein_model_sse_neighbors.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_sse_linear_loop_proximity.cpp
  //!
  //! @author mendenjl
  //! @date Feb 07, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinModelSSELinearLoopProximity :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinModelSSELinearLoopProximity *Clone() const
    {
      return new ExampleScoreProteinModelSSELinearLoopProximity( *this);
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
      const std::string native_pdb_filename( AddExampleInputPathToFilename( e_Biology, "T0705_native.pdb"));
      const std::string folded_pdb_filename( AddExampleInputPathToFilename( e_Biology, "T0705_bad_loops.pdb"));

      // get the protein model
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( native_pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );
      // get the protein model
      assemble::ProteinModel bad_loops_model
      (
        Proteins::GetModel( folded_pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test constructor from nr_excluded and width and scheme");
      score::ProteinModelSSELinearLoopProximity loop_norm_proximity( true, false);

      BCL_MessageStd( "test constructor from nr_excluded and width and scheme");
      score::ProteinModelSSELinearLoopProximity loop_dist_along_sse_unnorm( false, true);

      BCL_MessageStd( "test clone");
      util::ShPtr< score::ProteinModelSSELinearLoopProximity> sp_loop( loop_norm_proximity.Clone());
      BCL_ExampleCheck( sp_loop->GetScheme(), loop_norm_proximity.GetScheme());

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // initialize two loop sequence and euclidean distance pairs
      const storage::Pair< size_t, double> seq_euc_a( 4, 20);
      const storage::Pair< size_t, double> seq_euc_b( 4, 8);

      // test with an unclosable loop
      BCL_ExampleCheckWithinAbsTolerance( loop_norm_proximity( native_model), 0.65, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( loop_norm_proximity( bad_loops_model), 9.15, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( loop_dist_along_sse_unnorm( native_model), 0.797, 0.01);
      BCL_ExampleCheckWithinAbsTolerance( loop_dist_along_sse_unnorm( bad_loops_model), 9.905, 0.01);

      loop_dist_along_sse_unnorm.WriteDetailedSchemeAndValues( native_model, util::GetLogger());
      loop_dist_along_sse_unnorm.WriteDetailedSchemeAndValues( bad_loops_model, util::GetLogger());

      util::Stopwatch timer( "Time to score 1000 times", util::Time( 10, 0), util::Message::e_Standard);
      for( size_t i( 0), n_times( 1000); i < n_times; ++i)
      {
        loop_dist_along_sse_unnorm( bad_loops_model);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinModelSSELinearLoopProximity

  const ExampleClass::EnumType ExampleScoreProteinModelSSELinearLoopProximity::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelSSELinearLoopProximity())
  );

} // namespace bcl
