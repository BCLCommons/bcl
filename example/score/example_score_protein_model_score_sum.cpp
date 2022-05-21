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
#include "score/bcl_score_protein_model_score_sum.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "math/bcl_math_template_instantiations.h"
#include "score/bcl_score_aa_neighborhood_exposure.h"
#include "score/bcl_score_aa_pair_clash.h"
#include "score/bcl_score_aa_pair_distance.h"
#include "score/bcl_score_aa_sequence_pair.h"
#include "score/bcl_score_loop.h"
#include "score/bcl_score_loop_closure.h"
#include "score/bcl_score_phi_psi.h"
#include "score/bcl_score_protein_model.h"
#include "score/bcl_score_protein_model_aa_neighborhood.h"
#include "score/bcl_score_protein_model_sse.h"
#include "score/bcl_score_protein_model_sse_neighbors.h"
#include "score/bcl_score_protein_model_sse_pairs.h"
#include "score/bcl_score_radius_of_gyration.h"
#include "score/bcl_score_sse_pair_packing.h"
#include "score/bcl_score_sse_pairs_fragments.h"
#include "score/bcl_score_strand_pairing.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_score_sum.cpp
  //!
  //! @author karakam
  //! @date 11/28/10
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinModelScoreSum :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinModelScoreSum *Clone() const
    { return new ExampleScoreProteinModelScoreSum( *this);}

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
      // construct min sse sizes map
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 4;
      min_sse_sizes[ biol::GetSSTypes().COIL] = 999;

      // get 1UBI model
      const std::string filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel native_model
      (
        Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // get the vector of SSEs
      util::SiPtrVector< const assemble::SSE> sse_vector( native_model.GetSSEs());

      // change min_sse_sizes to exclude everything
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 999;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 999;

      // initialize the scores
      // aa distance score
      util::ShPtr< score::ProteinModel> sp_aadist
      (
        new score::ProteinModelSSEPairs
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
          (
            new score::AASequencePair( score::AAPairDistance(), false)
          ),
          false // no normalization
        )
      );

      // aa distance clash
      util::ShPtr< score::ProteinModel> sp_aaclash
      (
        new score::ProteinModelSSEPairs
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
          (
            new score::AASequencePair( score::AAPairClash(), false)
          ),
          false // no normalization
        )
      );

      // loop scores
      util::ShPtr< score::ProteinModel> sp_loop
      (
        new score::ProteinModelSSEPairs
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >( new score::Loop()),
          true // normalization
        )
      );

      // loop closure score
      util::ShPtr< score::ProteinModel> sp_loop_closure( new score::ProteinModelSSENeighbors( util::CloneToShPtr( score::LoopClosure( 1, 1.0, 1.0)), false));

      // phi psi
      util::ShPtr< score::ProteinModel> sp_phi_psi
      (
        new score::ProteinModelSSE
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >
          (
            new score::PhiPsi( score::PhiPsi::GetDefaultScheme())
          ),
          false
        )
      );

    ///////////////
    // operators //
    ///////////////

      // sse packing
      util::ShPtr< score::ProteinModel> sp_ssepack
      (
        new score::ProteinModelSSEPairs
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
          (
            new score::SSEPairsFragments
            (
              assemble::GetSSEGeometryPackingListPickers().e_BestInteractionWeight,
              score::SSEPairPacking( "ssepack_fr", "sse_fragment_angle_distance.histograms2D"),
              false
            )
          ),
          false
        )
      );

      // strand pairing
      util::ShPtr< score::ProteinModel> sp_strand
      (
        new score::ProteinModelSSEPairs
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
          (
            new score::SSEPairsFragments
            (
              assemble::GetSSEGeometryPackingListPickers().e_BestInteractionWeight,
              score::StrandPairing( "strand_fr", "strand_fragment_angle_distance.histograms2D"), false
            )
          ),
          false
        )
      );

      // radius of gyration
      util::ShPtr< score::ProteinModel> sp_rgyr( new score::RadiusOfGyration());

      // aa exposure scores
      util::ShPtr< score::ProteinModel> sp_aaneigh
      (
        new score::ProteinModelAANeighborhood
        (
          util::CloneToShPtr( score::AANeighborhoodExposure( assemble::AANeighborCount()))
        )
      );

      // construct score weight set
      storage::Map< util::ShPtr< score::ProteinModel>, double> score_weight_set;
      score_weight_set[ sp_aaclash] = 1000;
      score_weight_set[ sp_aadist] = 0.1;
      score_weight_set[ sp_aaneigh] = 11.2;
      score_weight_set[ sp_loop] = 20;
      score_weight_set[ sp_loop_closure] = 1000;
      score_weight_set[ sp_phi_psi] = 0.1;
      score_weight_set[ sp_rgyr] = 4.0;
      score_weight_set[ sp_ssepack] = 5.8;
      score_weight_set[ sp_strand] =  2.5;

      // initialize function schemes
      storage::Vector< std::string> function_schemes_a;
      // iterate over the map
      for
      (
        storage::Map< util::ShPtr< score::ProteinModel>, double>::const_iterator
        itr( score_weight_set.Begin()), itr_end( score_weight_set.End()); itr != itr_end; ++itr
      )
      {
        function_schemes_a.PushBack( itr->first->GetScheme());
      }
      function_schemes_a.PushBack( "sum");

      // construct score map
      storage::Map< std::string, util::ShPtr< score::ProteinModel> > score_map;
      score_map[ sp_aadist->GetScheme()] = sp_aadist;
      score_map[ sp_aaneigh->GetScheme()] = sp_aaneigh;
      score_map[ sp_loop->GetScheme()] = sp_loop;
      score_map[ sp_rgyr->GetScheme()] = sp_rgyr;

      // construct weight set
      storage::Map< std::string, double> weight_set;
      weight_set[ sp_aadist->GetScheme()] = 0.1;
      weight_set[ sp_aaneigh->GetScheme()] = 11.2;
      weight_set[ sp_loop->GetScheme()] = 20;
      weight_set[ sp_rgyr->GetScheme()] = 4.0;

      // construct corresponding function schemes
      storage::Vector< std::string> function_schemes_b;
      function_schemes_b.PushBack( sp_aadist->GetScheme());
      function_schemes_b.PushBack( sp_aaneigh->GetScheme());
      function_schemes_b.PushBack( sp_loop->GetScheme());
      function_schemes_b.PushBack( sp_rgyr->GetScheme());
      function_schemes_b.PushBack( "sum");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test constructors");

      // test default constructor
      score::ProteinModelScoreSum score_default;

      // test constructor from a score weight map
      score::ProteinModelScoreSum score_a( score_weight_set);
      BCL_ExampleCheck( score_a.GetFunctionSchemes(), function_schemes_a);

      // test constructor from a score map and a weight set
      score::ProteinModelScoreSum score_b( score_map, weight_set);
      BCL_ExampleCheck( score_b.GetFunctionSchemes(), function_schemes_b);

      // test copy constructor
      score::ProteinModelScoreSum score_copy( score_a);
      BCL_ExampleCheck( score_copy.GetFunctionSchemes(), function_schemes_a);

      // test clone constructor
      util::ShPtr< score::ProteinModelScoreSum> sp_score_b( score_b.Clone());
      BCL_ExampleCheck( sp_score_b->GetFunctionSchemes(), function_schemes_b);

    /////////////////
    // data access //
    /////////////////

      // test GetStaticClassName
      BCL_ExampleCheck( GetStaticClassName< score::ProteinModelScoreSum>(), "bcl::score::ProteinModelScoreSum");

      // test CreateValueTableHorizontal()
      BCL_MessageStd( "testing CreateValueTableHorizontal()");

      // now create a new model where there are no SSEs
      assemble::ProteinModel model( native_model);

      // initialize expected score_sum
      storage::Table< double> table_horizontal( score_b.CreateValueTableHorizontal( model));
      table_horizontal.WriteFormatted( util::GetLogger());
      const double score_horizontal( table_horizontal[ "weighted_value"].GetData().LastElement());
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        score_horizontal, -227.368, 0.001, "score_horizontal: " + util::Format()( score_horizontal)
      );

      // test CreateValueTableVertical()
      BCL_MessageStd( "testing CreateValueTableVertical()");
      storage::Table< double> table_vertical( score_b.CreateValueTableVertical( model));
      table_vertical.WriteFormatted( util::GetLogger());
      const double score_vertical( table_vertical[ "sum"].GetData().LastElement());
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        score_horizontal, score_vertical, 0.001, "score_vertical: " + util::Format()( score_vertical)
      );

    ///////////////
    // operators //
    ///////////////

      // test operator()
      BCL_MessageStd( "testing operator()");

      // reset the model
      model.FilterByMinSSESizes( min_sse_sizes);

      // expected score of empty model
      const double expected_score_empty( 0.0);
      // evaluate the score
      double this_score( score_copy( model));
      // display score and check with expected score
      BCL_MessageStd( "Score for empty model: " + util::Format()( this_score));
      BCL_ExampleCheckWithinTolerance( expected_score_empty, this_score, 0.01);

      // vector of expected scores
      const storage::Vector< double> expected_scores
      (
        storage::Vector< double>::Create( 89.0623, 61.7865, 11.0729, -5.35768, -510.442)
      );

      // insert first SSE, evaluate and check the result
      model.Insert( util::ShPtr< assemble::SSE>( sse_vector( 0)->Clone()));
      this_score = score_copy( model);
      BCL_MessageStd( "Score for model with 1 SSE : " + util::Format()( this_score));
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        expected_scores( 0), this_score, 0.01, "result: " + util::Format()( this_score)
      );

      // insert second SSE, evaluate and check the result
      model.Insert( util::ShPtr< assemble::SSE>( sse_vector( 1)->Clone()));
      this_score = score_copy( model);
      BCL_MessageStd( "Score for model with 2 SSEs : " + util::Format()( this_score));
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        expected_scores( 1), this_score, 0.01, "result: " + util::Format()( this_score)
      );

      // insert third SSE, evaluate and check the result
      model.Insert( util::ShPtr< assemble::SSE>( sse_vector( 2)->Clone()));
      this_score = score_copy( model);
      BCL_MessageStd( "Score for model with 3 SSEs : " + util::Format()( this_score));
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        expected_scores( 2), this_score, 0.01, "result: " + util::Format()( this_score)
      );

      // insert fourth SSE, evaluate and check the result
      model.Insert( util::ShPtr< assemble::SSE>( sse_vector( 3)->Clone()));
      this_score = score_copy( model);
      BCL_MessageStd( "Score for model with 4 SSEs : " + util::Format()( this_score));
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        expected_scores( 3), this_score, 0.01, "result: " + util::Format()( this_score)
      );

      // insert fifth SSE, evaluate and check the result
      model.Insert( util::ShPtr< assemble::SSE>( sse_vector( 4)->Clone()));
      this_score = score_copy( model);
      BCL_MessageStd( "Score for model with 5 SSEs : " + util::Format()( this_score));
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        expected_scores( 4), this_score, 0.01, "result: " + util::Format()( this_score)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test WriteDetailedSchemeAndValues()
      BCL_MessageStd( "testing WriteDetailedSchemeAndValues()");
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        score_b.WriteDetailedSchemeAndValues( model, util::GetLogger());
      }

      // test ReadWrite
      BCL_MessageStd( "testing ReadWrite");

      // construct a score with just radius of gyration and remember the score
      storage::Map< util::ShPtr< score::ProteinModel>, double> weight_set_rgyr;
      weight_set_rgyr[ sp_rgyr] = 1.0;
      score::ProteinModelScoreSum score_rgyr( weight_set_rgyr);
      const double expected_score_rgyr( score_rgyr( model));
      BCL_MessageStd( "score rgyr: " + util::Format()( expected_score_rgyr));

      // write the score and then read it back in
      WriteBCLObject( score_rgyr);
      score::ProteinModelScoreSum score_read;
      ReadBCLObject( score_read);

      // calculate and compare the values
      const double this_score_rgyr( score_read( model));
      BCL_MessageStd( "score read rgyr: " + util::Format()( this_score_rgyr));
      BCL_ExampleCheckWithinTolerance( expected_score_rgyr, this_score_rgyr, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinModelScoreSum

  const ExampleClass::EnumType ExampleScoreProteinModelScoreSum::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelScoreSum())
  );

} // namespace bcl
