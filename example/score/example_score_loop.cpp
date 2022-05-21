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
#include "score/bcl_score_loop.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "score/bcl_score_protein_model_sse_pairs.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_loop.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreLoop :
    public ExampleInterface
  {
  public:

    ExampleScoreLoop *Clone() const
    {
      return new ExampleScoreLoop( *this);
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
      //instantiate pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      //instantiate protein models of chains
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      ssetype_min_size[ biol::GetSSTypes().COIL] = 999;
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor uses the default histogram files
      score::Loop loop_score;

      // clone
      util::ObjectInterface *ptr( loop_score.Clone());

      // build a score_model
      score::ProteinModelScoreSum score_model;

      //loop score between ends of sses using a 2d spline residues-distance
      score_model += score::ProteinModelSSEPairs( util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >( loop_score.Clone()), false);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( GetStaticClassName< score::Loop>(), ptr->GetClassIdentifier());

      // access the scheme
      BCL_ExampleCheck( loop_score.GetScheme(), "loop");

      // access the the energy functions
      BCL_ExampleCheck( loop_score.GetEnergyFunctions().GetSize(), loop_score.GetDefaultMaxLoopLength() + 1);

      // access the max loop length
      BCL_ExampleCheck( loop_score.GetMaxLoopLength(), loop_score.GetDefaultMaxLoopLength());

    ///////////////
    // operators //
    ///////////////

      // operator that score a pair of sses
      BCL_ExampleCheckWithinAbsTolerance( loop_score( *model.GetSSEs()( 0), *model.GetSSEs()( 1)), -2.42917, 1e-5);

      // usage in conjunction with ProteinModelScoreSum
      BCL_ExampleCheckWithinAbsTolerance( score_model( model), -1226.87, 0.01);

    ////////////////
    // operations //
    ////////////////

      // directly score a seq and euc distance
      const storage::Pair< size_t, double> seq_euc_dist( 2, 5.0);
      BCL_MessageStd( "loop length: " + util::Format()( seq_euc_dist.First()));
      BCL_ExampleCheckWithinTolerance( loop_score.Score( seq_euc_dist), -2.68063, 0.0001);

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( loop_score);
      // read from file
      score::Loop loop_score_read;
      ReadBCLObject( loop_score_read);

      BCL_ExampleCheckWithinTolerance
      (
        loop_score( *model.GetSSEs()( 0), *model.GetSSEs()( 1)),
        loop_score_read( *model.GetSSEs()( 0), *model.GetSSEs()( 1)),
        1e-6
      );

      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        storage::Vector< std::string> tics_x;
        tics_x.AllocateMemory( loop_score.GetEnergyFunctions().GetSize());
        for( size_t i( 0); i < loop_score.GetEnergyFunctions().GetSize(); ++i)
        {
          tics_x.PushBack( util::Format()( i));
        }
        // write function to gnuplot file
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, "loop.gnuplot");
        {
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromCubicSplines
          (
            util::ConvertToConstSiPtrVector< math::CubicSplineDamped>( loop_score.GetEnergyFunctions()),
            true,
            false
          );
          heatmap.SetTitleAndLabel( "loop", "number residues", "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "-log(p/b)");
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetPixelAndRatio( 1200, 800, -1.0);
          heatmap.SetTicsX( tics_x, true, 1);
          heatmap.SetMinMaxZ( -2.5, 2.5);
          heatmap.SetFilename( "loop");
          heatmap.WriteScript( write);
        }
        io::File::CloseClearFStream( write);
        BCL_ExampleMustOpenOutputFile( write, "loop_long.gnuplot");
        {
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromCubicSpline( loop_score.GetEnergyFunctionAboveMaxLoopLength(), false, false);
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetTitleAndLabel
          (
            "energy for loops > " + util::Format()( loop_score.GetMaxLoopLength()) + " residues",
            "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "", "-log(p/b)"
          );
          heatmap.SetFilename( "loop_long");
          heatmap.WriteScript( write);
        }
        io::File::CloseClearFStream( write);
      }

// Following code is needed for generating the thresholds to be used in loop closure score
//      {
//        const std::string histogram_filename
//        (
//          "/home/karakam/folding_benchmark/score/loop_closure/5_3/loop_closure.stats"
//        );
//        const size_t max_nr_loop( 29);
//
//        storage::Vector< double> max_all( score::Loop::CalculateMaximumObservedDistances( histogram_filename, max_nr_loop, 1.0));
//        storage::Vector< double> max_99( score::Loop::CalculateMaximumObservedDistances( histogram_filename, max_nr_loop, 0.99));
//        storage::Vector< double> max_95( score::Loop::CalculateMaximumObservedDistances( histogram_filename, max_nr_loop, 0.95));
//
//        std::cout << "counts\tall\t99\t95" << std::endl;
//        for( size_t i( 0); i <= max_nr_loop; ++i)
//        {
//          std::cout << i << "\t" << max_all( i) << "\t" << max_99( i) << "\t" << max_95( i) << "\t" << std::endl;
//        }
//        io::File::CloseClearFStream( write);
//      }

    //////////////////////
    // helper functions //
    //////////////////////

      delete ptr;

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreLoop

  const ExampleClass::EnumType ExampleScoreLoop::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreLoop())
  );

} // namespace bcl

