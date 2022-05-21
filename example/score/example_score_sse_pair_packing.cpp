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
#include "score/bcl_score_sse_pair_packing.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "coord/bcl_coord_move_rotate_defined.h"
#include "coord/bcl_coord_move_translate_defined.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "score/bcl_score_protein_model_sse_pairs.h"
#include "score/bcl_score_sse_pairs_fragments.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_sse_pair_packing.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSSEPairPacking :
    public ExampleInterface
  {
  public:

    ExampleScoreSSEPairPacking *Clone() const
    { return new ExampleScoreSSEPairPacking( *this);}

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

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // get the protein model
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes)
      );

      // locate helices
      util::SiPtr< const assemble::SSE> helix_a( assemble::LocatorSSE( 'A', 245, 252).Locate( model));
      util::SiPtr< const assemble::SSE> helix_b( assemble::LocatorSSE( 'A', 267, 277).Locate( model));
      util::SiPtr< const assemble::SSE> strand_a( assemble::LocatorSSE( 'A', 12, 18).Locate( model));
      util::SiPtr< const assemble::SSE> strand_b( assemble::LocatorSSE( 'A', 23, 31).Locate( model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test default constructor");
      score::SSEPairPacking sse_packing_def;

      BCL_MessageStd( "test constructor from nr_excluded and width and scheme");
      score::SSEPairPacking sse_packing( "ssepack_fr", "sse_fragment_angle_distance.histograms2D");

      BCL_MessageStd( "test copy constructor");
      score::SSEPairPacking sse_packing_copy( sse_packing);

      BCL_MessageStd( "test clone");
      util::ShPtr< score::SSEPairPacking> sp_sse_packing( sse_packing.Clone());
      BCL_Example_Check
      (
        sp_sse_packing->GetScheme() == sse_packing.GetScheme(),
        "Cloned object should be \n" + util::Format()( sse_packing) + "\nnot\n" + util::Format()( *sp_sse_packing)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( sse_packing.GetClassIdentifier(), GetStaticClassName< score::SSEPairPacking>());

      // test GetMinimalInterfaceLength
      BCL_ExampleCheck( sse_packing.GetMinimalInterfaceLength(), 4.0);

      // test GetScheme
      BCL_ExampleCheck( sse_packing.GetScheme(), "ssepack_fr");

      // test access to distance range map
      BCL_MessageStd( "test GetDistanceRangeMap()");
      BCL_Example_Check
      (
        sse_packing.GetDistanceRangeMap().GetSize() == 7,
        "GetDistanceRangeMap() should have size 7 not " + util::Format()( sse_packing.GetDistanceRangeMap().GetSize())
      );

      // test access to emergy functions
      BCL_MessageStd( "test GetEnergyFunctions()");
      BCL_Example_Check
      (
        sse_packing.GetEnergyFunctions().GetSize() == 7,
        "GetEnergyFunctions() should have size 7 not " + util::Format()( sse_packing.GetEnergyFunctions().GetSize())
      );

    ////////////////
    // operations //
    ////////////////

      // test are valid sses
      BCL_MessageStd( "test AreValidSSEs( arg)");
      BCL_Example_Check
      (
        sse_packing.AreValidSSEs( *helix_a, *strand_a) && sse_packing.AreValidSSEs( *strand_a, *strand_b),
        "AreValidSSEs() should return true for all sses pairs"
      );

    ///////////////
    // operators //
    ///////////////

      score::ProteinModelScoreSum score_model;
      score_model +=
        score::ProteinModelSSEPairs
        (
          util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, assemble::SSE, double> >
          (
            new score::SSEPairsFragments( assemble::GetSSEGeometryPackingListPickers().e_BestInteractionWeight, sse_packing, false)
          ),
          false
        );

      util::ShPtrVector< assemble::SSE> helixpair;
      helixpair.PushBack( util::ShPtr< assemble::SSE>( helix_a->Clone()));
      helixpair.PushBack( util::ShPtr< assemble::SSE>( helix_b->Clone()));
      helixpair( 0)->SetToIdealConformationAtOrigin();
      helixpair( 1)->SetToIdealConformationAtOrigin();

      util::ShPtrVector< assemble::SSE> helixstrandpair;
      helixstrandpair.PushBack( util::ShPtr< assemble::SSE>( helix_a->Clone()));
      helixstrandpair.PushBack( util::ShPtr< assemble::SSE>( strand_a->Clone()));
      helixstrandpair( 0)->SetToIdealConformationAtOrigin();
      helixstrandpair( 1)->SetToIdealConformationAtOrigin();

      util::ShPtrVector< assemble::SSE> strandpair;
      strandpair.PushBack( util::ShPtr< assemble::SSE>( strand_a->Clone()));
      strandpair.PushBack( util::ShPtr< assemble::SSE>( strand_b->Clone()));
      strandpair( 0)->SetToIdealConformationAtOrigin();
      strandpair( 1)->SetToIdealConformationAtOrigin();

      // translate
      const size_t max_translation( 5);
      const size_t max_rot( 7);
      const size_t max_roty( 5);
      const double shift( 7.0 / max_translation);
      const double rotation( ( 2 * math::g_Pi) / max_rot);
      const double rotationy( ( 2 * math::g_Pi) / max_roty);

      // internal rotation
      coord::MoveRotateDefined rot_internal_z( rotation, coord::GetAxes().e_Z, true);

      // external rotation
      coord::MoveRotateDefined rot_external_z( rotation, coord::GetAxes().e_Z, false);

      // external rotation
      coord::MoveRotateDefined rot_internal_y( rotationy, coord::GetAxes().e_Y, false);

      // y translation
      coord::MoveTranslateDefined trans_external_y( linal::Vector3D( 0.0, shift, 0.0), false);
      util::Format format;
      format.Fill( '0').W( 4).FFP( 4);

      // statistics for sse_packing scores
      math::RunningAverageSD< double> mean_sd_helix;
      math::RunningAverageSD< double> mean_sd_helixstrand;
      math::RunningAverageSD< double> mean_sd_strand;
      math::RunningMinMax< double> min_max_helix;
      math::RunningMinMax< double> min_max_helixstrand;
      math::RunningMinMax< double> min_max_strand;

      util::ShPtr< assemble::Chain> sp_chain_helix
      (
        new assemble::Chain( model.GetChains()( 0)->GetSequence(), helixpair)
      );
      util::ShPtr< assemble::Chain> sp_chain_helixstrand
      (
        new assemble::Chain( model.GetChains()( 0)->GetSequence(), helixstrandpair)
      );
      util::ShPtr< assemble::Chain> sp_chain_strand
      (
        new assemble::Chain( model.GetChains()( 0)->GetSequence(), strandpair)
      );
      const assemble::ProteinModel tmp_model_helix( sp_chain_helix);
      const assemble::ProteinModel tmp_model_helixstrand( sp_chain_helixstrand);
      const assemble::ProteinModel tmp_model_strand( sp_chain_strand);

      // initial translation
      for( size_t i( 0); i < 7; ++i)
      {
        trans_external_y.Move( *helixpair( 1));
        trans_external_y.Move( *helixstrandpair( 1));
        trans_external_y.Move( *strandpair( 1));
      }

      // translation
      for( size_t trans( 0); trans < max_translation; ++trans)
      {
        trans_external_y.Move( *helixpair( 1));
        trans_external_y.Move( *helixstrandpair( 1));
        trans_external_y.Move( *strandpair( 1));
        for( size_t rot_e( 0); rot_e < max_rot; ++rot_e)
        {
          rot_external_z.Move( *helixpair( 1));
          rot_external_z.Move( *helixstrandpair( 1));
          rot_external_z.Move( *strandpair( 1));
          for( size_t rot_iz( 0); rot_iz < max_rot; ++rot_iz)
          {
            rot_internal_z.Move( *helixpair( 1));
            rot_internal_z.Move( *helixstrandpair( 1));
            rot_internal_z.Move( *strandpair( 1));
            for( size_t rot_iy( 0); rot_iy < max_roty; ++rot_iy)
            {
              rot_internal_y.Move( *helixpair( 1));
              rot_internal_y.Move( *helixstrandpair( 1));
              rot_internal_y.Move( *strandpair( 1));

              const double sse_packing_helix( score_model( tmp_model_helix));
              const double sse_packing_helixstrand( score_model( tmp_model_helixstrand));
              const double sse_packing_strand( score_model( tmp_model_strand));

              mean_sd_helix +=  sse_packing_helix;
              mean_sd_helixstrand += sse_packing_helixstrand;
              mean_sd_strand +=  sse_packing_strand;
              min_max_helix += sse_packing_helix;
              min_max_helixstrand += sse_packing_helixstrand;
              min_max_strand += sse_packing_strand;

              BCL_MessageVrb
              (
                "translation: " + format( trans * shift) +
                " rot ext: " + format( rot_e * rotation) +
                " rot intz: " + format( rot_iz * rotation) +
                " rot inty: " + format( rot_iy * rotation) +
                " sse_packing score helix: " + format( sse_packing_helix) +
                " sse_packing score helixstrand: " + format( sse_packing_helixstrand) +
                " sse_packing score strand: " + format( sse_packing_strand)
              );
//              if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
//              {
//                Proteins::WriteModelToPDB( tmp_model_helix, AddExampleOutputPathToFilename( sse_packing, "hh" + util::Format().Fill( '0').W( 2).FFP( 0)( trans) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_e) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_iz) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_iy) + format( sse_packing_helix)+ ".pdb"));
//                Proteins::WriteModelToPDB( tmp_model_helixstrand, AddExampleOutputPathToFilename( sse_packing, "hs" + util::Format().Fill( '0').W( 2).FFP( 0)( trans) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_e) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_iz) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_iy) + format( sse_packing_helixstrand)+ ".pdb"));
//                Proteins::WriteModelToPDB( tmp_model_strand, AddExampleOutputPathToFilename( sse_packing, "ss" + util::Format().Fill( '0').W( 2).FFP( 0)( trans) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_e) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_iz) + util::Format().Fill( '0').W( 2).FFP( 0)( rot_iy) + format( sse_packing_strand) + ".pdb"));
//              }
            }
          }
        }
      }

      // output mean, sd and min,max for each of the score
      BCL_MessageStd( "mean and sd for helix        sse_packing: " + format( mean_sd_helix.GetAverage())       + ' ' + format( mean_sd_helix.GetStandardDeviation()));
      BCL_MessageStd( "mean and sd for helix strand sse_packing: " + format( mean_sd_helixstrand.GetAverage()) + ' ' + format( mean_sd_helixstrand.GetStandardDeviation()));
      BCL_MessageStd( "mean and sd for strand       sse_packing: " + format( mean_sd_strand.GetAverage())      + ' ' + format( mean_sd_strand.GetStandardDeviation()));
      BCL_MessageStd( "min and max for helix        sse_packing: " + format( min_max_helix.GetMin())        + ' ' + format( min_max_helix.GetMax()));
      BCL_MessageStd( "min and max for helix strand sse_packing: " + format( min_max_helixstrand.GetMin())  + ' ' + format( min_max_helixstrand.GetMax()));
      BCL_MessageStd( "min and max for strand       sse_packing: " + format( min_max_strand.GetMin())       + ' ' + format( min_max_strand.GetMax()));

      const double expected_mean_helix(       -1.9758);
      const double expected_sd_helix(          1.9838);
      const double expected_mean_helixstrand( -2.2474);
      const double expected_sd_helixstrand(    3.2738);
      const double expected_mean_strand(      -0.3228);
      const double expected_sd_strand(         0.7447);

      const double expected_min_helix(        -7.4520);
      const double expected_max_helix(         0.8705);
      const double expected_min_helixstrand( -10.4553);
      const double expected_max_helixstrand(   4.0886);
      const double expected_min_strand(       -6.1820);
      const double expected_max_strand(        0.4263);

      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_mean_helix, mean_sd_helix.GetAverage()) && math::EqualWithinTolerance( expected_sd_helix, mean_sd_helix.GetStandardDeviation()),
        "expected mean and sd for helix sse_packing different from actual " + format( mean_sd_helix.GetAverage()) + " " + format( mean_sd_helix.GetStandardDeviation())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_mean_helixstrand, mean_sd_helixstrand.GetAverage()) && math::EqualWithinTolerance( expected_sd_helixstrand, mean_sd_helixstrand.GetStandardDeviation()),
        "expected mean and sd for helixstrand sse_packing different from actual " + format( mean_sd_helixstrand.GetAverage()) + " " + format( mean_sd_helixstrand.GetStandardDeviation())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_mean_strand, mean_sd_strand.GetAverage()) && math::EqualWithinTolerance( expected_sd_strand, mean_sd_strand.GetStandardDeviation()),
        "expected mean and sd for strand sse_packing different from actual " + format( mean_sd_strand.GetAverage()) + " " + format( mean_sd_strand.GetStandardDeviation())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_min_helix, min_max_helix.GetMin()) && math::EqualWithinTolerance( expected_max_helix, min_max_helix.GetMax()),
        "expected min and max for helix sse_packing different from actual " + format( min_max_helix.GetMin()) + " " + format( min_max_helix.GetMax())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_min_helixstrand, min_max_helixstrand.GetMin()) && math::EqualWithinTolerance( expected_max_helixstrand, min_max_helixstrand.GetMax()),
        "expected min and max for helixstrand sse_packing different from actual " + format( min_max_helixstrand.GetMin()) + " " + format( min_max_helixstrand.GetMax())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_min_strand, min_max_strand.GetMin()) && math::EqualWithinTolerance( expected_max_strand, min_max_strand.GetMax()),
        "expected min and max for strand sse_packing different from actual " + format( min_max_strand.GetMin()) + " " + format( min_max_strand.GetMax())
      );

      // print heatmap for sse packing potentials different contact types
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        for
        (
          storage::Map< contact::Type, util::ShPtr< math::BicubicSpline> >::const_iterator
            itr( sse_packing.GetEnergyFunctions().Begin()),
            itr_end( sse_packing.GetEnergyFunctions().End());
          itr != itr_end;
          ++itr
        )
        {
          // write function to gnuplot file
          io::OFStream write;
          BCL_ExampleMustOpenOutputFile( write, "sse_packing_" + itr->first.GetName() + ".gnuplot");
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromBicubicSpline( *itr->second, false, false);
          heatmap.SetTitleAndLabel
          (
            "secondary structure element packing potential: " + itr->first.GetName(),
            "angle [" + math::Angle::s_DegreeSymbolGnuplot + "]",
            "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]",
            "-log(p/b)"
          );
          const linal::Vector< double> x_binning
          (
            linal::FillVector< double>
            (
              ( itr->second->GetValues().GetNumberCols() + 1) / 8,
              math::Angle::Degree( itr->second->GetStartY() - 0.5 * itr->second->GetDeltaY()),
              math::Angle::Degree( itr->second->GetDeltaY()) * 8
            )
          );
          heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( x_binning, 8, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 8);
          heatmap.SetFont( "arialbd", 20);
          heatmap.SetFilename( "sse_packing_" + itr->first.GetName());
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test read write
      BCL_MessageStd( "testing read write")
      WriteBCLObject( sse_packing);
      score::SSEPairPacking sse_packing_read;
      ReadBCLObject( sse_packing_read);
      BCL_MessageStd( "read in ");
      BCL_ExampleIndirectCheck( sse_packing_read.GetScheme(), sse_packing.GetScheme(), "I/O");
      BCL_ExampleIndirectCheck( sse_packing_read.GetMinimalInterfaceLength(), sse_packing.GetMinimalInterfaceLength(), "I/O");
      BCL_ExampleIndirectCheck( sse_packing_read.GetHistogramFilename(), sse_packing.GetHistogramFilename(), "I/O");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreSSEPairPacking

  const ExampleClass::EnumType ExampleScoreSSEPairPacking::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSSEPairPacking())
  );

} // namespace bcl
