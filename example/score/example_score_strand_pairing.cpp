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
#include "score/bcl_score_strand_pairing.h"

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
  //! @example example_score_strand_pairing.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreStrandPairing :
    public ExampleInterface
  {
  public:

    ExampleScoreStrandPairing *Clone() const
    {
      return new ExampleScoreStrandPairing( *this);
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
      util::SiPtr< const assemble::SSE> strand_a( assemble::LocatorSSE( 'A', 12, 18).Locate( model));
      util::SiPtr< const assemble::SSE> strand_b( assemble::LocatorSSE( 'A', 23, 31).Locate( model));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      BCL_MessageStd( "test default constructor");
      score::StrandPairing strand_pairing_def;

      BCL_MessageStd( "test constructor from nr_excluded and width and scheme");
      score::StrandPairing strand_pairing( "strand_fr", "strand_fragment_angle_distance.histograms2D");

      BCL_MessageStd( "test copy constructor");
      score::StrandPairing strand_pairing_copy( strand_pairing);

      BCL_MessageStd( "test clone");
      util::ShPtr< score::StrandPairing> sp_strand_pairing( strand_pairing.Clone());
      BCL_Example_Check
      (
        sp_strand_pairing->GetScheme() == strand_pairing.GetScheme(),
        "Cloned object should be \n" + util::Format()( strand_pairing) + "\nnot\n" + util::Format()( *sp_strand_pairing)
      );

    /////////////////
    // data access //
    /////////////////

      // test the GetClassIdentifier
      BCL_ExampleCheck( strand_pairing.GetClassIdentifier(), GetStaticClassName< score::StrandPairing>());

      // test GetMinimalInterfaceLength
      BCL_ExampleCheck( strand_pairing.GetMinimalInterfaceLength(), 4.0);

      // test GetScheme
      BCL_ExampleCheck( strand_pairing.GetScheme(), "strand_fr");

    ////////////////
    // operations //
    ////////////////

      // test are valid sses
      BCL_MessageStd( "test AreValidSSEs( arg)");
      BCL_Example_Check
      (
        !strand_pairing.AreValidSSEs( *helix_a, *strand_a) && strand_pairing.AreValidSSEs( *strand_a, *strand_b),
        "AreValidSSEs() should return true for strand pair and false for helix strand pair"
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
            new score::SSEPairsFragments
            (
              assemble::GetSSEGeometryPackingListPickers().e_BestInteractionWeight, strand_pairing, false
            )
          ),
          false
        );

      util::ShPtrVector< assemble::SSE> strandpair;
      strandpair.PushBack( util::ShPtr< assemble::SSE>( strand_a->Clone()));
      strandpair.PushBack( util::ShPtr< assemble::SSE>( strand_b->Clone()));
      strandpair( 0)->SetToIdealConformationAtOrigin();
      strandpair( 1)->SetToIdealConformationAtOrigin();

      // translate
      const size_t max_translation( 9);
      const size_t max_rot( 7);
      const size_t max_roty( 5);
      const double shift( 5.5 / max_translation);
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
      format.Fill( '0').W( 4).FFP( 3);

      // statistics for sse_packing scores
      math::RunningAverageSD< double> mean_sd_strand;
      math::RunningMinMax< double> min_max_strand;

      util::ShPtr< assemble::Chain> sp_chain_strand
      (
        new assemble::Chain( model.GetChains()( 0)->GetSequence(), strandpair)
      );
      const assemble::ProteinModel tmp_model_strand( sp_chain_strand);

      // translation
      for( size_t trans( 0); trans < max_translation; ++trans)
      {
        trans_external_y.Move( *strandpair( 1));
        for( size_t rot_e( 0); rot_e < max_rot; ++rot_e)
        {
          rot_external_z.Move( *strandpair( 1));
          for( size_t rot_iz( 0); rot_iz < max_rot; ++rot_iz)
          {
            rot_internal_z.Move( *strandpair( 1));
            for( size_t rot_iy( 0); rot_iy < max_roty; ++rot_iy)
            {
              rot_internal_y.Move( *strandpair( 1));

              const double score_strand_pairing( score_model( tmp_model_strand));

              mean_sd_strand +=  score_strand_pairing;
              min_max_strand += score_strand_pairing;

              BCL_MessageVrb
              (
                "translation: " + format( trans * shift) +
                " rot ext: " + format( rot_e * rotation) +
                " rot intz: " + format( rot_iz * rotation) +
                " rot inty: " + format( rot_iy * rotation) +
                " strand_pairing score strand: " + format( score_strand_pairing)
              );
//              io::OFStream write;
//              BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( sse_sse_packing, "ss" + util::Format().Fill( '0').FFP( 2,0)( trans) + util::Format().Fill( '0').FFP( 2,0)( rot_e) + util::Format().Fill( '0').FFP( 2,0)( rot_iz) + util::Format().Fill( '0').FFP( 2,0)( rot_iy) + format( sse_packing_strand) + ".pdb"));
//              pdb::Factory( biol::GetAAClasses().e_AABackBone).WriteModelToPDB( tmp_model_strand, write);
//              io::File::CloseClearFStream( write);
            }
          }
        }
      }

      // output mean, sd and min,max for each of the score
      BCL_MessageStd( "mean and sd for strand       sse_packing: " + format( mean_sd_strand.GetAverage())      + " " + format( mean_sd_strand.GetStandardDeviation()));
      BCL_MessageStd( "min and max for strand       sse_packing: " + format( min_max_strand.GetMin())       + " " + format( min_max_strand.GetMax()));

      const double expected_mean_strand(  8.86472);
      const double expected_sd_strand(    8.48446);

      const double expected_min_strand( -17.4934);
      const double expected_max_strand(  47.5694);

      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_mean_strand, mean_sd_strand.GetAverage(), 0.02) && math::EqualWithinTolerance( expected_sd_strand, mean_sd_strand.GetStandardDeviation(), 0.02),
        "expected mean and sd for strand sse_packing score " + format( expected_mean_strand) + " " + format( expected_sd_strand) + " is different from actual: " + util::Format()( mean_sd_strand.GetAverage()) + " " + util::Format()( mean_sd_strand.GetStandardDeviation())
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_min_strand, min_max_strand.GetMin()) && math::EqualWithinTolerance( expected_max_strand, min_max_strand.GetMax(), 0.12),
        "expected min and max for strand sse_packing score " + format( expected_min_strand) + " " + format( expected_max_strand) + " is different from actual: " + util::Format()( min_max_strand.GetMin()) + " " + util::Format()( min_max_strand.GetMax())
      );

    //////////////////////
    // input and output //
    //////////////////////

      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // write function to gnuplot file
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, contact::GetTypes().STRAND_STRAND.GetName() + ".gnuplot");
        math::GnuplotHeatmap heatmap;
        heatmap.SetFromBicubicSpline( strand_pairing.GetEnergyFunction(), false, false);
        heatmap.SetTitleAndLabel( "strand pairing potential", "angle [" + math::Angle::s_DegreeSymbolGnuplot + "]", "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "-log(p/b)");
        heatmap.SetFont( "arialbd", 20);
        const linal::Vector< double> x_binning
        (
          linal::FillVector< double>
          (
            strand_pairing.GetEnergyFunction().GetValues().GetNumberCols() + 1,
            math::Angle::Degree( strand_pairing.GetEnergyFunction().GetStartY() - 0.5 * strand_pairing.GetEnergyFunction().GetDeltaY()),
            math::Angle::Degree( strand_pairing.GetEnergyFunction().GetDeltaY())
          )
        );
        const linal::Vector< double> y_binning
        (
          linal::FillVector< double>
          (
            strand_pairing.GetEnergyFunction().GetValues().GetNumberRows() + 1,
            strand_pairing.GetEnergyFunction().GetStartX() - 0.5 * strand_pairing.GetEnergyFunction().GetDeltaX(),
            strand_pairing.GetEnergyFunction().GetDeltaX()
          )
        );
        heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( x_binning, 2, util::Format().FFP( 0).W( 4).Fill( ' ').R()), false, 2);
        heatmap.SetTicsY( math::GnuplotHeatmap::TicsFromBinning( y_binning, 2, util::Format().FFP( 1).W( 3).Fill( ' ').R()), false, 2);
        heatmap.SetFilename( "sse_packing_" + contact::GetTypes().STRAND_STRAND.GetName());
        heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);
      }

      // test read write
      BCL_MessageStd( "testing read write")
      WriteBCLObject( strand_pairing);
      score::StrandPairing strand_pairing_read;
      ReadBCLObject( strand_pairing_read);
      BCL_Example_Check
      (
        strand_pairing_read.GetScheme()                 == strand_pairing.GetScheme() &&
        strand_pairing_read.GetMinimalInterfaceLength() == strand_pairing.GetMinimalInterfaceLength() &&
        strand_pairing_read.GetHistogramFilename()      == strand_pairing.GetHistogramFilename(),
        "The read sse packing score is different\n" +
        strand_pairing_read.GetScheme()                                  + " vs " + strand_pairing.GetScheme() + "\n" +
        util::Format()( strand_pairing_read.GetMinimalInterfaceLength()) + " vs " + util::Format()( strand_pairing.GetMinimalInterfaceLength()) + "\n" +
        strand_pairing_read.GetHistogramFilename()                       + " vs " + strand_pairing.GetHistogramFilename()
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreStrandPairing

  const ExampleClass::EnumType ExampleScoreStrandPairing::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreStrandPairing())
  );

} // namespace bcl
