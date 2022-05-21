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
#include "score/bcl_score_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "score/bcl_score_protein_model_sse.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_phi_psi.cpp
  //!
  //! @author rouvelgh
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by karakam on Dec 13, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScorePhiPsi :
    public ExampleInterface
  {
  public:

    ExampleScorePhiPsi *Clone() const
    {
      return new ExampleScorePhiPsi( *this);
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

      //instantiate pdb
      const std::string filename( AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
      //instantiate protein model of chains
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 4;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;
      // make a model with strands and helices only of min size 4 residues each and backbone atoms only
      const assemble::ProteinModel
        model( Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size));

      const biol::Membrane sp_membrane;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from scheme and histogram filename
      const score::PhiPsi score_phi_psi( "sse_bend", score::PhiPsi::GetDefaultHistogramFilename());

      BCL_Example_Check
      (
        score_phi_psi.GetScheme() == "sse_bend",
        "scheme was not initialize to \"sse_bend\" but is: " + score_phi_psi.GetScheme()
      );

      BCL_ExampleCheck( score_phi_psi.GetScheme(), "sse_bend");

    /////////////////
    // data access //
    /////////////////

      // output scheme
      BCL_MessageStd( "scheme for constructed score::PhiPsi: " + score_phi_psi.GetScheme());

      // class identifier and static class name
      BCL_ExampleCheck( GetStaticClassName< score::PhiPsi>(), "bcl::score::PhiPsi");

      // histogram filename
      BCL_MessageStd( "histogram file used for generating energies: " + score_phi_psi.GetHistogramFilename());

      // access to the energy function map
      BCL_ExampleIndirectCheck
      (
        score_phi_psi.GetEnergyFunctions().GetSize(), 3,
        "there should be at least 3 energy functions for 3 SSTypes in the map!"
      );

      // access to the energy function map
      BCL_ExampleIndirectCheck
      (
        score_phi_psi.GetAATypeEnergyFunctions().GetSize() >= 20, true,
        "there should be at least 20 energy functions for 20 natural AAs in the map!"
      );

      // for debug write gnuplots for spline to files
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // iterate over all energy functions
        for
        (
          storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> >::const_iterator
            map_itr( score_phi_psi.GetEnergyFunctions().Begin()),
            map_itr_end( score_phi_psi.GetEnergyFunctions().End());
          map_itr != map_itr_end;
          ++map_itr
        )
        {
          // iterate over all energy functions
          for
          (
            storage::Map< biol::AAType, math::BicubicSpline>::const_iterator
              itr( map_itr->second.Begin()), itr_end( map_itr->second.End());
            itr != itr_end;
            ++itr
          )
          {
            // write function to gnuplot file
            io::OFStream write;
            BCL_ExampleMustOpenOutputFile( write, "phi_psi_" + itr->first.GetName() + "_" + map_itr->first.GetName() + ".gnuplot");
            math::GnuplotHeatmap heatmap;
            heatmap.SetFromBicubicSpline( itr->second, false, false);
            heatmap.SetTitleAndLabel( "phi psi potential " + map_itr->first.GetName() + " " + itr->first.GetName(), "phi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "psi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "-log(p/b)");
            heatmap.SetFont( "arialbd", 16);
            heatmap.SetPixelAndRatio( 1080, 800, -1);
            const linal::Vector< double> binning
            (
              linal::FillVector< double>
              (
                itr->second.GetValues().GetNumberCols() + 1,
                math::Angle::Degree( itr->second.GetStartY() - 0.5 * itr->second.GetDeltaY()),
                math::Angle::Degree( itr->second.GetDeltaY())
              )
            );
            heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
            heatmap.SetTicsY( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
            heatmap.SetFilename( "phi_psi_" + itr->first.GetName() + "_" + map_itr->first.GetName());
            heatmap.WriteScript( write);
            io::File::CloseClearFStream( write);
          }
        }
        // iterate over all energy functions
        for
        (
          storage::Map< biol::SSType, storage::Map< biol::AAType, math::BicubicSpline> >::const_iterator
            map_itr( score_phi_psi.GetEnergyFunctionsMembrane().Begin()),
            map_itr_end( score_phi_psi.GetEnergyFunctionsMembrane().End());
          map_itr != map_itr_end;
          ++map_itr
        )
        {
          // iterate over all energy functions
          for
          (
            storage::Map< biol::AAType, math::BicubicSpline>::const_iterator
              itr( map_itr->second.Begin()), itr_end( map_itr->second.End());
            itr != itr_end;
            ++itr
          )
          {
            // write function to gnuplot file
            io::OFStream write;
            BCL_ExampleMustOpenOutputFile( write, "phi_psi_" + itr->first.GetName() + "_" + map_itr->first.GetName() + ".membrane..gnuplot");
            math::GnuplotHeatmap heatmap;
            heatmap.SetFromBicubicSpline( itr->second, false, false);
            heatmap.SetTitleAndLabel( "phi psi potential " + map_itr->first.GetName() + " " + itr->first.GetName(), "phi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "psi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "-log(p/b)");
            heatmap.SetFont( "arialbd", 16);
            heatmap.SetPixelAndRatio( 1080, 800, -1);
            const linal::Vector< double> binning
            (
              linal::FillVector< double>
              (
                itr->second.GetValues().GetNumberCols() + 1,
                math::Angle::Degree( itr->second.GetStartY() - 0.5 * itr->second.GetDeltaY()),
                math::Angle::Degree( itr->second.GetDeltaY())
              )
            );
            heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
            heatmap.SetTicsY( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
            heatmap.SetFilename( "phi_psi_" + itr->first.GetName() + "_" + map_itr->first.GetName() + ".membrane");
            heatmap.WriteScript( write);
            io::File::CloseClearFStream( write);
          }
        }
        // iterate over all energy functions
        for
        (
          storage::Map< biol::AAType, math::BicubicSpline>::const_iterator
            itr( score_phi_psi.GetAATypeEnergyFunctions().Begin()), itr_end( score_phi_psi.GetAATypeEnergyFunctions().End());
          itr != itr_end;
          ++itr
        )
        {
          // write function to gnuplot file
          io::OFStream write;
          BCL_ExampleMustOpenOutputFile( write, "phi_psi_" + itr->first.GetName() + ".gnuplot");
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromBicubicSpline( itr->second, false, false);
          heatmap.SetTitleAndLabel( "phi psi potential " + itr->first.GetName(), "phi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "psi [" + math::Angle::s_DegreeSymbolGnuplot + "]", "-log(p/b)");
          heatmap.SetPixelAndRatio( 1080, 800, -1);
          heatmap.SetFont( "arialbd", 16);
          const linal::Vector< double> binning
          (
            linal::FillVector< double>
            (
              itr->second.GetValues().GetNumberCols() + 1,
              math::Angle::Degree( itr->second.GetStartY() - 0.5 * itr->second.GetDeltaY()),
              math::Angle::Degree( itr->second.GetDeltaY())
            )
          );
          heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
          heatmap.SetTicsY( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 0).W( 3).Fill( ' ').R()), false, 1);
          heatmap.SetFilename( "phi_psi_" + itr->first.GetName());
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }

      }

    ///////////////
    // operators //
    ///////////////

      // iterate over all sses in the first chain and score them
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( model.GetChains()( 0)->GetData().Begin()), sse_itr_end( model.GetChains()( 0)->GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        const storage::Pair< double, size_t> score( score_phi_psi( **sse_itr, sp_membrane));
        BCL_MessageStd
        (
          ( *sse_itr)->GetIdentification() + ":\t" +
          util::Format()( score.First())
        );
      }

      // test clone
      util::ShPtr< score::PhiPsi> sp_score( score_phi_psi.Clone());

      // score the entire protein model and report its score
      const double score_of_model( score::ProteinModelSSE( sp_score, true)( model));
      const double expected_model_score( -5.23718);
      BCL_MessageStd( "phi_psi score of 1IE9.pdb " + util::Format()( score_of_model));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( expected_model_score, score_of_model), true,
        "expected phi_psi score for 1A0E.pdb was not met: " +
        util::Format()( score_of_model) + " != " + util::Format()( expected_model_score)
      );

      // test with coil
      {
        //instantiate protein model of chains
        storage::Map< biol::SSType, size_t> ssetype_min_size_coil;
        ssetype_min_size_coil[ biol::GetSSTypes().HELIX] = 4;
        ssetype_min_size_coil[ biol::GetSSTypes().STRAND] = 4;
        ssetype_min_size_coil[ biol::GetSSTypes().COIL] = 0;
        // make a model with strands and helices only of min size 4 residues each and backbone atoms only
        const assemble::ProteinModel model_coil
        (
          Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size_coil)
        );
        sp_score->SetSSTypes( storage::Set< biol::SSType>( biol::GetSSTypes().COIL));
        // score the entire protein model and report its score
        const double score
        (
          score::ProteinModelSSE
          (
            util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >( sp_score), true
          ).operator()( model_coil)
        );
        BCL_MessageStd( "phi_psi score of 1IE9.pdb " + util::Format()( score));
        BCL_ExampleIndirectCheckWithinTolerance( score, -4.35425, 0.001, "phi_psi score for 1IE9.pdb");
      }

      // test with helix
      {
        //instantiate protein model of chains
        storage::Map< biol::SSType, size_t> ssetype_min_size_coil;
        ssetype_min_size_coil[ biol::GetSSTypes().HELIX] = 4;
        ssetype_min_size_coil[ biol::GetSSTypes().STRAND] = 4;
        ssetype_min_size_coil[ biol::GetSSTypes().COIL] = 0;
        // make a model with strands and helices only of min size 4 residues each and backbone atoms only
        const assemble::ProteinModel model_coil
        (
          Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size_coil)
        );
        sp_score->SetSSTypes( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX));
        // score the entire protein model and report its score
        const double score
        (
          score::ProteinModelSSE
          (
            util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >( sp_score), true
          ).operator()( model_coil)
        );
        const double expected_score( -3.828568);
        BCL_MessageStd( "phi_psi score of 1IE9.pdb " + util::Format()( score));
        BCL_ExampleIndirectCheckWithinTolerance( score, -3.74025, 0.001, "phi_psi score for 1IE9.pdb in membrane");
      }

    //////////////////////
    // input and output //
    //////////////////////

      // write score_phi_psi to file
      WriteBCLObject( score_phi_psi);
      // read score::PhiPsi
      score::PhiPsi score_read;
      ReadBCLObject( score_read);

      BCL_ExampleCheck( score_phi_psi.GetScheme(), score_read.GetScheme());
      BCL_ExampleCheck( score_phi_psi.GetHistogramFilename(), score_read.GetHistogramFilename());

      // create ShPtr
      util::ShPtr< score::PhiPsi> sp_score_read( score_read.Clone());

      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( expected_model_score, score::ProteinModelSSE( sp_score_read, true)( model)),
        true,
        "score::PhiPsi was not read properly!"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScorePhiPsi

  const ExampleClass::EnumType ExampleScorePhiPsi::s_Instance
  (
    GetExamples().AddEnum( ExampleScorePhiPsi())
  );

} // namespace bcl
