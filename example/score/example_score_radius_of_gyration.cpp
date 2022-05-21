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
#include "score/bcl_score_radius_of_gyration.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_radius_of_gyration.cpp
  //!
  //! @author rouvelgh
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRadiusOfGyration :
    public ExampleInterface
  {
  public:

    ExampleScoreRadiusOfGyration *Clone() const
    {
      return new ExampleScoreRadiusOfGyration( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // score for chain and models histogram
      score::RadiusOfGyration score_chains( false, false, "rgyr_chain", "radius_of_gyration_chains.histogram");
      score::RadiusOfGyration score_model( false, false, "rgyr_model", "radius_of_gyration_model.histogram", "radius_of_gyration_model.histogram");
      util::ShPtr< score::RadiusOfGyration> clone_constr( score_chains.Clone());

      // create protein model (with bcl sidechains) to test SquareRadiusOfGyration on
      biol::AASideChainFactory this_object( false, true);
      const std::string pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel this_model( Proteins::GetModel( pdb, biol::GetAAClasses().e_AAComplete));
      util::ShPtr< assemble::ProteinModel> model_ptr( this_object.ProteinModelWithSideChains( this_model));

      util::ShPtr< assemble::ProteinModel> membrane_protein_ptr( model_ptr->HardCopy());
      assemble::ProteinModel &membrane_protein( *membrane_protein_ptr);
      util::ShPtr< assemble::ProteinModelData> sp_data( membrane_protein.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, util::ShPtr< biol::Membrane>( new biol::Membrane()));
      membrane_protein.SetProteinModelData( sp_data);

      storage::Set< biol::AtomType> atom_set
      (
        storage::Set< biol::AtomType>::Create
        (
          biol::GetAtomTypes().CB, biol::GetAtomTypes().N, biol::GetAtomTypes().CA
        )
      );

      // construct splines
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // output the gnuplot heatmap files
        io::OFStream write;
        {
          const linal::Vector< double> binning( score_chains.GetEnergyFunctionSoluble().GetXValues());
          BCL_ExampleMustOpenOutputFile( write, "radius_of_gyration_chains.gnuplot");
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromCubicSpline( score_chains.GetEnergyFunctionSoluble(), false, false);
          heatmap.SetTitleAndLabel( "radius of gyration potential chains", "square radius of gyration / number of residues [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "^2]", "", "-log(p)");
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetPixelAndRatio( 1080, 800, -2.0);
          heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 8, util::Format().FFP( 1).W( 3).Fill( ' ').R()), false, 2);
          heatmap.SetRotationXTics( 90);
          heatmap.SetFilename( "radius_of_gyration_chains");
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }
        {
          const linal::Vector< double> binning( score_model.GetEnergyFunctionSoluble().GetXValues());
          BCL_ExampleMustOpenOutputFile( write, "radius_of_gyration_models.gnuplot");
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromCubicSpline( score_model.GetEnergyFunctionSoluble(), false, false);
          heatmap.SetTitleAndLabel( "radius of gyration potential models", "square radius of gyration / number of residues [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "^2]", "", "-log(p)");
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetPixelAndRatio( 1080, 800, -2.0);
          heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 8, util::Format().FFP( 1).W( 3).Fill( ' ').R()), false, 2);
          heatmap.SetRotationXTics( 90);
          heatmap.SetFilename( "radius_of_gyration_models");
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }
        {
          const linal::Vector< double> binning( score_chains.GetEnergyFunctionMembrane().GetXValues());
          BCL_ExampleMustOpenOutputFile( write, "radius_of_gyration_membrane.gnuplot");
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromCubicSpline( score_chains.GetEnergyFunctionMembrane(), false, false);
          heatmap.SetTitleAndLabel( "radius of gyration potential membrane", "square radius of gyration / number of residues [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "^2]", "", "-log(p)");
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetPixelAndRatio( 1080, 800, -2.0);
          heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 8, util::Format().FFP( 1).W( 3).Fill( ' ').R()), false, 2);
          heatmap.SetRotationXTics( 90);
          heatmap.SetFilename( "radius_of_gyration_membrane");
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }

        // load table of raw contact order statistics
        const std::string filename( AddExampleInputPathToFilename( e_Biology, "radius_of_gyration.table"));
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, filename);
        storage::Table< double> table;
        table.ReadFormatted( read);
        io::File::CloseClearFStream( read);

        // create histograms
        const storage::Map< std::string, math::Histogram> histograms( score::RadiusOfGyration::HistogramsFromTable( table));

        // iterate over histograms
        for( storage::Map< std::string, math::Histogram>::const_iterator itr( histograms.Begin()), itr_end( histograms.End()); itr != itr_end; ++itr)
        {
          io::OFStream write;
          BCL_ExampleMustOpenOutputFile( write, itr->first + ".histogram");
          write << itr->second;
          io::File::CloseClearFStream( write);
        }
      }

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::score::RadiusOfGyration");
      const std::string given_static_class_name( GetStaticClassName< score::RadiusOfGyration>());
      BCL_Example_Check
      (
        given_static_class_name == correct_static_class_name,
        "GetStaticClassName gives " + given_static_class_name + " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        correct_static_class_name == clone_constr->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test GetScheme
      BCL_MessageStd( "Testing GetScheme");
      BCL_Example_Check
      (
        score_chains.GetScheme() == "rgyr_chain",
        "The scheme for this object is " + util::Format()( score_chains.GetScheme()) + " but should be "
        + "rgyr_chain"
      );

      // test GetHistogramFileName
      BCL_MessageStd( "Testing GetHistogramFileName");
      BCL_Example_Check
      (
        score_chains.GetHistogramFilenameSoluble() == "radius_of_gyration_chains.histogram",
        "The histogram file used for score_chains object is " + util::Format()( score_chains.GetHistogramFilenameSoluble()) +
        " but should be " + "radius_of_gyration_chains.histogram"
      );

      // test GetHistogramFileName
      BCL_Example_Check
      (
        score_model.GetHistogramFilenameSoluble() == "radius_of_gyration_model.histogram",
        " The histogram file used for the score_model object is " + util::Format()( score_model.GetHistogramFilenameSoluble()) +
        " but should be " + "radius_of_gyration_model.histogram"
      );

      // test GetEnergyFunction
      BCL_ExampleCheckWithinTolerance( score_model.GetEnergyFunctionSoluble().operator ()( 0.785), 0.405724, 0.0001);

    ////////////////
    // operations //
    ////////////////

      // test SquareRadiusOfGyration
      BCL_ExampleCheckWithinTolerance( score_model.SquareRadiusOfGyration( *model_ptr), 139.417, 0.0001);

      // test SquareRadiusOfGyrationCollapsed
      BCL_ExampleCheckWithinTolerance( score_model.SquareRadiusOfGyrationCollapsed( membrane_protein, atom_set), 116.012, 0.0001);

    ///////////////
    // operators //
    ///////////////

      // test operator
      BCL_ExampleCheckWithinTolerance( score_model( membrane_protein), -50.1814, 1e-4);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for score::RadiusOfGyration");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( score_model);
      score::RadiusOfGyration read_obj;
      BCL_MessageVrb( "read object");
      ReadBCLObject( read_obj);
      BCL_ExampleCheckWithinTolerance( score_model( *model_ptr), read_obj( *model_ptr), 1e-6);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRadiusOfGyration

  const ExampleClass::EnumType ExampleScoreRadiusOfGyration::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRadiusOfGyration())
  );

} // namespace bcl

