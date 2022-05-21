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
#include "score/bcl_score_contact_order.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_contact_order.cpp
  //!
  //! @author woetzen
  //! @date Apr 2, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreContactOrder :
    public ExampleInterface
  {
  public:

    ExampleScoreContactOrder *Clone() const
    {
      return new ExampleScoreContactOrder( *this);
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
      score::ContactOrder score_model_rel_aas( contact::Order::e_RelativeAAsUsed, false);
      score::ContactOrder score_model_rel_seq( contact::Order::e_RelativeSequenceLength, false);
      util::ShPtr< score::ContactOrder> clone_constr( score_model_rel_aas.Clone());

      storage::Vector< score::ContactOrder> contact_order_scores;
      for( size_t i( contact::Order::e_RelativeAAsUsed); i <= contact::Order::e_RelativeSqrSequenceLength; ++i)
      {
        contact_order_scores.PushBack( score::ContactOrder( contact::Order::NormalizationType( i), false));
      }

      // create protein model
      const std::string pdb( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      assemble::ProteinModel this_model( Proteins::GetModel( pdb, biol::GetAAClasses().e_AACaCb));
      storage::Map< biol::SSType, size_t> sse_min_size;
      sse_min_size[ biol::GetSSTypes().HELIX] = 5;
      sse_min_size[ biol::GetSSTypes().STRAND] = 3;
      assemble::ProteinModel this_model_sse( Proteins::GetModel( pdb, biol::GetAAClasses().e_AACaCb, sse_min_size));

      // construct splines
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // iterate through scores
        for( storage::Vector< score::ContactOrder>::const_iterator itr( contact_order_scores.Begin()), itr_end( contact_order_scores.End()); itr != itr_end; ++itr)
        {
          const score::ContactOrder &score_model( *itr);
          const std::string &normalization( contact::Order::GetNormalizationTypeDescriptor( score_model.GetContactOrderFunction().GetNormalization()));
          // output the gnuplot heatmap files
          io::OFStream write;
          {
            const linal::Vector< double> binning( score_model.GetEnergyFunction()->GetXValues());
            BCL_ExampleMustOpenOutputFile( write, "co_order_" + normalization + ".gnuplot");
            math::GnuplotHeatmap heatmap;
            heatmap.SetFromCubicSpline( *score_model.GetEnergyFunction(), false, false);
            heatmap.SetTitleAndLabel( "contact order potential " + normalization, "contact order / #AAs", "", "-log(p)");
            heatmap.SetPixelAndRatio( 1080, 800, -1.0);
            heatmap.SetTicsX( math::GnuplotHeatmap::TicsFromBinning( binning, 1, util::Format().FFP( 2).W( 4).Fill( ' ').R()), false, 1);
            heatmap.SetRotationXTics( 90);
            heatmap.SetFont( "arialbd", 16);
            heatmap.SetFilename( "co_order_" + normalization);
            heatmap.WriteScript( write);
            io::File::CloseClearFStream( write);
          }
        }

        // load table of raw contact order statistics
        const std::string filename( AddExampleInputPathToFilename( e_Biology, "contact_order.table"));
        io::IFStream read;
        BCL_ExampleMustOpenInputFile( read, filename);
        storage::Table< double> table;
        table.ReadFormatted( read);
        io::File::CloseClearFStream( read);

        // create histograms
        const storage::Map< std::string, math::Histogram> histograms( score::ContactOrder::HistogramsFromColumns( table));

        // iterate over histograms
        for( storage::Map< std::string, math::Histogram>::const_iterator itr( histograms.Begin()), itr_end( histograms.End()); itr != itr_end; ++itr)
        {
          if( itr->first == "nr_aas" || itr->first == "nr_aas_sses")
          {
            continue;
          }
          io::OFStream write;

          BCL_ExampleMustOpenOutputFile( write, itr->first + ".histogram");
          write << itr->second;
          io::File::CloseClearFStream( write);
        }
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< score::ContactOrder>(), clone_constr->GetClassIdentifier());

      // test GetScheme
      BCL_MessageStd( "Testing GetScheme");
      BCL_ExampleCheck( score_model_rel_aas.GetScheme(), "co_score");

      // test GetHistogramFileName
      BCL_MessageStd( "Testing GetHistogramFileName");
      BCL_ExampleCheck
      (
        score_model_rel_aas.GetHistogramFilename(),
        score::ContactOrder::GetDefaultHistogramFilename( contact::Order::e_RelativeAAsUsed)
      );
      BCL_ExampleCheck
      (
        score_model_rel_seq.GetHistogramFilename(),
        score::ContactOrder::GetDefaultHistogramFilename( contact::Order::e_RelativeSequenceLength)
      );

      // test GetEnergyFunction
      BCL_MessageStd( "Testing GetEnergyFunction");
      BCL_ExampleCheckWithinTolerance( score_model_rel_aas.GetEnergyFunction()->operator ()( 0.25), -0.935785, 0.0001);
      BCL_ExampleCheckWithinTolerance( score_model_rel_seq.GetEnergyFunction()->operator ()( 0.25), -0.78382, 0.0001);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // test operator
      BCL_MessageStd( "Testing Operator");
      BCL_ExampleCheckWithinTolerance( score_model_rel_aas( this_model), -105.156, 0.0001);
      BCL_ExampleCheckWithinTolerance( score_model_rel_seq( this_model_sse), 11.231, 0.0001);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for score::ContactOrder");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( score_model_rel_aas);
      score::ContactOrder read_obj( contact::Order::e_RelativeSequenceLength);
      BCL_MessageVrb( "read object");
      ReadBCLObject( read_obj);
      BCL_ExampleIndirectCheck( score_model_rel_aas( this_model), read_obj( this_model), "read and write");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreContactOrder

  const ExampleClass::EnumType ExampleScoreContactOrder::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreContactOrder())
  );

} // namespace bcl

