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
#include "score/bcl_score_sse_membrane_alignment.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "score/bcl_score_protein_model_sse.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_sse_membrane_alignment.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSSEMembraneAlignment :
    public ExampleInterface
  {
  public:

    ExampleScoreSSEMembraneAlignment *Clone() const
    {
      return new ExampleScoreSSEMembraneAlignment( *this);
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
      //instantiate pdb
      io::IFStream read;
      const std::string filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      //instantiate proteinmodels of chains
      const assemble::ProteinModel model( Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone));

      assemble::LocatorSSE locate_helix( 'A', 245, 252);

      util::ShPtrVector< assemble::SSE> spv_helix;
      spv_helix.PushBack( util::ShPtr< assemble::SSE>( locate_helix.Locate( model)->Clone()));
      // move to origin
      spv_helix.FirstElement()->Translate( -spv_helix.FirstElement()->GetCenter());

      util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( model.GetChain( 'A')->GetSequence(), spv_helix.HardCopy()));
      assemble::ProteinModel tmp_model( sp_chain);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      score::SSEMembraneAlignment score_default;

      // construct from membrane object
      util::ShPtr< biol::Membrane> sp_membrane( new biol::Membrane( 5.0, 5.0, 1.0));
      score::SSEMembraneAlignment score_a;

      // apply random transformation, moving both membrane and model (should not change scoring)
      math::TransformationMatrix3D random_transform
      (
        coord::OrientationInterface::GenerateRandomTransformationAroundCenter
        (
          10.0,
          math::g_Pi / 2.0,
          linal::Vector3D( 0.0, 0.0, 0.0)
        )
      );
      spv_helix.FirstElement()->Transform( random_transform);
      sp_membrane->Transform( random_transform);
      tmp_model.Transform( random_transform);

      // add membrane to model
      util::ShPtr< assemble::ProteinModelData> sp_data( tmp_model.GetProteinModelData());
      sp_data->Insert( assemble::ProteinModelData::e_Membrane, sp_membrane);
      tmp_model.SetProteinModelData( sp_data);

      // clone
      util::ShPtr< util::ObjectInterface> ptr( score_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // test cloned
      BCL_ExampleIndirectCheck( ptr->GetClassIdentifier(), GetStaticClassName< score::SSEMembraneAlignment>(), "Clone");

      // test scheme
      BCL_Example_Check
      (
        score_a.GetScheme() == "ssealign",
        "The scheme should be ssealign but instead it is : " + score_a.GetScheme()
      );

      // test default histogram filename
      BCL_Example_Check
      (
        score::SSEMembraneAlignment::GetDefaultHistogramFilename() == "sse_membrane_alignment.histograms",
        "The default histogram filename is different: " + score::SSEMembraneAlignment::GetDefaultHistogramFilename()
      );

      // test default scheme
      BCL_Example_Check
      (
        score::SSEMembraneAlignment::GetDefaultScheme() == "ssealign",
        "The default scheme should be ssealign but instead it is : " + score::SSEMembraneAlignment::GetDefaultScheme()
      );

      // access to energy function
      BCL_Example_Check
      (
        score_a.GetEnergyFunctions().GetSize() == 2,
        "there should be 2 energy functions maps (helix and strand) but found: " +
        util::Format()( score_a.GetEnergyFunctions().GetSize())
      );

      // access to energy function
      BCL_Example_Check
      (
        score_a.GetEnergyFunctions().GetValue( biol::GetSSTypes().STRAND).GetSize() == 3,
        "there should be 3 spline functions for strand for three environments but found: " +
        util::Format()( score_a.GetEnergyFunctions().GetValue( biol::GetSSTypes().STRAND).GetSize())
      );

    ////////////////
    // operations //
    ////////////////

      // angle to membrane plane
      const double expected_angle( 0.569801);
      const double calculated_angle( score_a.AngleToMembranePlane( *spv_helix.FirstElement(), coord::GetAxes().e_Z, *sp_membrane));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_angle, calculated_angle),
        "calculated angle of sse z-axis to membrane plane is incorrect: " + util::Format()( calculated_angle) + " != " +
        util::Format()( expected_angle)
      );

      // WeightXAxis
      const double expected_weight_helix( 1.0);
      const double calculated_weight_helix( score_a.WeightXAxis( *spv_helix.FirstElement(), *sp_membrane));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_weight_helix, calculated_weight_helix),
        "calculated weight for HELIX x-axis is incorrect: " + util::Format()( calculated_weight_helix) + " != " +
        util::Format()( expected_weight_helix)
      );

    ///////////////
    // operators //
    ///////////////

      const double shift( 2.75);
      const double rotation( math::g_Pi / 5);

      // operator
      const storage::Pair< double, size_t> expected_score_nr1( 9.05445, 4);
      const storage::Pair< double, size_t> calculated_score_nr1( score_a( *spv_helix.FirstElement(), *sp_membrane));
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_score_nr1.First(), calculated_score_nr1.First())
        && expected_score_nr1.Second() == calculated_score_nr1.Second(),
        "calculated score and nr entities is wrong: " + util::Format()( calculated_score_nr1) + " != " +
        util::Format()( expected_score_nr1)
      );

      for( size_t trans( 1); trans < 15; ++trans)
      {
        spv_helix.FirstElement()->Translate( linal::Vector3D( 0.0, 0.0, shift));

        for( size_t angle( 1); angle <= 10; ++angle)
        {
          spv_helix.FirstElement()->Rotate( math::RotationMatrix3D( coord::GetAxes().e_Y, rotation));

          BCL_MessageStd
          (
            "trans: " + util::Format().W( 8).FFP( 3)( trans * shift) +
            " rotation: " + util::Format().W( 8).FFP( 3)( math::Angle::Degree( angle * rotation)) +
            " score: " + util::Format()( score::ProteinModelSSE( util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >( score_default.Clone()), false)( tmp_model))
          );
        }
      }

      BCL_ExampleCheckWithinAbsTolerance( score_a( *spv_helix.FirstElement(), *sp_membrane).First(), 7.81834, 0.0001);

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( score_default);
      // read from file
      score::SSEMembraneAlignment score_read;
      ReadBCLObject( score_read);

      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          score_default( *spv_helix.FirstElement(), *sp_membrane).First(),
          score_read( *spv_helix.FirstElement(), *sp_membrane).First()
        ),
        "written sse membrane alignment score returns a different score"
      );

      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // iterate over sse types
        for
        (
          biol::SSTypes::const_iterator
            ss_type_itr( biol::GetSSTypes().Begin()), ss_type_itr_end( biol::GetSSTypes().COIL.GetIterator());
          ss_type_itr != ss_type_itr_end;
          ++ss_type_itr
        )
        {
          // write function to gnuplot file
          io::OFStream write;
          BCL_ExampleMustOpenOutputFile( write, ss_type_itr->GetName() + "_membrane_alignment.gnuplot");
          util::SiPtrVector< const math::CubicSplineDamped> energy_functions;
          storage::Vector< std::string> descriptors;
          // iterate over environment type to collect energy functions
          for
          (
            auto env_itr( score_default.GetEnergyFunctions().GetValue( *ss_type_itr).Begin()),
              env_itr_end( score_default.GetEnergyFunctions().GetValue( *ss_type_itr).End());
            env_itr != env_itr_end;
            ++env_itr
          )
          {
            energy_functions.PushBack( &env_itr->second);
            descriptors.PushBack( env_itr->first.GetName());
          }
          math::GnuplotHeatmap heatmap;
          heatmap.SetFromCubicSplines( energy_functions, false, false);
          heatmap.SetTitleAndLabel( "angle potential " + ss_type_itr->GetName() + " main axis to membrane normal", "angle [rad]", "", "-log(p)");
          heatmap.SetFont( "arialbd", 16);
          heatmap.SetTicsY( descriptors, true, 1);
          heatmap.SetFilename( "membrane_alignment_" + ss_type_itr->GetName());
          heatmap.WriteScript( write);
          io::File::CloseClearFStream( write);
        }
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreSSEMembraneAlignment

  const ExampleClass::EnumType ExampleScoreSSEMembraneAlignment::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSSEMembraneAlignment())
  );

} // namespace bcl
