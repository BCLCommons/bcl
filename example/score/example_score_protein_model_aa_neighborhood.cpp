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
#include "score/bcl_score_protein_model_aa_neighborhood.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "score/bcl_score_aa_neighborhood_exposure.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_exposures.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreProteinModelAANeighborhood :
    public ExampleInterface
  {
  public:

    ExampleScoreProteinModelAANeighborhood *Clone() const
    {
      return new ExampleScoreProteinModelAANeighborhood( *this);
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

    //! write a gnuplot heatmap for soluble exposure potentials
    void WriteGnuplot
    (
      const storage::Map< biol::AAType, math::CubicSplineDamped> &MAP,
      const std::string &IDENTIFIER,
      std::ostream &OSTREAM
    ) const
    {
      util::SiPtrVector< const math::CubicSplineDamped> splines;
      storage::Vector< std::string> spline_descriptors;

      // iterate over all aa types
      for
      (
        storage::Map< biol::AAType, math::CubicSplineDamped>::const_iterator
          itr( MAP.Begin()), itr_end( MAP.End());
        itr != itr_end;
        ++itr
      )
      {
        splines.PushBack( itr->second);
        spline_descriptors.PushBack( itr->first->GetThreeLetterCode());
      }

      // write splines to gnuplot file
      math::GnuplotHeatmap heatmap;
      const math::CubicSplineDamped &template_spline( *splines.FirstElement());
      heatmap.SetFromFunctions
      (
        util::SiPtrVector< const math::FunctionInterfaceSerializable< double, double> >( splines),
        18 * 4,
        template_spline.GetXValues().First(),
        template_spline.GetDelta() / 4.0,
        false,
        false,
        false
      );
      heatmap.SetTitleAndLabel( "amino acid " + IDENTIFIER + " potentials", IDENTIFIER, "", "-log(p/b)");
      heatmap.SetPixelAndRatio( 800, 800, -1);
      heatmap.SetFont( "arialbd", 16);
      heatmap.SetTicsY( spline_descriptors, true, 1);
      heatmap.SetRotationXTics( 90.0);
      heatmap.SetFilename( IDENTIFIER + "_potentials");
      heatmap.WriteScript( OSTREAM);
    }

    int Run() const
    {
      // instantiate pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      //instantiate proteinmodels of chains
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 4;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;
      assemble::ProteinModel model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      util::ShPtrVector< assemble::SSE> helixpair;
      for
      (
        storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
          sse_itr( model.GetChains()( 0)->GetData().Begin()), sse_itr_end( model.GetChains()( 0)->GetData().End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        if( ( *sse_itr)->GetFirstAA()->GetPdbID() == 245 || ( *sse_itr)->GetFirstAA()->GetPdbID() == 267)
        {
          helixpair.PushBack( *sse_itr);
          helixpair.LastElement()->SetToIdealConformationAtOrigin();
        }
      }

      const double shift( 1.0);
      const double rotation( math::g_Pi / 2);

      util::ShPtr< biol::Membrane> membrane( new biol::Membrane( 10.0, 5.0, 5.0));

      util::ShPtr< score::AANeighborhoodExposure> sp_neighbor_count
      (
        new score::AANeighborhoodExposure( assemble::AANeighborCount())
      );
      util::ShPtr< score::AANeighborhoodExposure> sp_neighbor_vector
      (
        new score::AANeighborhoodExposure( assemble::AANeighborVector())
      );
      util::ShPtr< score::AANeighborhoodExposure> sp_ols_sasa
      (
        new score::AANeighborhoodExposure( assemble::AASasaOLS())
      );

      const score::ProteinModel::Type score_type( score::ProteinModel::e_Sequence);
      const std::string readable_scheme( "AA Environment");

      storage::List< score::ProteinModelAANeighborhood> exposure_scores;
      exposure_scores.PushBack
      (
        score::ProteinModelAANeighborhood
        (
          sp_neighbor_count,
          score::ProteinModelAANeighborhood::e_None,
          true,
          score_type,
          readable_scheme
        )
      );
      exposure_scores.PushBack( score::ProteinModelAANeighborhood( sp_neighbor_vector));
      exposure_scores.PushBack( score::ProteinModelAANeighborhood( sp_ols_sasa));

      BCL_ExampleCheck( exposure_scores.FirstElement().GetType(), score_type);
      BCL_ExampleCheck( exposure_scores.FirstElement().GetReadableScheme(), readable_scheme);

      for( size_t trans( 1); trans < 10; trans += 3)
      {
        BCL_MessageStd( "translational position: " + util::Format().W( 8).FFP( 3)( trans * shift));
        helixpair.LastElement()->Translate( linal::Vector3D( 0.0, 0.0, shift));

        for( size_t angle( 1); angle <= 4; ++angle)
        {
          BCL_MessageStd
          (
            "rotation angle: " + util::Format().W( 8).FFP( 3)( math::Angle::Degree( angle * rotation))
          );
          helixpair.LastElement()->Rotate( math::RotationMatrix3D( coord::GetAxes().e_Y, rotation));

          util::ShPtr< assemble::Chain> sp_chain( new assemble::Chain( model.GetChains()( 0)->GetSequence(), helixpair));
          assemble::ProteinModel tmp_model( sp_chain);

          // score soluble protein
          for
          (
            storage::List< score::ProteinModelAANeighborhood>::const_iterator
              score_itr( exposure_scores.Begin()), score_itr_end( exposure_scores.End());
            score_itr != score_itr_end; ++score_itr)
          {
            BCL_MessageStd
            (
              "soluble:  " + score_itr->GetScheme() + ":\t" + util::Format()( score_itr->operator()( tmp_model))
            );
          }
          assemble::ProteinModel tmp_model_membrane( sp_chain);
          util::ShPtr< assemble::ProteinModelData> sp_data( tmp_model_membrane.GetProteinModelData());
          sp_data->Insert( assemble::ProteinModelData::e_Membrane, membrane);
          tmp_model_membrane.SetProteinModelData( sp_data);

          // score membrane protein
          for
          (
            storage::List< score::ProteinModelAANeighborhood>::const_iterator
              score_itr( exposure_scores.Begin()), score_itr_end( exposure_scores.End());
            score_itr != score_itr_end; ++score_itr
          )
          {
            BCL_MessageStd
            (
              "membrane: " + score_itr->GetScheme() + ":\t" + util::Format()( score_itr->operator()( tmp_model_membrane))
            );
          }
        }
      }

      // print gnuplot heat map
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, sp_neighbor_count->GetScheme() + "_potentials.gnuplot");
        WriteGnuplot( sp_neighbor_count->GetSolublePotentials(), sp_neighbor_count->GetScheme(), write);
        BCL_ExampleMustOpenOutputFile( write, sp_neighbor_vector->GetScheme() + "_potentials.gnuplot");
        WriteGnuplot( sp_neighbor_vector->GetSolublePotentials(), sp_neighbor_vector->GetScheme(), write);
        BCL_ExampleMustOpenOutputFile( write, sp_ols_sasa->GetScheme() + "_potentials.gnuplot");
        WriteGnuplot( sp_ols_sasa->GetSolublePotentials(), sp_ols_sasa->GetScheme(), write);

        // iterate over environments
        biol::GetEnvironmentTypes().GetReducedTypes();
        for
        (
          biol::EnvironmentTypes::const_iterator env_itr( biol::GetEnvironmentTypes().GetReducedTypes().Begin()),
            env_itr_end( biol::GetEnvironmentTypes().GetReducedTypes().End());
          env_itr != env_itr_end; ++env_itr
        )
        {
          BCL_ExampleMustOpenOutputFile( write, env_itr->GetName() + "_potentials.gnuplot");
          WriteGnuplot( sp_neighbor_count->GetMembranePotentials( *env_itr), env_itr->GetName(), write);
        }
        io::File::CloseClearFStream( write);
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreProteinModelAANeighborhood

  const ExampleClass::EnumType ExampleScoreProteinModelAANeighborhood::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelAANeighborhood())
  );

} // namespace bcl

