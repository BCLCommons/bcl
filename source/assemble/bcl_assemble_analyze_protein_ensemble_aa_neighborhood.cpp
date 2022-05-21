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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_analyze_protein_ensemble_aa_neighborhood.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_aa_sasa_ols.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_membrane.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_min_max.h"
#include "score/bcl_score_aa_neighborhood_distances.h"
#include "score/bcl_score_aa_neighborhood_exposure.h"
#include "score/bcl_score_aa_pair_clash.h"
#include "score/bcl_score_aa_pair_distance.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! @brief names of allowed scores
    const std::string AnalyzeProteinEnsembleAANeighborhood::s_ScoreNames[ AnalyzeProteinEnsembleAANeighborhood::s_NumberScores] =
    {
      "neighbor_count", "neighbor_vector", "overlapping_spheres", "neighbor_distances", "neighbor_clash"
    };

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeProteinEnsembleAANeighborhood::s_Instance
    (
      util::Enumerated< AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeProteinEnsembleAANeighborhood())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeProteinEnsembleAANeighborhood::AnalyzeProteinEnsembleAANeighborhood() :
      m_OutFilePostFix( ".score.pml")
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeProteinEnsembleAANeighborhood
    AnalyzeProteinEnsembleAANeighborhood *AnalyzeProteinEnsembleAANeighborhood::Clone() const
    {
      return new AnalyzeProteinEnsembleAANeighborhood( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeProteinEnsembleAANeighborhood::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeProteinEnsembleAANeighborhood::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeProteinEnsembleAANeighborhood::GetAlias() const
    {
      static const std::string s_name( "ProteinEnsembleAANeighborhood");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeProteinEnsembleAANeighborhood::operator()( const ProteinEnsemble &ENSEMBLE) const
    {
      std::stringstream script;

      const util::Format model_number_format( util::Format().W( size_t( std::ceil( std::log10( ENSEMBLE.GetSize())))).R().Fill( '0'));
      size_t count( 0);

      std::string first_model_name;

      // iterate over all proteins in the ensemble
      for( ProteinEnsemble::const_iterator ens_itr( ENSEMBLE.Begin()), ens_itr_end( ENSEMBLE.End()); ens_itr != ens_itr_end; ++ens_itr, ++count)
      {
        std::string pymol_model_name( "model" + model_number_format( count));
        // reference on protein model
        const ProteinModel &model( **ens_itr);

        // check for model id
        const util::ShPtr< util::Wrapper< std::string> > sp_id( model.GetProteinModelData()->GetData( ProteinModelData::e_Identification));
        if( sp_id.IsDefined())
        {
          pymol_model_name = *sp_id + model_number_format( count);
        }

        // check if protein model has a membrane
        const util::SiPtr< const biol::Membrane> sp_membrane( model.GetProteinModelData()->GetData( ProteinModelData::e_Membrane));

        // get the filename
        const std::string filename( *util::ShPtr< util::Wrapper< std::string> >( model.GetProteinModelData()->GetData( ProteinModelData::e_PDBFile)));

        // add load to the pymol script
        script << "#color protein model "     << pymol_model_name << '\n';
        script << "load " << filename << ", " << pymol_model_name << '\n';
        script << "hide everything, "         << pymol_model_name << '\n';
        script << "show cartoon, "            << pymol_model_name << "\n\n";

        // create a neighbor list container for the ProteinModel
        const AANeighborListContainer neighbor_list_container( model.GetAminoAcids(), m_Score->GetDistanceCutoff(), m_Score->GetMinimalSequenceSeparation(), true);

        // storage map for amino acids and their score
        storage::Map< util::SiPtr< const biol::AABase>, double> aa_scores;

        // store range of scores
        math::RunningMinMax< double> score_min_max;

        // iterate over all amino neighbor list
        for
        (
          AANeighborListContainer::const_iterator
            aa_neigh_itr( neighbor_list_container.Begin()), aa_neigh_itr_end( neighbor_list_container.End());
          aa_neigh_itr != aa_neigh_itr_end;
          ++aa_neigh_itr
        )
        {
          // score this amino acids neighborlist
          const double score( m_Score->operator ()( aa_neigh_itr->second, sp_membrane));

          if( !util::IsDefined( score))
          {
            continue;
          }

          // consider the score
          aa_scores[ aa_neigh_itr->second.GetCenterAminoAcid()] = score;
          score_min_max += score;
        }

        // range to use
        const math::Range< double> range
        (
          m_ScoreRange.IsEmpty()
          ? math::Range< double>( score_min_max.GetMin(), score_min_max.GetMax())
          : m_ScoreRange
        );

        // set all b factors to half of the score range
        script << "alter " << pymol_model_name << ", b=0.0" << "\n\n";

        // iterate over all scores
        for
        (
          storage::Map< util::SiPtr< const biol::AABase>, double>::const_iterator
            itr( aa_scores.Begin()), itr_end( aa_scores.End());
          itr != itr_end;
          ++itr
        )
        {
          script << "alter " << pymol_model_name
                 << " and chain " << itr->first->GetChainID() << " and resi " << itr->first->GetPdbID()
                 << ", b=" << itr->second << '\n';
        }

        // create a color spectrum over the model
        script << "spectrum b, blue_white_red, " << pymol_model_name << ", " << range.GetMin() << ", " << range.GetMax() << '\n';

        // now color the map based on the underlying protein
        const std::string ramp_name( m_Score->GetScheme() + model_number_format( count));
        script << "ramp_new " << ramp_name << ", " << pymol_model_name << ", [" << range.GetMin() << "," << range.GetMiddle() << "," << range.GetMax() << "], [blue,white,red]\n\n";

        // create a selection of all positively and negatively scoring residues
        const std::string select_p_name( "positive" + model_number_format( count));
        const std::string select_n_name( "negative" + model_number_format( count));
        script << "select " << select_p_name << ", " << pymol_model_name << " and not name c+n+o and b >  0.0001\n";
        script << "select " << select_n_name << ", " << pymol_model_name << " and not name c+n+o and b < -0.0001\n\n";
        script << "group sidechains" << model_number_format( count) << ", " << select_p_name << ' ' << select_n_name << '\n';

        if( count == 0)
        {
          first_model_name = pymol_model_name;
        }
        else
        {
          script << "#disable all for any model except the first\n";
          script << "disable " << pymol_model_name << '\n';
          script << "disable " << ramp_name        << '\n';
          script << "disable sidechains" << model_number_format( count) << '\n';
        }

        script << "\n#end " << pymol_model_name << "\n\n";
      }
      script << "zoom " << first_model_name << '\n';

      // end
      return script.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeProteinEnsembleAANeighborhood::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeProteinEnsembleAANeighborhood::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return a score from the score name
    //! @param SCORE_NAME
    //! @return ShPtr to exposure score
    util::ShPtr< score::AANeighborhoodInterface> AnalyzeProteinEnsembleAANeighborhood::ScoreFromScoreName( const std::string &SCORE_NAME)
    {
      if( SCORE_NAME == s_ScoreNames[ 0])
      {
        return util::CloneToShPtr( score::AANeighborhoodExposure( AANeighborCount()));
      }

      if( SCORE_NAME == s_ScoreNames[ 1])
      {
        return util::CloneToShPtr( score::AANeighborhoodExposure( AANeighborVector()));
      }

      if( SCORE_NAME == s_ScoreNames[ 2])
      {
        return util::CloneToShPtr( score::AANeighborhoodExposure( AASasaOLS()));
      }

      if( SCORE_NAME == s_ScoreNames[ 3])
      {
        return util::CloneToShPtr( score::AANeighborhoodDistances( score::AAPairDistance()));
      }

      if( SCORE_NAME == s_ScoreNames[ 4])
      {
        return util::CloneToShPtr( score::AANeighborhoodDistances( score::AAPairClash()));
      }

      // end
      return util::ShPtr< score::AANeighborhoodInterface>();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeProteinEnsembleAANeighborhood::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        ""
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".score.pml"
      );

      parameters.AddInitializer
      (
        "score",
        "name of the score to be used",
        io::Serialization::GetAgentWithCheck
        (
          &m_ScoreName,
          command::ParameterCheckAllowed( storage::Vector< std::string>( s_NumberScores, s_ScoreNames))
        ),
        s_ScoreNames[ 0]
      );

      parameters.AddInitializer
      (
        "range",
        "the range of over which the residues are colored, if not given, the determined dynamic range over the protein model is used"
        " (must be quoted if range string contains () or ,)",
        io::Serialization::GetAgent( &m_ScoreRange),
        "\"(0.0,0.0)\""
      );

      return parameters;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool AnalyzeProteinEnsembleAANeighborhood::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      m_Score = ScoreFromScoreName( m_ScoreName);
      return true;
    }

  } // namespace assemble
} // namespace bcl
