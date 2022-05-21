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
#include "restraint/bcl_restraint_analyze_atom_distance_score_heatmap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "score/bcl_score_restraint_distance_spin_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceScoreHeatmap::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceScoreHeatmap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceScoreHeatmap::AnalyzeAtomDistanceScoreHeatmap() :
      m_Score(),
      m_OutFilePostFix( ".AtomDistanceScoreHeatmap"),
      m_HistogramMinimum( -1.0),
      m_HistogramBinSize( 0.1),
      m_HistogramNumberOfBins( 10),
      m_RestraintType( EPRDistanceData::GetDefaultHandler())
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceScoreHeatmap
    AnalyzeAtomDistanceScoreHeatmap *AnalyzeAtomDistanceScoreHeatmap::Clone() const
    {
      return new AnalyzeAtomDistanceScoreHeatmap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceScoreHeatmap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceScoreHeatmap::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceScoreHeatmap::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceScoreHeatmap");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistanceScoreHeatmap::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> data( m_RestraintType->ReadRestraintsFromFile());

      const storage::Vector< math::Histogram> score_histograms( GetScoreHistograms( data, ENSEMBLE));

      const storage::Vector< std::string> restraint_tics( GetRestraintNameTics( data));

      math::GnuplotHeatmap heatmap( GetHeatMap( score_histograms, restraint_tics));

      // write heat map to string stream
      std::stringstream stream;
      heatmap.WriteScript( stream);

      // return the string
      return stream.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistanceScoreHeatmap::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistanceScoreHeatmap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceScoreHeatmap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "creates heat map showing how frequently a given score is achieved by each atom distance restraint"
      );

      parameters.AddInitializer
      (
        "restraint_score_type",
        "the atom distance restraint type where the score will come from",
        io::Serialization::GetAgent( &m_Score),
        score::RestraintDistanceSpinLabel::GetDefaultScheme()
      );

      parameters.AddInitializer
      (
        "restraint_type",
        "the restraint type; used to load the actual restraints",
        io::Serialization::GetAgent( &m_RestraintType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScoreHeatmap"
      );

      parameters.AddInitializer
      (
        "histogram_minimum",
        "the minimal value representing the left boundary of the score histogram",
        io::Serialization::GetAgent( &m_HistogramMinimum),
        "-1.0"
      );

      parameters.AddInitializer
      (
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "histogram_num_bins",
        "the number of bins in the score histograms",
        io::Serialization::GetAgent( &m_HistogramNumberOfBins),
        "10"
      );

      return parameters;
    }

    //! @brief creates histograms of how frequently a score is achieved by each restraint
    //! @param DATA the restraints that will be scored
    //! @param ENSEMBLE the ensemble of models that will be scored according to the atom distance restraints
    //! @return vector of histograms - one histogram for each restraint + one for average sum
    storage::Vector< math::Histogram> AnalyzeAtomDistanceScoreHeatmap::GetScoreHistograms
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // takes into account adding columns for the average score
      static const size_t s_additional_columns( 1);

      // hold the score data for each restraint and score sum and total score
      storage::Vector< math::Histogram> score_histograms
      (
        DATA.GetSize() + s_additional_columns,
        math::Histogram( m_HistogramMinimum, m_HistogramBinSize, m_HistogramNumberOfBins)
      );

      // score to score the individual restraints
      const math::FunctionInterfaceSerializable< AtomDistanceAssignment, double> &individual_score( *m_Score->GetScore());

      const double data_size( DATA.GetSize());

      // iterate through the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator protein_itr( ENSEMBLE.Begin()), protein_itr_end( ENSEMBLE.End());
        protein_itr != protein_itr_end;
        ++protein_itr
      )
      {
        const assemble::ProteinModel &model( **protein_itr);

        // iterators for the statistics vectors
        storage::Vector< math::Histogram>::iterator
          score_histograms_itr( score_histograms.Begin()), score_histograms_itr_end( score_histograms.End());

        // to hold sum of restraint scores for this protein
        double score_sum( 0);

        // iterate through the atom distance restraints
        for
        (
          util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
          restraint_itr != restraint_itr_end && score_histograms_itr != score_histograms_itr_end;
          ++restraint_itr, ++score_histograms_itr
        )
        {
          // get the assignment
          const AtomDistanceAssignment assignment( ( *restraint_itr)->GenerateAssignment( model));

          // get the current score
          const double score( individual_score( assignment));

          score_sum += score;

          // add the score to the two statistics objects
          score_histograms_itr->PushBack( score);
        } // end iterate through restraints

        // add the score sum and total score to the two histograms objects
        score_histograms_itr->PushBack( score_sum / data_size);
      } // end iterate through ensemble

      // normalize all the histograms
      for
      (
        storage::Vector< math::Histogram>::iterator
          histogram_itr( score_histograms.Begin()), histogram_itr_end( score_histograms.End());
        histogram_itr != histogram_itr_end;
        ++histogram_itr
      )
      {
        histogram_itr->Normalize();
        BCL_MessageDbg( "histogram " + util::Format()( histogram_itr_end - histogram_itr) + "\n" + util::Format()( *histogram_itr));
      }

      return score_histograms;
    }

    //! @brief creates the tics that will be used in the heatmap - converts each restraint into a string to use as tic
    //! @param DATA the restraints that will be used to create tics
    //! @return vector of strings which are the tics representing each restraint
    storage::Vector< std::string> AnalyzeAtomDistanceScoreHeatmap::GetRestraintNameTics
    (
      const util::ShPtrVector< AtomDistance> &DATA
    )
    {
      storage::Vector< std::string> tics;

      // iterate through the data to get the tics from the restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_a( *( *data_itr)->GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_b( *( *data_itr)->GetData().Second());
        tics.PushBack( assemble::LocatorAtomCoordinatesInterface::GetNameFromPair( atom_locator_a, atom_locator_b));
      }

      tics.PushBack( "mean");

      return tics;
    }

    //! @brief creates gnuplot heat map object from the score histograms and the restraint names as tics
    //! @param SCORE_HISTOGRAMS the histograms that will be used to make the heat map
    //! @param TICS the names of the restraints
    //! @return heat map object which represents the distribution of scores for each restraint
    math::GnuplotHeatmap AnalyzeAtomDistanceScoreHeatmap::GetHeatMap
    (
      const storage::Vector< math::Histogram> &SCORE_HISTOGRAMS, const storage::Vector< std::string> &TICS
    )
    {
      const util::SiPtrVector< const math::Histogram> histograms
      (
        util::ConvertToConstSiPtrVector( SCORE_HISTOGRAMS)
      );

      math::GnuplotHeatmap heatmap;

      heatmap.SetFromHistograms( histograms, true, true);
      heatmap.SetTitleAndLabel( "Title", "restraint", "score", "Frequency (fraction of models)");
      heatmap.SetTicsX( TICS, true, 1);
      heatmap.SetRotationXTics( 90);
      heatmap.SetFilename( "AnalyzeAtomDistanceScoreHeatmap.gnuplot");
      heatmap.SetFont( "/usr/share/fonts/dejavu-lgc/DejaVuLGCSansMono.ttf", 10);

      return heatmap;
    }

  } // namespace restraint
} // namespace bcl
