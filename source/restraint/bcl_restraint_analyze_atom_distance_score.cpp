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
#include "restraint/bcl_restraint_analyze_atom_distance_score.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "restraint/bcl_restraint_analyze_atom_distance_mean_sd.h"
#include "restraint/bcl_restraint_epr_distance_data.h"
#include "score/bcl_score_restraint_distance_spin_label.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceScore::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceScore())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceScore::AnalyzeAtomDistanceScore() :
      m_Score(), // EPRDistanceData().GetAlias()
      m_PrintAllModelScores( true),
      m_OutFilePostFix( ".AnalyzeAtomDistanceScore"),
      m_PrintStatistics( true),
      m_RestraintType( EPRDistanceData::GetDefaultHandler())
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceScore
    AnalyzeAtomDistanceScore *AnalyzeAtomDistanceScore::Clone() const
    {
      return new AnalyzeAtomDistanceScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceScore::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceScore::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceScore");
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
    std::string AnalyzeAtomDistanceScore::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> data( m_RestraintType->ReadRestraintsFromFile());

      // format for the line names
      static const size_t s_width( 17);
      util::Format line_name_format( util::Format().W( s_width));

      // true if scores for each model will be printed -  need to get length of pdb filename
      if( m_PrintAllModelScores)
      {
        // get the name of the first model in order to determine its length
        util::ShPtr< util::Wrapper< std::string> > pdb_name
        (
          ( *ENSEMBLE.Begin())->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );
        BCL_Assert( pdb_name.IsDefined(), "could not get pdb name from protein model");

        line_name_format = util::Format().W( pdb_name->length());
      }

      static const util::Format s_format( util::Format().W( s_width).FFP( 3));
      // add the first line of the analysis
      {
        analysis += AnalyzeAtomDistanceMeanSD::GetRestraintHeader( data, line_name_format, s_format);
        analysis += ( s_format( "sum") + s_format( "total_score") + "\n");
      }

      // add the score analysis lines
      analysis += GetScoreAnalysis( data, ENSEMBLE, line_name_format, s_format);

      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeAtomDistanceScore::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Score, ISTREAM);
      io::Serialize::Read( m_PrintAllModelScores, ISTREAM);
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_PrintStatistics, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeAtomDistanceScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Score, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_PrintAllModelScores, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_PrintStatistics, OSTREAM, INDENT) ;

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceScore::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates atom distance scores for ensemble of models and gives mean and sd of score of each restraint."
        "Can also optionally write the individual restraint score for every model in the ensemble"
      );

      parameters.AddInitializer
      (
        "score",
        "the atom distance score",
        io::Serialization::GetAgent( &m_Score),
        score::RestraintDistanceSpinLabel::GetDefaultScheme()
      );

      parameters.AddInitializer
      (
        "print_all_scores",
        "one if scores should be printed for every model - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintAllModelScores),
        "1"
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AnalyzeAtomDistanceScore"
      );

      parameters.AddInitializer
      (
        "print_statistics",
        "one if the mean, std dev, min, and max statistics should be printed for each atom distance restraint - 0 otherwise",
        io::Serialization::GetAgent( &m_PrintStatistics),
        "1"
      );

      parameters.AddInitializer
      (
        "restraint_type",
        "actual type of restraint that holds the data",
        io::Serialization::GetAgent( &m_RestraintType),
        EPRDistanceData().GetAlias()
      );
      return parameters;
    }

    //! @brief gives the string that is the scoring analysis
    //! @param DATA the list of atom distance objects which are the restraints
    //! @param ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    std::string AnalyzeAtomDistanceScore::GetScoreAnalysis
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    ) const
    {
      std::string analysis_string;

      // get statistics of atom distance restraint scores and add formatted analysis for each model if desired
      storage::Pair
      <
        storage::Vector< math::RunningAverageSD< double> >, storage::Vector< math::RunningMinMax< double> >
      > statistics( GetScoreStatistics( DATA, ENSEMBLE, LINE_NAME_FORMAT, FORMAT, analysis_string));

      if( m_PrintStatistics)
      {
        analysis_string +=
        (
          "\n" + GetScoreStatisticsFormattedAnalysis
          (
            statistics.First(), statistics.Second(), LINE_NAME_FORMAT, FORMAT
          )
        );
      }

      return analysis_string;
    }

    //! @brief gives statistics for each atom distance restraint and add individual model analysis if necessary
    //! @param DATA the list of atom distance objects which are the restraints
    //! @PARAM ENSEMBLE the ensemble from which distance mean and std devs will be calculated
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @param ANALYSIS_STRING string holding formatted analysis
    //! @return pair of vectors with a statistic object for each atom distance restraint
    storage::Pair
    <
      storage::Vector< math::RunningAverageSD< double> >, storage::Vector< math::RunningMinMax< double> >
    > AnalyzeAtomDistanceScore::GetScoreStatistics
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT, std::string &ANALYSIS_STRING
    ) const
    {
      // takes into account adding columns for the total score and the simple sum of the individual restraint scores
      static const size_t s_additional_columns( 2);

      // to hold mean and sd of scores
      storage::Vector< math::RunningAverageSD< double> > mean_sd
      (
        DATA.GetSize() + s_additional_columns, math::RunningAverageSD< double>()
      );

      // to hold min and max of scores
      storage::Vector< math::RunningMinMax< double> > min_max
      (
        DATA.GetSize() + s_additional_columns, math::RunningMinMax< double>()
      );

      // score to score the individual restraints
      const math::FunctionInterfaceSerializable< AtomDistanceAssignment, double> &individual_score( *m_Score->GetScore());

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
        storage::Vector< math::RunningAverageSD< double> >::iterator
          mean_sd_itr( mean_sd.Begin()), mean_sd_itr_end( mean_sd.End());
        storage::Vector< math::RunningMinMax< double> >::iterator
          min_max_itr( min_max.Begin()), min_max_itr_end( min_max.End());

        // true if scores for each model will be printed -  need to print of pdb filename
        if( m_PrintAllModelScores)
        {
          // get the name of the first model in order to determine its length
          util::ShPtr< util::Wrapper< std::string> > pdb_name
          (
            ( *protein_itr)->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
          );
          BCL_Assert( pdb_name.IsDefined(), "could not get atom distance restraints from protein model");
          if( protein_itr != ENSEMBLE.Begin())
          {
            ANALYSIS_STRING += "\n";
          }
          ANALYSIS_STRING += LINE_NAME_FORMAT( std::string( *pdb_name));
        }

        // to hold sum of restraint scores for this protein
        double score_sum( 0);

        // iterate through the atom distance restraints
        for
        (
          util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
          restraint_itr != restraint_itr_end && mean_sd_itr != mean_sd_itr_end && min_max_itr != min_max_itr_end;
          ++restraint_itr, ++mean_sd_itr, ++min_max_itr
        )
        {
          // get the assignment
          const AtomDistanceAssignment assignment( ( *restraint_itr)->GenerateAssignment( model));

          // get the current score
          const double score( individual_score( assignment));

          score_sum += score;

          // add the score to the two statistics objects
          *mean_sd_itr += score;
          *min_max_itr += score;

          // true if scores for each model will be printed -  need to print current score
          if( m_PrintAllModelScores)
          {
            ANALYSIS_STRING += FORMAT( score);
          }
        } // end iterate through restraints

        // add the score sum and total score to the two statistics objects
        const double score_total( m_Score->operator()( model));
        *mean_sd_itr += score_sum;
        *++mean_sd_itr += score_total;
        *min_max_itr += score_sum;
        *++min_max_itr += score_total;

        // true if scores for each model will be printed -  need to print score sum and total score
        if( m_PrintAllModelScores)
        {
          ANALYSIS_STRING += FORMAT( score_sum);
          ANALYSIS_STRING += FORMAT( score_total);
        }
      } // end iterate through ensemble

      return storage::Pair
        <
          storage::Vector< math::RunningAverageSD< double> >, storage::Vector< math::RunningMinMax< double> >
        >( mean_sd, min_max);
    }

    //! @brief gives formatted analysis of the statistics of the atom distance restraints
    //! @param MEAN_SD the mean and standard deviation statistics of the atom distance restraints
    //! @param MIN_MAX the min and max statistics of the atom distance restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @return string which has the formatted statistics analysis of atom restraint scores
    std::string AnalyzeAtomDistanceScore::GetScoreStatisticsFormattedAnalysis
    (
      const storage::Vector< math::RunningAverageSD< double> > &MEAN_SD,
      const storage::Vector< math::RunningMinMax< double> > &MIN_MAX,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    ) const
    {
      std::string analysis;
      analysis += ( GetMeanSDLines( MEAN_SD, LINE_NAME_FORMAT, FORMAT) + "\n");
      analysis += GetMinMaxLines( MIN_MAX, LINE_NAME_FORMAT, FORMAT);

      return analysis;
    }

    //! @brief gives formatted analysis of the mean and sd statistics of the atom distance restraints
    //! @param MEAN_SD the mean and standard deviation statistics of the atom distance restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @return string which has the formatted mean and sd statistics analysis of atom restraint scores
    std::string AnalyzeAtomDistanceScore::GetMeanSDLines
    (
      const storage::Vector< math::RunningAverageSD< double> > &MEAN_SD,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {
      std::string mean_line( LINE_NAME_FORMAT( "Mean"));
      std::string sdev_line( LINE_NAME_FORMAT( "StdDev"));

      // iterate through the mean and std dev statistics
      for
      (
        storage::Vector< math::RunningAverageSD< double> >::const_iterator
          stat_itr( MEAN_SD.Begin()), stat_itr_end( MEAN_SD.End());
        stat_itr != stat_itr_end; ++stat_itr
      )
      {
        mean_line += FORMAT( stat_itr->GetAverage());
        sdev_line += FORMAT( stat_itr->GetStandardDeviation());
      }

      std::string analysis( mean_line + "\n" + sdev_line);

      return analysis;
    }

    //! @brief gives formatted analysis of the min and max statistics of the atom distance restraints
    //! @param MIN_MAX the min and max statistics of the atom distance restraints
    //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
    //! @param FORMAT the format object to format the scores
    //! @return string which has the formatted min and max statistics analysis of atom restraint scores
    std::string AnalyzeAtomDistanceScore::GetMinMaxLines
    (
      const storage::Vector< math::RunningMinMax< double> > &MIN_MAX,
      const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
    )
    {
      std::string min_line( LINE_NAME_FORMAT( "Min"));
      std::string max_line( LINE_NAME_FORMAT( "Max"));

      // iterate through the mean and std dev statistics
      for
      (
        storage::Vector< math::RunningMinMax< double> >::const_iterator
          stat_itr( MIN_MAX.Begin()), stat_itr_end( MIN_MAX.End());
        stat_itr != stat_itr_end; ++stat_itr
      )
      {
        min_line += FORMAT( stat_itr->GetMin());
        max_line += FORMAT( stat_itr->GetMax());
      }

      std::string analysis( min_line + "\n" + max_line);

      return analysis;
    }

  } // namespace restraint

} // namespace bcl
