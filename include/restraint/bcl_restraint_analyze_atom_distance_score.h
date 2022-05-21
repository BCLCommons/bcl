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

#ifndef BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_SCORE_H_
#define BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_SCORE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "io/bcl_io.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "math/bcl_math_running_min_max.h"
#include "score/bcl_score_restraint_atom_distance.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeAtomDistanceScore
    //! @brief calculates atom distance scores for ensemble of models and gives mean and sd of score of each restraint
    //! @details Can also optionally write the individual restraint score for every model in the ensemble
    //!
    //! @see @link example_restraint_analyze_atom_distance_score.cpp @endlink
    //! @author alexanns
    //! @date Jul 31, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeAtomDistanceScore :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the object that will be used to score the models in the ensemble
      util::Implementation< score::RestraintAtomDistance> m_Score;

      //! true if the restraint scores for every model should be printed
      bool m_PrintAllModelScores;

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! true if the mean, std dev, min, and max statistics should be printed for each atom distance restraint
      bool m_PrintStatistics;

      //! type of distance restraints used
      util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > > m_RestraintType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeAtomDistanceScore();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistanceScore
      AnalyzeAtomDistanceScore *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief gives the string that is the scoring analysis
      //! @param DATA the list of atom distance objects which are the restraints
      //! @PARAM ENSEMBLE the ensemble from which distance mean and std devs will be calculated
      //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
      //! @param FORMAT the format object to format the scores
      std::string GetScoreAnalysis
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
        const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
      ) const;

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
      > GetScoreStatistics
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE,
        const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT, std::string &ANALYSIS_STRING
      ) const;

      //! @brief gives formatted analysis of the statistics of the atom distance restraints
      //! @param MEAN_SD the mean and standard deviation statistics of the atom distance restraints
      //! @param MIN_MAX the min and max statistics of the atom distance restraints
      //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
      //! @param FORMAT the format object to format the scores
      //! @return string which has the formatted statistics analysis of atom restraint scores
      std::string GetScoreStatisticsFormattedAnalysis
      (
        const storage::Vector< math::RunningAverageSD< double> > &MEAN_SD,
        const storage::Vector< math::RunningMinMax< double> > &MIN_MAX,
        const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
      ) const;

      //! @brief gives formatted analysis of the mean and sd statistics of the atom distance restraints
      //! @param MEAN_SD the mean and standard deviation statistics of the atom distance restraints
      //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
      //! @param FORMAT the format object to format the scores
      //! @return string which has the formatted mean and sd statistics analysis of atom restraint scores
      static std::string GetMeanSDLines
      (
        const storage::Vector< math::RunningAverageSD< double> > &MEAN_SD,
        const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
      );

      //! @brief gives formatted analysis of the min and max statistics of the atom distance restraints
      //! @param MIN_MAX the min and max statistics of the atom distance restraints
      //! @param LINE_NAME_FORMAT the format to use to format the name of the line i.e. the first column
      //! @param FORMAT the format object to format the scores
      //! @return string which has the formatted min and max statistics analysis of atom restraint scores
      static std::string GetMinMaxLines
      (
        const storage::Vector< math::RunningMinMax< double> > &MIN_MAX,
        const util::Format &LINE_NAME_FORMAT, const util::Format &FORMAT
      );

    }; // class AnalyzeAtomDistanceScore

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_SCORE_H_
