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

#ifndef BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_SCORE_HEATMAP_H_
#define BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_SCORE_HEATMAP_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "score/bcl_score_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeAtomDistanceScoreHeatmap
    //! @brief creates heat map showing how frequently a given score is achieved by each atom distance restraint
    //! @details creates heat map showing how frequently a given score is achieved by each atom distance restraint
    //!
    //! @see @link example_restraint_analyze_atom_distance_score_heatmap.cpp @endlink
    //! @author alexanns
    //! @date Jul 31, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeAtomDistanceScoreHeatmap :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the object that will be used to score the models in the ensemble
      util::Implementation< score::RestraintAtomDistance> m_Score;

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the minimal value representing the left boundary of the score histogram
      double m_HistogramMinimum;

      //! the width of one bin of the score histograms
      double m_HistogramBinSize;

      //! the number of bins in the score histograms
      size_t m_HistogramNumberOfBins;

      //! actual restraint type
      util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > > m_RestraintType;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeAtomDistanceScoreHeatmap();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistanceScoreHeatmap
      AnalyzeAtomDistanceScoreHeatmap *Clone() const;

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

      //! @brief creates histograms of how frequently a score is achieved by each restraint
      //! @param DATA the restraints that will be scored
      //! @param ENSEMBLE the ensemble of models that will be scored according to the atom distance restraints
      //! @return vector of histograms - one histogram for each restraint + one for average sum
      storage::Vector< math::Histogram> GetScoreHistograms
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
      ) const;

    public:

      //! @brief creates the tics that will be used in the heatmap - converts each restraint into a string to use as tic
      //! @param DATA the restraints that will be used to create tics
      //! @return vector of strings which are the tics representing each restraint
      static storage::Vector< std::string> GetRestraintNameTics( const util::ShPtrVector< AtomDistance> &DATA);

      //! @brief creates gnuplot heat map object from the score histograms and the restraint names as tics
      //! @param SCORE_HISTOGRAMS the histograms that will be used to make the heat map
      //! @param TICS the names of the restraints
      //! @return heat map object which represents the distribution of scores for each restraint
      static math::GnuplotHeatmap GetHeatMap
      (
        const storage::Vector< math::Histogram> &SCORE_HISTOGRAMS, const storage::Vector< std::string> &TICS
      );

    }; // class AnalyzeAtomDistanceScoreHeatmap

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_SCORE_HEATMAP_H_
