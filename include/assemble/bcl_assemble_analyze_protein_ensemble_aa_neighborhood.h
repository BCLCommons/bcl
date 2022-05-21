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

#ifndef BCL_ASSEMBLE_ANALYZE_PROTEIN_ENSEMBLE_AA_NEIGHBORHOOD_H_
#define BCL_ASSEMBLE_ANALYZE_PROTEIN_ENSEMBLE_AA_NEIGHBORHOOD_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_analyze_protein_ensemble_interface.h"
#include "math/bcl_math_range.h"
#include "score/bcl_score_aa_neighborhood_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeProteinEnsembleAANeighborhood
    //! @brief generate a pymol file which colors residues according to their aa neighborhood score
    //!
    //! @see @link example_assemble_analyze_protein_ensemble_aa_neighborhood.cpp @endlink
    //! @author woetzen
    //! @date Sep 20, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeProteinEnsembleAANeighborhood :
      public AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! score range for coloring
      math::Range< double> m_ScoreRange;

      //! the name of the score that is used
      std::string m_ScoreName;

      //! environment score, to score an amino acid within a protein model
      util::ShPtr< score::AANeighborhoodInterface> m_Score;

      //! @brief names of allowed scores
      static const size_t s_NumberScores = 5;
      static const std::string s_ScoreNames[ s_NumberScores];

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeProteinEnsembleAANeighborhood();

      //! @brief Clone function
      //! @return pointer to new AnalyzeProteinEnsembleAANeighborhood
      AnalyzeProteinEnsembleAANeighborhood *Clone() const;

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

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const ProteinEnsemble &ENSEMBLE) const;

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

      //! @brief return a score from the score name
      //! @param SCORE_NAME
      //! @return ShPtr to exposure score
      static util::ShPtr< score::AANeighborhoodInterface> ScoreFromScoreName( const std::string &SCORE_NAME);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class AnalyzeProteinEnsembleAANeighborhood

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_ANALYZE_PROTEIN_ENSEMBLE_AA_NEIGHBORHOOD_H_
