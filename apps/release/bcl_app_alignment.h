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

#ifndef BCL_APP_ALIGNMENT_H_
#define BCL_APP_ALIGNMENT_H_

// include header of this class
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "function/bcl_function.fwd.hh"
#include "sspred/bcl_sspred.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "score/bcl_score_aa_assignments.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    // TODO: align 1EPW.fasta 1HN0.fasta -open_gap -0.5 -extend_gap -0.15 -open_boundary_gap -.2 -extend_boundary_gap -0.1 -pam100 -blosum45
    // above command line will not work because the blast profile is not read since blast is not scored, but blast profile is needed

    // TODO: jufo is unnecessarily calculated
    // align 1EPW.fasta 1HN0.fasta -open_gap -0.5 -extend_gap -0.15 -open_boundary_gap -.2 -extend_boundary_gap -0.1 -pam100 -blosum45

    // align.exe  Z:\dunbrack\1ABA_  Z:\dunbrack\1BEA_  Z:\dunbrack\1GSA_  Z:\dunbrack\1MLA_  Z:\dunbrack\1SRA_  Z:\dunbrack\1TCA_  Z:\dunbrack\2BAA_  Z:\dunbrack\3CLA_  Z:\dunbrack\1ARB_  Z:\dunbrack\1CFB_  Z:\dunbrack\1SVB_  Z:\dunbrack\2A0B_  Z:\dunbrack\3SEB_  Z:\dunbrack\3VUB_  -identity -blast -sam -jufo -psipred

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Alignment
    //! @brief TODO: add brief comment
    //!
    //! @see @link example_app_alignment.cpp @endlink
    //! @author heinzes1
    //! @date
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Alignment :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! compute parameters to convert PairScores to ZScores (set this option to true and call it with just a single pairscore of interest and a list of represenative sequences (~20)
      static const bool s_ZScoreParameter = false;

      // number gap penalties
      static const size_t s_NumberGapPenalties = 4;

      // option strings gap penalties
      static const std::string s_GapPenaltiesDescription[];

      // default gap penalties
      static const std::string s_DefaultGapPenalties[];

      //! input ascii and fasta files
      util::ShPtr< command::FlagInterface> m_FastaListFlag;

      //! prefix to the output alignment file
      util::ShPtr< command::FlagInterface> m_OutputprefixFlag;

      // variables for weights for pair scores: commandlineflagwithparams
      storage::Map< score::AAAssignment, util::ShPtr< command::FlagInterface> > m_PairScoreFlags;

      // variables for gap scores: commandlineflagwithparams
      util::ShPtrVector< command::FlagInterface> m_GapPenaltyFlags;

      //! aligner to calculate the alignments
      util::ShPtr< command::FlagInterface> m_FlagAligner;

      //! input alignment for scoring only
      util::ShPtr< command::FlagInterface> m_ScoreAlignmentFlag;

    public:

      // instantiate enumerator for this class
      static const ApplicationType Alignment_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      Alignment();

    public:

      //! @brief Clone function
      //! @return pointer to new Alignment
      Alignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the bcl::commons name
      //! @return string for the bcl::commons name of that application
      std::string GetBCLScopedName() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! Main
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! read sequence profile data
      biol::AASequence ReadAllSequenceData
      (
        const std::string &CODE,
        const storage::Set< score::AAAssignment> &READ,
        std::ostream &OUTSTREAM = util::GetLogger()
      ) const;

      //! AddScore
      bool AddScore
      (
        const score::AAAssignment &SCORE,
        const double SCORE_WEIGHT,
        function::BinarySum< const biol::AABase, const biol::AABase, double> &SCORE_FCT,
        std::ostream &OUTSTREAM = util::GetLogger()
      ) const;

      //! estimate parameters to convert pair scores to Z-scores
      void ComputeParametersToConvertPairScoresToZScores
      (
        const storage::Vector< biol::AASequence> &SEQUENCES,
        const storage::Set< score::AAAssignment> &SELECTED_SCORES,
        std::ostream &OUTSTREAM = util::GetLogger()
      ) const;

      //! read secondary structure prediction
      biol::AASequence &ReadSSPrediction
      (
        biol::AASequence &SEQUENCE,
        const std::string &CODE,
        const sspred::Method &METHOD
      ) const;

    }; // class Alignment

  } // namespace app
} // namespace bcl
#endif // BCL_APP_ALIGNMENT_H_
