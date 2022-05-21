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

#ifndef BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_MOVIE_H_
#define BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_MOVIE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_movie_printer_interface.h"
#include "mc/bcl_mc_print_interface.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterProteinModelMovie
    //! @brief Prints out protein models and a script for generating a movie from those models
    //! @details Writes out protein models during a minimization.  A script is then generated to visualize the
    //!          minimization using the given MC movie printer.
    //!
    //! @see @link example_assemble_printer_protein_model_movie.cpp @endlink
    //! @author weinerbe
    //! @date Nov 23, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterProteinModelMovie :
      public mc::PrintInterface< ProteinModel, double>
    {

    private:

    //////////
    // data //
    //////////

      //! prefix
      std::string m_Prefix;

      //! tag
      std::string m_Tag;

      //! movie printer to use
      mutable util::ShPtr< mc::MoviePrinterInterface> m_MoviePrinter;

      //! scoring function for a more detailed score output
      util::ShPtr< score::ProteinModelScoreSum> m_ScoringFunction;

      //! m_StepStatusSet collection of step statuses (like accepted, improved, ...) for which info is printed
      storage::Set< opti::StepStatusEnum> m_StepStatusSet;

      //! superimposition measure to use
      quality::SuperimposeMeasure m_Superimpose;

      //! set of qualities to be calculated and printed
      storage::Set< quality::Measure> m_Qualities;

      //! round number
      size_t m_RoundNumber;

      //! stage number
      size_t m_StageNumber;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterProteinModelMovie();

      //! @brief construct from data members and information required to initialize movie printer
      //! @param PREFIX prefix as absolute path
      //! @param MOVIE_PRINTER mc movie printer to use
      //! @param SCORING_FUNCTION scoring function
      //! @param STEP_STATUS_SET step statuses to print
      //! @param SUPERIMPOSE superimpose measure to use
      //! @param QUALITIES quality measures to calculate
      PrinterProteinModelMovie
      (
        const std::string &PREFIX,
        const util::ShPtr< mc::MoviePrinterInterface> &MOVIE_PRINTER,
        const util::ShPtr< score::ProteinModelScoreSum> &SCORING_FUNCTION,
        const storage::Set< opti::StepStatusEnum> &STEP_STATUS_SET,
        const quality::SuperimposeMeasure &SUPERIMPOSE,
        const storage::Set< quality::Measure> &QUALITIES
      );

      //! @brief Clone function
      //! @return pointer to new PrinterProteinModelMovie
      PrinterProteinModelMovie *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return const reference to the prefix member variable
      //! @return prefix as string
      const std::string &GetPrefix() const;

      //! @brief required by interface but not implemented since the movie printer has its own prefix that
      //!        should not change
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX);

    ////////////////
    // operations //
    ////////////////

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER);

      //! @brief reset and initialize the printer with the given round and stage number
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER);

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< ProteinModel, double> &TRACKER) const;

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

    private:

      //! @brief static function to write given model to using the given FILENAME and TAG
      //! @param ARGUMENT_RESULT_PAIR Pair of argument and corresponding result
      //! @param FILENAME filename of the file to be written
      //! @param TAG tag to be used
      //! @return whether writing succeeded
      bool WriteToFile
      (
        const util::ShPtr< storage::Pair< ProteinModel, double> > &ARGUMENT_RESULT_PAIR,
        const std::string &FILENAME,
        const std::string &TAG
      ) const;

      //! @brief creates filename from minimization information
      //! @param TAG tag for the step status if applicable
      //! @param ITERATION number of iteration
      //! @return the full filename of a model pdb
      std::string CreateFileName( const std::string &TAG, const size_t ITERATION = util::GetUndefinedSize_t()) const;

    }; // class PrinterProteinModelMovie

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PRINTER_PROTEIN_MODEL_MOVIE_H_
