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

#ifndef BCL_ASSEMBLE_PRINTER_TRACKER_HISTORY_H_
#define BCL_ASSEMBLE_PRINTER_TRACKER_HISTORY_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "mc/bcl_mc.fwd.hh"
#include "score/bcl_score.fwd.hh"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_print_interface.h"
#include "mc/bcl_mc_temperature_interface.h"
#include "quality/bcl_quality_measures.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterTrackerHistory
    //! @brief printer to be used that allows printing a detailed tracker history for debugging
    //! @details This class is derived from PrinterProteinModelScores and allows printing iteration by iteration detailed
    //! information for debugging purposes:
    //! - Iteration #
    //! - Step Status
    //! - Scheme of the mutate applied
    //! - Number of SSEs in the model
    //! - Temperature
    //! - Individual scores
    //! - Sum score
    //! - Quality measures
    //!
    //! @see @link example_assemble_printer_tracker_history.cpp @endlink
    //! @author karakam
    //! @date Apr 10, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PrinterTrackerHistory :
      public mc::PrintInterface< ProteinModel, double>
    {

    //////////
    // data //
    //////////

    private:

      //! Path to be used
      std::string m_Path;

      //! Prefix to be used (contains no path)
      std::string m_Prefix;

      //! ShPtr to the scoring function to be used
      util::ShPtr< score::ProteinModelScoreSum> m_ScoringFunction;

      //! ShPtr to the temperature used
      util::ShPtr< mc::TemperatureInterface> m_Temperature;

      //! Qualities
      storage::Set< quality::Measure> m_Qualities;

      //! vector to store column names
      storage::Vector< std::string> m_ColumnNames;

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
      PrinterTrackerHistory();

      //! @brief constructor from a prefix, ShPtrs to scoring function, tracker, temperature and a native model
      //! @param PREFIX Prefix to be added to output filename
      //! @param SP_SCORING_FUNCTION ShPtr to scoring function to be used
      //! @param SP_TRACKER ShPtr to tracker that is used
      //! @param SP_TEMPERATURE ShPtr to temperature used
      //! @param QUALITIES set of quality measures to calculate
      PrinterTrackerHistory
      (
        const std::string &PATH,
        const std::string &PREFIX,
        const util::ShPtr< score::ProteinModelScoreSum> &SP_SCORING_FUNCTION,
        const util::ShPtr< mc::TemperatureInterface> &SP_TEMPERATURE,
        const storage::Set< quality::Measure> &QUALITIES
      );

      //! @brief Clone function
      //! @return pointer to new PrinterTrackerHistory
      PrinterTrackerHistory *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return prefix
      //! @return prefix
      const std::string &GetPrefix() const;

      //! @brief set prefix to given PREFIX
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX);

    ////////////////
    // operations //
    ////////////////

      //! @brief the filename the history is written to
      //! @return the filename the history is written to
      std::string Filename() const;

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

    public:

      //! @brief collects column names in a vector of strings and returns it
      //! @return vector of string that holds the column names
      storage::Vector< std::string> CollectColumnNames() const;

      //! @brief function to write to stream
      //! @param TRACKER holds the status of the approximation process
      //! @return whether writing succeeded
      void WriteToFile( const opti::Tracker< ProteinModel, double> &TRACKER) const;

    }; // class PrinterTrackerHistory

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PRINTER_TRACKER_HISTORY_H_
