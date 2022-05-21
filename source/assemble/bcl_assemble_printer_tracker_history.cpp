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
#include "assemble/bcl_assemble_printer_tracker_history.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_quality_batch.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_temperature_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterTrackerHistory::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterTrackerHistory())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default construactor
    PrinterTrackerHistory::PrinterTrackerHistory() :
      m_Path(),
      m_Prefix(),
      m_ScoringFunction(),
      m_Temperature(),
      m_ColumnNames()
    {
    }

    //! @brief constructor from a prefix, ShPtrs to scoring function, tracker, temperature and a native model
    //! @param PATH path to be used for writing output files
    //! @param PREFIX Prefix to be added to output filename
    //! @param SP_SCORING_FUNCTION ShPtr to scoring function to be used
    //! @param SP_TRACKER ShPtr to tracker that is used
    //! @param SP_TEMPERATURE ShPtr to temperature used
    //! @param QUALITIES set of quality measures to calculate
    PrinterTrackerHistory::PrinterTrackerHistory
    (
      const std::string &PATH,
      const std::string &PREFIX,
      const util::ShPtr< score::ProteinModelScoreSum> &SP_SCORING_FUNCTION,
      const util::ShPtr< mc::TemperatureInterface> &SP_TEMPERATURE,
      const storage::Set< quality::Measure> &QUALITIES
    ) :
      m_Path( PATH),
      m_Prefix( PREFIX),
      m_ScoringFunction( SP_SCORING_FUNCTION),
      m_Temperature( SP_TEMPERATURE),
      m_Qualities( QUALITIES),
      m_ColumnNames( CollectColumnNames())
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterTrackerHistory
    PrinterTrackerHistory *PrinterTrackerHistory::Clone() const
    {
      return new PrinterTrackerHistory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &PrinterTrackerHistory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return prefix
    //! @return prefix
    const std::string &PrinterTrackerHistory::GetPrefix() const
    {
      return m_Prefix;
    }

    //! @brief set prefix to given PREFIX
    //! @param PREFIX new prefix
    void PrinterTrackerHistory::SetPrefix( const std::string &PREFIX)
    {
      m_Prefix = PREFIX;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the filename the history is written to
    //! @return the filename the history is written to
    std::string PrinterTrackerHistory::Filename() const
    {
      const std::string filename( m_Path + util::GetRuntimeEnvironment().GetPathSeperator() + m_Prefix
        + GetRoundNumberFormat()( m_RoundNumber) + "_" + GetStageNumberFormat()( m_StageNumber) + ".tracker_history");

      return filename;
    }

    //! @brief reset and initialize the printer
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @return true if initialization was successful
    void PrinterTrackerHistory::Initialize( const size_t &ROUND_NUMBER)
    {
      // set the round number
      m_RoundNumber = ROUND_NUMBER;
    }

    //! @brief reset and initialize the printer with the given round and stage number
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
    //! @return true if initialization was successful
    void PrinterTrackerHistory::Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
    {
      // set the round number and stage number
      m_RoundNumber = ROUND_NUMBER;
      m_StageNumber = STAGE_NUMBER;
    }

    //! @brief prints information concerning the approximation process based on the status of the given tracker
    //! @param TRACKER holds the status of the approximation process
    void PrinterTrackerHistory::Print( const opti::Tracker< ProteinModel, double> &TRACKER) const
    {
      // Write history to file
      WriteToFile( TRACKER);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterTrackerHistory::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterTrackerHistory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief collects column names in a vector of strings and returns it
    //! @return vector of string that holds the column names
    storage::Vector< std::string> PrinterTrackerHistory::CollectColumnNames() const
    {
      // initialize column names
      storage::Vector< std::string> column_names
      (
        storage::Vector< std::string>::Create
        (
          "iter_no", "step_status", "nr_aas", "nr_sses", "nr_helix", "nr_strand", "temperature", "mutate"
        )
      );

      // now get the scheme names from the scoring function and append them to this vector
      column_names.Append( m_ScoringFunction->GetFunctionSchemes());

      // if the native model is defined with coordinates
      if( !m_Qualities.IsEmpty())
      {
        column_names.Append( QualityBatch::ColumnNamesFromQualities( m_Qualities));
      }
      // end
      return column_names;
    }

    //! @brief function to write to stream
    //! @param TRACKER holds the status of the approximation process
    //! @return whether writing succeeded
    void PrinterTrackerHistory::WriteToFile( const opti::Tracker< ProteinModel, double> &TRACKER) const
    {
      // static column seperator
      static const std::string col_sep( " ");

      // initialize formats
      static const util::Format column_format( util::Format().W( 12).R());
      static const util::Format number_format( util::Format().W( 12).FFP( 3).R());
      static const util::Format mutate_format( util::Format().W( 40).R());

      // initialize filename
      const std::string filename( Filename());

      // get the scheme for the last step
      const std::string &mutate_scheme( TRACKER.GetStepScheme());

      // get a reference on the last argument
      const ProteinModel &this_model( TRACKER.GetCurrent()->First());

      // initialize write
      io::OFStream write;

      // if this is the first iteration
      if( TRACKER.GetIteration() == 1)
      {
        // open with write mode
        io::File::MustOpenOFStream( write, filename);

        // write out the header for the table
        // iterate over the column names
        for
        (
          storage::Vector< std::string>::const_iterator
            name_itr( m_ColumnNames.Begin()), name_itr_end( m_ColumnNames.End());
          name_itr != name_itr_end; ++name_itr
        )
        {
          // if this is mutate column
          if( *name_itr == "mutate")
          {
            write << mutate_format( *name_itr) << col_sep;
          }
          else
          {
            write << column_format( *name_itr) << col_sep;
          }
        }
        write << '\n';
      }
      else
      {
        // otherwise open with append mode
        io::File::MustOpenOFStream( write, filename, std::ios::app);
      }

      // get the step status name
      const std::string step_status_name( TRACKER.GetStatusOfLastStep().GetString());

      // now start printing the values
      // iteration number
      write << column_format( TRACKER.GetIteration()) << col_sep;
      // step status
      write << column_format( step_status_name) << col_sep;
      // number amino acids
      write << column_format( this_model.GetNumberAAs()) << col_sep;
      // number SSEs
      write << column_format( this_model.GetNumberSSEs()) << col_sep;
      // number helices
      write << column_format( this_model.GetNumberSSE( biol::GetSSTypes().HELIX)) << col_sep;
      // number strands
      write << column_format( this_model.GetNumberSSE( biol::GetSSTypes().STRAND)) << col_sep;
      // temperature
      write << number_format( m_Temperature->GetLastCalculatedTemperature()) << col_sep;
      // write the mutate name
      write << mutate_format( mutate_scheme) << col_sep;

      // now get the scores
      storage::Table< double> scores( m_ScoringFunction->CreateValueTableHorizontal( this_model));
      // get the row with weighted values and write with the same format
      storage::Vector< double> score_vector( scores[ "weighted_value"].GetData());

      // iterate over columns
      for
      (
        storage::Vector< double>::const_iterator score_itr( score_vector.Begin()), score_itr_end( score_vector.End());
        score_itr != score_itr_end; ++score_itr
      )
      {
        write << number_format( *score_itr) << col_sep;
      }

      // if any quality measures are given
      if( !m_Qualities.IsEmpty())
      {
        // now calculate the measures
        const storage::List< storage::Pair< std::string, double> > measures_list
        (
          QualityBatch( m_Qualities, biol::GetAtomTypes().CA)( this_model)
        );

        // iterate over the measures vector
        for
        (
          storage::List< storage::Pair< std::string, double> >::const_iterator
            measure_itr( measures_list.Begin()), measure_itr_end( measures_list.End());
          measure_itr != measure_itr_end; ++measure_itr
        )
        {
          // print the value
          write << number_format( measure_itr->Second()) << col_sep;
        }
      }

      // add the newline character and return true
      write << '\n';
      io::File::CloseClearFStream( write);
    }

  } // namespace assemble
} // namespace bcl
