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
#include "assemble/bcl_assemble_printer_protein_model_movie.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_multiplier.h"
#include "assemble/bcl_assemble_quality.h"
#include "assemble/bcl_assemble_quality_batch.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> PrinterProteinModelMovie::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterProteinModelMovie())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PrinterProteinModelMovie::PrinterProteinModelMovie() :
      m_Prefix(),
      m_Tag(),
      m_MoviePrinter(),
      m_ScoringFunction(),
      m_StepStatusSet(),
      m_Superimpose( quality::GetSuperimposeMeasures().e_NoSuperimpose),
      m_Qualities(),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>())
    {
    }

    //! @brief construct from data members and information required to initialize movie printer
    //! @param PREFIX prefix as absolute path
    //! @param MOVIE_PRINTER mc movie printer to use
    //! @param SCORING_FUNCTION scoring function
    //! @param STEP_STATUS_SET step statuses to print
    //! @param SUPERIMPOSE superimpose measure to use
    //! @param QUALITIES quality measures to calculate
    PrinterProteinModelMovie::PrinterProteinModelMovie
    (
      const std::string &PREFIX,
      const util::ShPtr< mc::MoviePrinterInterface> &MOVIE_PRINTER,
      const util::ShPtr< score::ProteinModelScoreSum> &SCORING_FUNCTION,
      const storage::Set< opti::StepStatusEnum> &STEP_STATUS_SET,
      const quality::SuperimposeMeasure &SUPERIMPOSE,
      const storage::Set< quality::Measure> &QUALITIES
    ) :
      m_Prefix( PREFIX),
      m_MoviePrinter( MOVIE_PRINTER),
      m_ScoringFunction( SCORING_FUNCTION),
      m_StepStatusSet( STEP_STATUS_SET),
      m_Superimpose( SUPERIMPOSE),
      m_Qualities( QUALITIES),
      m_RoundNumber( 0),
      m_StageNumber( util::GetUndefined< size_t>())
    {
    }

    //! @brief Clone function
    //! @return pointer to new PrinterProteinModelMovie
    PrinterProteinModelMovie *PrinterProteinModelMovie::Clone() const
    {
      return new PrinterProteinModelMovie( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PrinterProteinModelMovie::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return const reference to the prefix member variable
    //! @return prefix as string
    const std::string &PrinterProteinModelMovie::GetPrefix() const
    {
      return m_Prefix;
    }

    //! @brief required by interface but not implemented since the movie printer has its own prefix that
    //!        should not change
    //! @param PREFIX new prefix
    void PrinterProteinModelMovie::SetPrefix( const std::string &PREFIX)
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief reset and initialize the printer
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModelMovie::Initialize( const size_t &ROUND_NUMBER)
    {
      // call initialize with undefined stage number
      Initialize( m_RoundNumber, util::GetUndefined< size_t>());
    }

    //! @brief reset and initialize the printer with the given round and stage number
    //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
    //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
    //! @return true if initialization was successful
    void PrinterProteinModelMovie::Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
    {
      // update round number
      m_RoundNumber = ROUND_NUMBER;

      // set the stage number
      m_StageNumber = STAGE_NUMBER;
    }

    //! @brief prints information concerning the approximation process based on the status of the given tracker
    //! @param TRACKER holds the status of the approximation process
    void PrinterProteinModelMovie::Print( const opti::Tracker< ProteinModel, double> &TRACKER) const
    {
      // call WriteToFile function
      WriteToFile( TRACKER.GetBest(), CreateFileName( "final"), "final");

      // initialize output for script
      io::OFStream write;

      // construct the script filename
      std::string movie_script_filename( GetPrefix() + GetRoundNumberFormat()( m_RoundNumber));
      if( util::IsDefined( m_StageNumber))
      {
        movie_script_filename += "_" + GetStageNumberFormat()( m_StageNumber);
      }
      movie_script_filename += "_movie." + m_MoviePrinter->GetScriptFileExtension();
      io::File::MustOpenOFStream( write, movie_script_filename);

      // write the script
      m_MoviePrinter->WriteScript( write);
      io::File::CloseClearFStream( write);

      BCL_MessageCrt
      (
        "run the script " + movie_script_filename + " and use the following command line to encode the video: " +
        m_MoviePrinter->FFMpegCommandLine()
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PrinterProteinModelMovie::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &PrinterProteinModelMovie::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief static function to write given model to using the given FILENAME and TAG
    //! @param ARGUMENT_RESULT_PAIR Pair of argument and corresponding result
    //! @param FILENAME filename of the file to be written
    //! @param TAG tag to be used
    //! @return whether writing succeeded
    bool PrinterProteinModelMovie::WriteToFile
    (
      const util::ShPtr< storage::Pair< ProteinModel, double> > &ARGUMENT_RESULT_PAIR,
      const std::string &FILENAME,
      const std::string &TAG
    ) const
    {
      // open stream
      io::OFStream write;
      io::File::MustOpenOFStream( write, FILENAME);

      // copy the model
      util::ShPtr< ProteinModel> copy_ptr( ARGUMENT_RESULT_PAIR->First().HardCopy());
      ProteinModel &copy( *copy_ptr);

      // superimpose if specified
      if( m_Superimpose.IsDefined() && m_Superimpose != quality::GetSuperimposeMeasures().e_NoSuperimpose)
      {
        // superimpose the model and store the transformation to be used for the multiplier
        const math::TransformationMatrix3D transform
        (
          Quality::SuperimposeModel( m_Superimpose, copy, biol::GetAtomTypes().CA).Second()
        );

        // cast a pointer to the multiplier data if any
        util::ShPtr< ProteinModelMultiplier> sp_multiplier
        (
          copy.GetProteinModelData()->GetData( ProteinModelData::e_Multiplier)
        );

        // if the pointer is defined
        if( sp_multiplier.IsDefined())
        {
          // set the transformation
          sp_multiplier->Transform( math::Inverse( transform));

          // update the protein data
          util::ShPtr< ProteinModelData> sp_protein_model_data( copy.GetProteinModelData());
          sp_protein_model_data->Replace( ProteinModelData::e_Multiplier, sp_multiplier);
          copy.SetProteinModelData( sp_protein_model_data);
        }
      }

      // initialize header
      storage::Table< double> table
      (
        math::SumFunctionMixin< score::ProteinModel>::GetValueTableVerticalColumnNames()
      );

      // write score to file
      if( m_ScoringFunction.IsDefined())
      {
        table.Append( m_ScoringFunction->CreateSortedReadableTable( ARGUMENT_RESULT_PAIR->First()));
      }
      else
      {
        write << "score: " << ARGUMENT_RESULT_PAIR->Second() << '\n';
      }

      // if qualities were given
      if( !m_Qualities.IsEmpty())
      {
        // create quality section
        table.InsertRow
        (
          "Quality",
          storage::Vector< double>::Create
          (
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        );

        table.Append
        (
          QualityBatch( m_Qualities, biol::GetAtomTypes().CA, "", true).ConstructTable( ARGUMENT_RESULT_PAIR->First())
        );
      }

      // write model to file
      pdb::Factory().WriteModelToPDB( copy, write);
      BCL_MessageStd( "pdb written to " + util::Format()( FILENAME));

      // close and clear stream
      io::File::CloseClearFStream( write);

      // add to movie script
      if( TAG == "final")
      {
        // add frame
        m_MoviePrinter->AddFinalFrame( FILENAME, table);
      }
      else
      {
        opti::StepStatusEnum step_status;
        std::stringstream err_stream;
        if( step_status.TryRead( util::ObjectDataLabel( TAG), err_stream))
        {
          // add frame
          m_MoviePrinter->AddFrame( FILENAME, step_status, table);
        }
      }

      // return success
      return true;
    }

    //! @brief creates filename from minimization information
    //! @param TAG tag for the step status if applicable
    //! @param ITERATION number of iteration
    //! @return the full filename of a model pdb
    std::string PrinterProteinModelMovie::CreateFileName( const std::string &TAG, const size_t ITERATION) const
    {
      // initialize string
      std::string filename( m_Prefix + GetRoundNumberFormat()( m_RoundNumber) + "_");

      // if valid stage
      if( util::IsDefined( m_StageNumber))
      {
        // add stage prefix
        filename += GetStageNumberFormat()( m_StageNumber) + "_";
      }

      // if iteration is defined
      if( util::IsDefined( ITERATION))
      {
        filename += GetIterationNumberFormat()( ITERATION) + "_";
      }

      // add the tag
      filename += TAG + ".pdb";

      // end
      return filename;
    }

  } // namespace assemble
} // namespace bcl
