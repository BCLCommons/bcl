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

#ifndef BCL_MC_PRINTER_FILE_H_
#define BCL_MC_PRINTER_FILE_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_print_interface.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "opti/bcl_opti_step_status.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrinterFile
    //! @brief default printer to be used Monte Carlo minimization classes
    //! @details This class provides an implementation for opti::PrintInterface
    //!
    //! @see @link example_mc_printer_file.cpp @endlink
    //! @author karakam, fischea
    //! @date Feb 23, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class PrinterFile :
      public PrintInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! round number
      size_t m_RoundNumber;

      //! stage number
      size_t m_StageNumber;

      //! brief directory
      io::Directory m_Path;

      //! prefix for files
      std::string m_Prefix;

      //! m_StepStatusSet collection of step statuses (like accepted, improved, ...) for which info is printed
      storage::Set< opti::StepStatusEnum> m_StepStatusSet;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PrinterFile() :
        m_RoundNumber( 0),
        m_StageNumber( 0),
        m_Path(),
        m_Prefix()
      {
      }

      //! @brief construct with all member variables
      //! @param PREFIX prefix string
      //! @param STEP_STATUS_SET collection of step statuses (like accepted, improved, ...) for which info is printed
      //! @param TAG tag string
      PrinterFile
      (
        const std::string &PREFIX,
        const storage::Set< opti::StepStatusEnum> &STEP_STATUS_SET
      ) :
        m_RoundNumber( 0),
        m_StageNumber( util::GetUndefined< size_t>()),
        m_Path(),
        m_Prefix(),
        m_StepStatusSet( STEP_STATUS_SET)
      {
        SetPrefix( PREFIX);
      }

      //! @brief Clone function
      //! @return pointer to new PrinterFile< t_ArgumentType, t_ResultType>
      PrinterFile< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new PrinterFile< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return prefix
      //! @return prefix
      const std::string &GetPrefix() const
      {
        return m_Prefix;
      }

      //! @brief set prefix to given PREFIX
      //! @param PREFIX new prefix
      void SetPrefix( const std::string &PREFIX)
      {
        const storage::VectorND< 2, std::string> path_prefix( io::File::SplitToPathAndFileName( PREFIX));
        m_Path = io::Directory( path_prefix.First());
        BCL_Assert( m_Path.DoesExist(), "path component of prefix does not exists!");
        m_Prefix = path_prefix.Second();
      }

      //! @brief return const reference to the step status set
      //! @return const reference to the step status set
      const storage::Set< opti::StepStatusEnum> &GetStepStatusSet() const
      {
        return m_StepStatusSet;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER)
      {
        m_RoundNumber = ROUND_NUMBER;
      }

      //! @brief reset and initialize the printer
      //! @param ROUND_NUMBER for multiple optimizations, a different round number will be passed
      //! @param STAGE_NUMBER for multiple optimizations, a different stage number will be passed
      //! @return true if initialization was successful
      void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER)
      {
        m_RoundNumber = ROUND_NUMBER;
        m_StageNumber = STAGE_NUMBER;
      }

      //! @brief prints information concerning the approximation process based on the status of the given tracker
      //! @param TRACKER holds the status of the approximation process
      void Print( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER) const
      {
        // get the phase dependent tag
        const std::string tag
        (
          TRACKER.GetPhase() == opti::e_Start ? "start" : TRACKER.GetPhase() == opti::e_End ? "final" : ""
        );
        WriteToFile( TRACKER, CreateFilename( m_Prefix, tag));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_RoundNumber  , ISTREAM);
        io::Serialize::Read( m_StageNumber  , ISTREAM);
        io::Serialize::Read( m_Path         , ISTREAM);
        io::Serialize::Read( m_Prefix       , ISTREAM);
        io::Serialize::Read( m_StepStatusSet, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_RoundNumber  , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_StageNumber  , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Path         , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Prefix       , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_StepStatusSet, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief static function to write given model to using the given FILENAME and TAG
      //! @param TRACKER tracker to be printed
      //! @param FILENAME filename of the file to be written
      void WriteToFile( const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER, const std::string &FILENAME) const
      {
        // open stream, write, close
        io::OFStream write;
        io::File::MustOpenOFStream( write, FILENAME);
        write << TRACKER.GetBest()->First();
        io::File::CloseClearFStream( write);
        BCL_MessageStd( "argument written to: " + util::Format()( FILENAME));
        BCL_MessageStd( "result: " + util::Format()( TRACKER.GetBest()->Second()));
      }

      //! @brief static function to write to stream with given prefix
      //! @param PREFIX prefix
      //! @param TAG tag for the step status if applicable
      //! @return whether writing succeeded
      std::string CreateFilename( const std::string &PREFIX, const std::string &TAG) const
      {
        return CreateFilename( PREFIX, TAG, m_RoundNumber, m_StageNumber);
      }

      //! @brief static function to write to stream with given prefix
      //! @param PREFIX prefix
      //! @param TAG tag for the step status if applicable
      //! @param ROUND_NUMBER the current minimization round number
      //! @param STAGE_NUMBER the current stage number we are at
      //! @return whether writing succeeded
      std::string CreateFilename
      (
        const std::string &PREFIX,
        const std::string &TAG,
        const size_t ROUND_NUMBER,
        const size_t STAGE_NUMBER
      ) const
      {
        // construct filename
        std::string filename( PREFIX + PrintInterface< t_ArgumentType, t_ResultType>::GetRoundNumberFormat()( ROUND_NUMBER) + "_");

        // if valid stage
        if( util::IsDefined( STAGE_NUMBER))
        {
          // add stage prefix
          filename += PrintInterface< t_ArgumentType, t_ResultType>::GetStageNumberFormat()( STAGE_NUMBER) + "_";
        }

        // add the tag
        filename += TAG + ".bcl";

        // end
        return m_Path.AppendFilename( filename);
      }

    }; // template class PrinterFile< t_ArgumentType, t_ResultType>

    // instantiate s_Instance
    template< typename t_ArgumentType, typename t_ResultType>
    const util::SiPtr< const util::ObjectInterface> PrinterFile< t_ArgumentType, t_ResultType>::s_Instance
    (
      GetObjectInstances().AddInstance( new PrinterFile< t_ArgumentType, t_ResultType>())
    );

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_PRINTER_FILE_H_
