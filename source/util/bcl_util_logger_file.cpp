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

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "util/bcl_util_cleanable_interface.h"
#include "util/bcl_util_loggers.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoggerFile
    //! @brief LoggerInterface implementation printing to a file
    //! @details implementation of the LoggerInterface printing all messages to the commandline-set file
    //!
    //! @author heinzes1
    //! @date 01/12/2009
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoggerFile :
      public LoggerInterface,
      public CleanableInterface
    {

    private:

    //////////
    // data //
    //////////

      //! file stream
      io::OFStream m_LoggerStream;

      //! message identifier
      std::string m_Filename;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoggerFile() :
        LoggerInterface( std::cout.rdbuf()),
        m_Filename()
      {
      }

      //! @brief overwritten copy constructor
      //! @param LOGGER_FILE reference to a LoggerFile object
      LoggerFile( const LoggerFile &LOGGER_FILE) :
        LoggerInterface( LOGGER_FILE.rdbuf()),
        m_LoggerStream(),
        m_Filename()
      {
        SetLogIdentifier( LOGGER_FILE.m_Filename);
      }

      //! @brief Clone function
      //! @return pointer to new LoggerDefault
      virtual LoggerFile *Clone() const
      {
        return new LoggerFile( *this);
      }

      //! @brief virtual destructor
      virtual ~LoggerFile()
      {
        CleanUp();
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief virtual function that will be called if cleanable is registered
      virtual void CleanUp()
      {
        if( m_LoggerStream.is_open())
        {
          m_LoggerStream.close();
          m_LoggerStream.clear();
        }
      }

      //! @brief get the ideal maximum line width for the given logger
      //! @return the ideal maximum line width for the given logger
      size_t GetMaxLineWidth() const
      {
        return LoggerInterface::GetDefaultMaxLineWidth();
      }

      //! @brief set an identifier string for the messages
      //! @param IDENTIFIER string that appears in front of each message
      //! @return true on success
      bool SetLogIdentifier( const std::string &IDENTIFIER)
      {
        // was already initialized
        if( m_Filename == IDENTIFIER)
        {
          return true;
        }

        // stream with different filename is open - needs to be closed
        if( m_LoggerStream.is_open())
        {
          m_LoggerStream.close();
          m_LoggerStream.clear();
        }

        // set new filename and open file
        m_Filename = IDENTIFIER;
        io::File::MustOpenOFStream( m_LoggerStream, m_Filename);

        // set the buffer to the files buffer
        std::ostream::rdbuf( m_LoggerStream.rdbuf());

        // end
        return true;
      }

      //! @brief get the identifier
      //! @return the identifier
      const std::string &GetLogIdentifier() const
      {
        return m_Filename;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief pass message to the log
      //! @param MESSAGE reference to message string
      void LogMessage( const std::string &MESSAGE)
      {
        *this << "=" << MESSAGE << std::endl;
      }

      //! @brief pass error message to the log
      //! @param ERR_STRING reference to error message string
      void LogStatus( const std::string &STATUS)
      {
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! instance of the file logger added to the enum
      static Logger s_LoggerFileInstance;

    }; // class LoggerFile

    Logger LoggerFile::s_LoggerFileInstance
    (
      GetLoggers().AddEnum( "File", ShPtr< LoggerInterface>( new LoggerFile()))
    );

  } // namespace util
} // namespace bcl
