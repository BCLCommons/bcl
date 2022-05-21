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

#ifndef BCL_UTIL_LOGGER_INTERFACE_H_
#define BCL_UTIL_LOGGER_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_object_interface.h"

// external includes - sorted alphabetically
#include <iostream>

namespace bcl
{
  namespace util
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoggerInterface
    //! @brief interface for logging standard and error messages
    //! @details The interface implements std::ostream, so messages send to loggers via operator<< are going to standard
    //! output stream
    //!
    //! @remarks example unnecessary
    //! @author heinzes1
    //! @date Jan 12, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoggerInterface :
      public ObjectInterface,
      public std::ostream
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoggerInterface() :
        std::ostream( std::cout.rdbuf())
      {
      }

      //! @brief constructor using a streambuf
      //! @param STREAM_BUFFER pointer to a output stream buffer
      LoggerInterface( std::streambuf *STREAM_BUFFER) :
        std::ostream( STREAM_BUFFER)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Logger
      virtual LoggerInterface *Clone() const = 0;

      //! @brief set an identifier string for the messages
      //! @param IDENTIFIER string that can be used to initialize the logger
      //! @return bool - true on success
      virtual bool SetLogIdentifier( const std::string &IDENTIFIER) = 0;

      //! @brief get the identifier
      //! @return the identifier
      virtual const std::string &GetLogIdentifier() const = 0;

      //! @brief get the ideal maximum line width for the given logger
      //! @return the ideal maximum line width for the given logger
      virtual size_t GetMaxLineWidth() const = 0;

      //! @brief get the default maximum line width to be used for a particular logger
      //! this value should be used whenever the output is not to a terminal or the terminal window size is unavailable
      static size_t GetDefaultMaxLineWidth();

      //! @brief get the line width of the terminal
      //! this value should be used whenever the output is to a terminal
      static size_t GetTerminalLineWidth();

      //! @brief set a flag to ignore the terminal's line width when reporting the default max line width
      //! @param IGNORE_TERMINAL_LINE_WIDTH true to ignore the terminal's line width when reporting GetTerminalLineWidth
      static void SetIgnoreTerminalLineWidth( const bool &IGNORE_TERMINAL_LINE_WIDTH);

    ////////////////
    // operations //
    ////////////////

      //! @brief log message to standard output, abstract in interface
      //! @param MESSAGE reference of the message string
      virtual void LogMessage( const std::string &MESSAGE) = 0;

      //! @brief log a status message to standard output, abstract in interface
      //! @param STATUS reference of the status string
      //! @note consecutive status messages are overwritten, unlike conventional messages
      virtual void LogStatus( const std::string &STATUS)
      {
        // write to standard out as default
        LogMessage( STATUS);
      }

      //! @brief log message to error output
      //! @param ERR_STRING reference to the error string
      virtual void LogError( const std::string &ERR_STRING)
      {
        // write to standard output as default
        LogMessage( ERR_STRING);
      }

    private:

      //! @brief get whether the terminal's max line width should be ignored by all loggers.  This is useful when
      //!        running the examples, because we don't want output files to change based on terminal settings
      static bool &GetIgnoreTerminalLineWidth()
      {
        static bool s_ignore_terminal_line_width( false);
        return s_ignore_terminal_line_width;
      }

    }; // class Logger

    // forward declaration of access function to the current global logger
    BCL_API LoggerInterface &GetLogger();

  } // namespace util
} // namespace bcl

//! BCL Debug Macro
//! This macro prints off its argument (the actual code) and then the result of evaluating that code
//! It should never appear in committed code (outside examples that illustrate its use)
#define \
  BCL_Debug( CODE) \
  bcl::util::GetLogger() << "DEBUG: " << #CODE << " YIELDED: " << (CODE) << std::endl;

#endif // BCL_UTIL_LOGGER_INTERFACE_H_
