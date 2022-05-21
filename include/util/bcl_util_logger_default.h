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

#ifndef BCL_UTIL_LOGGER_DEFAULT_H_
#define BCL_UTIL_LOGGER_DEFAULT_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_logger_interface.h"
#include "bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoggerDefault
    //! @brief default LoggerInterface printing to the screen
    //! @details This class is the default implementation of the LoggerInterface printing all messages to std::out and
    //! error messages to std::err
    //!
    //! @see @link example_util_logger_default.cpp @endlink
    //! @author heinzes1
    //! @date Jan 12, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoggerDefault :
      public LoggerInterface
    {

    private:

    //////////
    // data //
    //////////

      //! message identifier, printed in front of the message
      std::string m_Identifier;

      //! size of last status update; or 0 if the last message was not a status
      size_t      m_LastStatusSize;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoggerDefault() :
        LoggerInterface( std::cout.rdbuf()),
        m_Identifier(),
        m_LastStatusSize( 0)
      {
      }

      //! @brief overwritten copy constructor
      //! @param LOGGER_DEFAULT reference to a LoggerDefault object
      LoggerDefault( const LoggerDefault &LOGGER_DEFAULT) :
        LoggerInterface( LOGGER_DEFAULT.rdbuf()),
        m_Identifier( LOGGER_DEFAULT.m_Identifier),
        m_LastStatusSize( LOGGER_DEFAULT.m_LastStatusSize)
      {
      }

      //! @brief Clone function
      //! @return pointer to new LoggerDefault
      LoggerDefault *Clone() const
      {
        return new LoggerDefault( *this);
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

      //! @brief set an identifier string for the messages
      //! @param IDENTIFIER string that appears in front of each message
      //! @return true on success
      bool SetLogIdentifier( const std::string &IDENTIFIER)
      {
        // set the identifier
        m_Identifier = IDENTIFIER;

        // end
        return true;
      }

      //! @brief get the identifier
      //! @return the identifier
      const std::string &GetLogIdentifier() const
      {
        return m_Identifier;
      }

      //! @brief get the ideal maximum line width for the given logger
      //! @return the ideal maximum line width for the given logger
      size_t GetMaxLineWidth() const
      {
        return LoggerInterface::GetTerminalLineWidth();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief pass message to the log
      //! @param MESSAGE reference to message string
      void LogMessage( const std::string &MESSAGE)
      {
        if( m_LastStatusSize) // if the last message was a status message, print a new line
        {
          *this << '\n';
          m_LastStatusSize = 0;
        }
        *this << m_Identifier << "=" << MESSAGE << std::endl;
      }

      //! @brief pass error message to the log
      //! @param ERR_STRING reference to error message string
      void LogError( const std::string &ERR_STRING)
      {
        if( m_LastStatusSize) // if the last message was a status message, print a new line
        {
          std::cerr << '\n';
          m_LastStatusSize = 0;
        }
        std::cerr << m_Identifier << "=" << ERR_STRING << std::endl;
      }

      //! @brief pass error message to the log
      //! @param ERR_STRING reference to error message string
      void LogStatus( const std::string &STATUS)
      {
        if( m_LastStatusSize)
        {
          // return to the beginning of the line if the last message was a status
          *this << '\r';
        }
        *this << m_Identifier << " Status: " << STATUS;

        // set the width of the next operation to overwrite the last status message, if it was longer
        if( STATUS.size() < m_LastStatusSize)
        {
          width( m_LastStatusSize - STATUS.size());
          *this << ' ';
        }
        *this << std::flush;
        m_LastStatusSize = STATUS.size() + 1;
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
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class LoggerDefault

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_LOGGER_DEFAULT_H_
