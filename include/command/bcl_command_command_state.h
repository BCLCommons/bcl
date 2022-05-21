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

#ifndef BCL_COMMAND_COMMAND_STATE_H_
#define BCL_COMMAND_COMMAND_STATE_H_

// includes namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_default_flag_types.h"
#include "signal/bcl_signal_signal.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CommandState
    //! @brief This class is a command line helper class derived from ObjectInterface
    //!
    //! @see @link example_command_command_state.cpp @endlink
    //! @author mendenjl
    //! @date Nov 21, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CommandState :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! map of flag names (key) to given vector of strings (parameters)
      //! These are the flags and values given over the command line that have not yet been used
      storage::Map< std::string, storage::Vector< std::string> > m_AvailableFlags;

      //! whether the command state is in a "dry-run" state, and thus no signal functions should be called
      bool m_DryRun;

      //! next flag type to handle
      FlagTypeEnum m_CurrentFlagType;

      //! signal handler for ParseArguments
      mutable signal::Signal1< const CommandState &> m_ParseArgumentsSignal;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! response file char
      static const char s_ResponseFileChar = '@';

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param DRY_RUN true if no flags should have their signal functions called
      CommandState( const bool &DRY_RUN = false);

      //! @brief virtual copy constructor
      //! @return pointer to CommandState object
      CommandState *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get all flags that remain in the state
      const storage::Map< std::string, storage::Vector< std::string> > &GetState() const;

      //! @brief get access to the global command state
      //! @return the command state object that will be used in int main() const
      static CommandState &GetGlobalCommandState();

      //! @brief get access to a bool that indicates whether this is the static initialization phase
      //! @return bool that indicates hether this is the static initialization phase
      static bool &IsInStaticInitialization();

      //! @brief get access to a bool that indicates whether help was requested anywhere on the command line
      //! @return bool that indicates whether help was requested anywhere on the command line
      static bool &GetWasHelpRequested();

      //! @brief get access to a bool that indicates whether help was given anywhere on the command line
      //! @return bool that indicates whether help was given anywhere on the command line
      static bool &GetWasHelpGiven();

      //! @brief get whether main() is currently parsing the command line
      //! @return bool that indicates whether main() is currently parsing the command line
      static bool &GetInMainCommandLineParsing();

    ////////////////
    // operations //
    ////////////////

      //! @brief sets the internal argument map from the given command line strings
      //! @param NUMBER_ARGUMENTS argc from the command line
      //! @param ARGUMENTS given arguments
      //! @param ERR_STREAM stream to write errors out to
      //! @return true on success
      bool ParseArguments( const int NUMBER_ARGUMENTS, const char **ARGUMENTS, std::ostream &ERR_STREAM);

      //! @brief sets the internal argument map from the given command line strings
      //! @param CMD_LINE command line strings
      //! @param ERR_STREAM stream to write errors out to
      //! @return true on success
      bool ParseArguments( const storage::Vector< std::string> &CMD_LINE, std::ostream &ERR_STREAM);

      //! @brief access to the parse arguments signal handler
      //! @return ref to the SignalHandler that emits when ParseArguments is called
      signal::Signal1< const CommandState &> &GetParseArgumentsSignal() const
      {
        return m_ParseArgumentsSignal;
      }

      //! @brief get the number of parameters remaining
      //! @return the number parameters that remain to be parsed
      size_t GetNumberRemainingParameters() const;

      //! @brief get the arguments for a particular flag (empty vector if the flag was not set)
      //! @param FLAG_NAME the flag of interest, should not have the - prefix
      //! @return arguments for the given flag, if any were given
      const storage::Vector< std::string> &GetArguments( const std::string &FLAG_NAME) const;

      //! @brief test whether a particular flag was given
      //! @param FLAG_NAME the flag of interest, should not have the - prefix
      //! @return true if the flag was given
      bool WasFlagGiven( const std::string &FLAG_NAME) const;

      //! @brief get the next flag type to consider
      //! @return the next flag type to handle
      const FlagType &GetFlagType() const;

      //! @brief set next flag type to consider
      //! @param TYPE the next type of flag to consider
      void SetFlagType( const FlagType &TYPE);

      //! @brief Use *this to handle a particular flag, if it was set; and then remove it from the state
      //! @param FLAG the flag to handle
      //! @param ERR_STREAM the stream to use for writing errors to
      //! @return true on success
      bool Update( FlagInterface &FLAG, std::ostream &ERR_STREAM);

      //! @brief Use *this to handle a particular parameter, if it was set; and then remove it from the state
      //! @param PARAMETER the parameter to handle
      //! @param ERR_STREAM the stream to use for writing errors to
      //! @return true on success
      bool Update( ParameterInterface &PARAMETER, std::ostream &ERR_STREAM);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief writes CommandState to ostream OSTREAM
      //! @param OSTREAM - the outstream to write to
      //! @param INDENT indentation
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read CommandState for std::istream ISTREAM
      //! @param ISTREAM - the instream to read from
      //! @return std::istream &ISTREAM - return the stream after reading
      std::istream &Read( std::istream &ISTREAM);

    private:

    ////////////////
    // operations //
    ////////////////

      //! @brief handle setting the signal for all flags
      void HandleDryRunState();

      //! @brief handle unsetting the signal for all flags
      void UnsetDryRunState();

      //! @brief reads in parameters until the next flag is reached
      //! @param ARG_ITR - an iterator over the arguments
      //! @param ARG_ITR_END - the end of the iterator
      //! @return storage::Vector< std::string> - the parameters
      static storage::Vector< std::string> GetParametersTillNextFlag
      (
        std::vector< std::string>::const_iterator &ARG_ITR,
        const std::vector< std::string>::const_iterator &ARG_ITR_END
      );

      //! @brief expand response files recursively
      //! @param ARGUMENTS list of arguments
      //! @param ERR_STREAM the stream to use for writing errors to
      //! @return arguments where arguments preceding '@' are treated as files, opened and inserted as a list of arguments
      static std::pair< storage::Vector< std::string>, bool> ExpandResponseFiles
      (
        const storage::Vector< std::string> &ARGUMENTS,
        std::ostream &ERR_STREAM
      );

    }; //class CommandState

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_COMMAND_STATE_H_
