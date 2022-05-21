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

#ifndef BCL_COMMAND_COMMAND_H_
#define BCL_COMMAND_COMMAND_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Command
    //! @brief TODO: document
    //!
    //! @see @link example_command_command.cpp @endlink
    //! @author woetzen
    //! @date Nov 17, 2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Command :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< FlagInterface>       m_Parameters;       //!< these are the parameters given after the executable
      util::ShPtrVector< FlagInterface> m_FlagsWithParams;  //!< these are all flags with a list of params
      util::ShPtrVector< FlagInterface> m_DefaultFlags;     //!< Generic bcl flags

      // map from flag name to corresponding flag
      storage::Map< std::string, util::SiPtr< FlagInterface> > m_FlagNamesToFlags;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! response file char
      static const char s_ResonseFileChar = '@';

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Command();

      //! @brief copy constructor
      Command( const Command &COMMAND);

      //! @brief virtual copy constructor
      //! @return pointer to Command object
      Command *Clone() const;

      //! @brief assignment operator
      Command &operator =( const Command &ORIGINAL);

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const;

      //! @brief returns the parameters
      //! @return the parameters
      const util::ShPtr< FlagInterface> &GetParameters() const;

      //! @brief sets the parameters
      //! @param PARAMETERS parameters to set
      void SetParameters( const util::ShPtr< FlagInterface> &PARAMETERS);

      //! @brief returns the flags with parameters that are application-specific (e.g. not bcl-wide)
      //! @return FlagInterfaces
      const util::ShPtrVector< FlagInterface> &GetAppFlagsWithParams() const;

      //! @brief returns the flags whose presence is independent of the application
      //! @return FlagInterfaces of all default flags that have been added to this command
      const util::ShPtrVector< FlagInterface> &GetBclFlagsWithParams() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief read from storage::Vector of strings
      //! @param ARGUMENT_LIST - the list of arguments
      //! @param ERROR_STREAM - the error stream to which errors should be written
      //! @param DRY_RUN - true is this is only a dry run, and no signals should be triggered (default false)
      //! @return bool - returns true on success and false if there was an error
      bool ReadArguments
      (
        const storage::Vector< std::string> &ARGUMENT_LIST,
        std::ostream &ERROR_STREAM,
        const bool DRY_RUN = false
      );

      //! @brief read from a command state
      //! @param STATE - the current command state
      //! @param ERROR_STREAM - the error stream to which errors should be written
      //! @return bool - returns true on success and false if there was an error
      bool SetFlags( CommandState &STATE, std::ostream &ERROR_STREAM);

      //! @brief checks whether a flag is set
      //! @param STRING the flag to check
      //! @return bool - returns true if the flag is set, false if it's not
      bool IsFlagSet( const std::string &STRING);

      //! @brief returns a flag with parameters with given NAME
      //! @param NAME - the name of the requested flag
      //! @return util::ShPtr< FlagInterface> - a pointer to the requested flag
      util::SiPtr< FlagInterface> GetFlagWithParams( const std::string &NAME);

      //! @brief add an additional parameter to the commandline
      //! @param COMMANDLINE_PARAMETER adds given ParameterInterface
      void AddParameter( const util::ShPtr< ParameterInterface> &COMMANDLINE_PARAMETER);

      //! @brief add an additional FlagWithParams to the commandline
      //! @param COMMANDLINE_FLAG_PARAMS add given FlagInterface
      void AddFlag( const util::ShPtr< FlagInterface> &COMMANDLINE_FLAG_PARAMS);

      //! @brief add additional FlagWithParams
      //! @param COMMANDLINE_FLAG_PARAMS adds all FlagInterfaces
      void PushBack( const util::ShPtrVector< FlagInterface> &COMMANDLINE_FLAG_PARAMS);

      //! @brief reset all flags and parameters
      void ResetFlagsAndParameters();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief returns the string seperating sections
      //! @return line break terminated separator string
      static const std::string &DefaultSectionSeparator();

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteHelp( std::ostream &OSTREAM) const;

      //! @brief writes the user provided commandline
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteUserCommand( std::ostream &OSTREAM) const;

      //! @brief writes the usage command, complete with required parameters and flags
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteUsage( std::ostream &OSTREAM) const;

    protected:

      //! @brief writes Command to ostream OSTREAM
      //! @param OSTREAM - the outstream to write to
      //! @param INDENT indentation
      //! @return std::ostream &OSTREAM - return the stream after writing
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read Command for std::istream ISTREAM
      //! @param ISTREAM - the instream to read from
      //! @return std::istream &ISTREAM - return the stream after reading
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief regenerate the flag to names map
      void RegenerateFlagNameMap();

    }; //class Command

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_COMMAND_H_
