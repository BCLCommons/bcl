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

#ifndef BCL_COMMAND_FLAG_SEPARATOR_H_
#define BCL_COMMAND_FLAG_SEPARATOR_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_flag_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FlagSeparator
    //! @brief Helper class that allows visual separation of flags when help is asked for
    //! @details This is a dummy FlagInterface derived class to be used as a visual separation between flags when
    //! help information is printed out.
    //!
    //! @see @link example_command_flag_separator.cpp @endlink
    //! @author karakam
    //! @date Nov 23, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FlagSeparator :
      public FlagInterface
    {

    private:

    //////////
    // data //
    //////////

      //! text to be used as separator
      std::string m_SeparatorText;

      //! static variable to undefined list of parameters
      static util::ShPtrVector< ParameterInterface> s_Parameters;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FlagSeparator();

      //! @brief constructor from a separator text
      FlagSeparator( const std::string &SEPARATOR_TEXT);

      //! @brief Clone function
      //! @return pointer to new FlagSeparator
      FlagSeparator *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get separator text
      //! @return separator text
      const std::string &GetSeparatorText() const
      {
        return m_SeparatorText;
      }

      //! @brief returns name of flag
      //! @return name of flag as std::string
      const std::string &GetName() const;

      //! @brief set name - only possible if there was no name given yet
      //! @param NAME - what you want to name the flag
      void SetName( const std::string &NAME);

      //! @brief returns description of flag
      //! @return description of flag as string
      const std::string &GetDescription() const;

      //! @brief set flag to true if it was false
      void SetFlag();

      //! @brief set flag to false if it was true
      void UnsetFlag();

      //! @brief reset the flag
      //! @detail resets all internal parameters, then removes all internally held parameters
      void ResetFlag();

      //! @brief returns true if was set in command line
      //! @return whether this flag was set in the command line
      bool GetFlag() const;

      //! @brief returns the function to be called whenever this flag is updated
      //! @return the function to be called whenever this flag is updated
      t_Signal GetSignal() const;

      //! @brief returns the number of parameters
      //! @return number of parameters after flag
      size_t GetSize() const;

      //! @brief returns m_ParameterList
      //! @return the parameter list as a ShPtrVector
      const util::ShPtrVector< ParameterInterface> &GetParameterList() const;

      //! @brief returns m_ParameterList
      //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
      util::ShPtrVector< ParameterInterface> &GetParameterList();

      //! @brief returns the first parameter
      //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
      const util::ShPtr< ParameterInterface> &GetFirstParameter() const;

      //! @brief returns the # of required parameters for this flag
      //! @return the effective # of required parameters for this flag
      size_t GetNumberRequiredParameters() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief fill the parameter list from storage vector of strings
      //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
      //! @param ERROR_STREAM the stream to which errors should be written
      //! @return true if successful, false otherwise
      virtual bool ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM);

      //! @brief checks that for every parameter was given a default value or it was passed through commandline
      //! @brief or if it is a dynamic list, if it meets the size specifications
      //! @param ERROR_STREAM the stream to which errors should be written
      //! @return true if successful, false otherwise
      virtual bool IsValidList( std::ostream &ERROR_STREAM) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM the stream to which the help is written to
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM) const;

      //! @brief writes the user provided commandline
      //! @param OSTREAM the stream to which the commandline is written
      //! @return given ostream to which the commandline was written
      std::ostream &WriteUserCommand( std::ostream &OSTREAM) const;

      //! @brief writes the usage command, complete with required parameters and flags
      //! @param OSTREAM - the outstream to write to
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteRequiredUsage( std::ostream &OSTREAM) const;

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

    }; // class FlagSeparator

  } // namespace command
} // namespace bcl

#endif // BCL_COMMAND_FLAG_SEPARATOR_H_
