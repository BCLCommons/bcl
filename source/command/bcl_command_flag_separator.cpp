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
#include "command/bcl_command_flag_separator.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FlagSeparator::s_Instance
    (
      GetObjectInstances().AddInstance( new FlagSeparator())
    );

    //! static variable to undefined list of parameters
    util::ShPtrVector< ParameterInterface> FlagSeparator::s_Parameters;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FlagSeparator::FlagSeparator() :
      m_SeparatorText()
    {
    }

    //! @brief constructor from a separator text
    FlagSeparator::FlagSeparator( const std::string &SEPARATOR_TEXT) :
      m_SeparatorText( SEPARATOR_TEXT)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FlagSeparator
    FlagSeparator *FlagSeparator::Clone() const
    {
      return new FlagSeparator( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FlagSeparator::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns name of flag
    //! @return name of flag as std::string
    const std::string &FlagSeparator::GetName() const
    {
      // initialize static empty string and return it
      static std::string s_name;
      return s_name;
    }

    //! @brief set name - only possible if there was no name given yet
    //! @param NAME - what you want to name the flag
    void FlagSeparator::SetName( const std::string &NAME)
    {

    }

    //! @brief returns description of flag
    //! @return description of flag as string
    const std::string &FlagSeparator::GetDescription() const
    {
      // initialize empty string and return it
      static const std::string s_empty;
      return s_empty;
    }

    //! @brief set flag to true if it was false
    void FlagSeparator::SetFlag()
    {

    }

    //! @brief set flag to false if it was true
    void FlagSeparator::UnsetFlag()
    {

    }

    //! @brief reset the flag
    //! @detail resets all internal parameters, then removes all internally held parameters
    void FlagSeparator::ResetFlag()
    {

    }

    //! @brief returns true if was set in command line
    //! @return whether this flag was set in the command line
    bool FlagSeparator::GetFlag() const
    {
      return true;
    }

    //! @brief returns the function to be called whenever this flag is updated
    //! @return the function to be called whenever this flag is updated
    FlagInterface::t_Signal FlagSeparator::GetSignal() const
    {
      return NULL;
    }

    //! @brief returns the number of parameters
    //! @return number of parameters after flag
    size_t FlagSeparator::GetSize() const
    {
      return size_t( 0);
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a ShPtrVector
    const util::ShPtrVector< ParameterInterface> &FlagSeparator::GetParameterList() const
    {
      return s_Parameters;
    }

    //! @brief returns m_ParameterList
    //! @return the parameter list as a util::ShPtrVector< ParameterInterface>
    util::ShPtrVector< ParameterInterface> &FlagSeparator::GetParameterList()
    {
      return s_Parameters;
    }

    //! @brief returns the first parameter
    //! @return the first parameter for this flag as a util::ShPtr< ParameterInterface>
    const util::ShPtr< ParameterInterface> &FlagSeparator::GetFirstParameter() const
    {
      // undefined first parameter
      static const util::ShPtr< ParameterInterface> s_undefined;

      // this function should never be called
      BCL_Assert( false, "this class should never have been asked for first parameter");

      // end
      return s_undefined;
    }

    //! @brief returns the # of required parameters for this flag
    //! @return the effective # of required parameters for this flag
    size_t FlagSeparator::GetNumberRequiredParameters() const
    {
      return 0;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief fill the parameter list from storage vector of strings
    //! @param PARAMETER_LIST a storage vector of strings containing the parameters as strings
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagSeparator::ReadFromList( const storage::Vector< std::string> &PARAMETER_LIST, std::ostream &ERROR_STREAM)
    {
      return true;
    }

    //! @brief checks that for every parameter was given a default value or it was passed through commandline
    //! @brief or if it is a dynamic list, if it meets the size specifications
    //! @param ERROR_STREAM the stream to which errors should be written
    //! @return true if successful, false otherwise
    bool FlagSeparator::IsValidList( std::ostream &ERROR_STREAM) const
    {
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the help was written to
    std::ostream &FlagSeparator::WriteHelp( std::ostream &OSTREAM) const
    {
      // write name and description
      OSTREAM << '\n' << m_SeparatorText << '\n';

      //end
      return OSTREAM;
    }

    //! @brief writes the user provided commandline
    //! @param OSTREAM the stream to which the commandline is written
    //! @return given ostream to which the commandline was written
    std::ostream &FlagSeparator::WriteUserCommand( std::ostream &OSTREAM) const
    {
      // end
      return OSTREAM;
    }

    //! @brief writes the usage command, complete with required parameters and flags
    //! @param OSTREAM - the outstream to write to
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &FlagSeparator::WriteRequiredUsage( std::ostream &OSTREAM) const
    {
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FlagSeparator::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SeparatorText, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FlagSeparator::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SeparatorText, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
