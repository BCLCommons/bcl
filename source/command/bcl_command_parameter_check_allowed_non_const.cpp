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
#include "command/bcl_command_parameter_check_allowed_non_const.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace command
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from allowed parameters
    //! @param ALLOWED_PARAMETERS vector of allowed strings
    ParameterCheckAllowedNonConst::ParameterCheckAllowedNonConst( const storage::Vector< std::string> &ALLOWED_PARAMETERS) :
      m_AllowedParameters( &ALLOWED_PARAMETERS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckAllowedNonConst object
    ParameterCheckAllowedNonConst *ParameterCheckAllowedNonConst::Clone() const
    {
      return new ParameterCheckAllowedNonConst( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckAllowedNonConst::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. is in the list of allowed strings
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckAllowedNonConst::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // iterate over all parameters in m_AllowedParameters
      if( m_AllowedParameters->Find( PARAMETER) < m_AllowedParameters->GetSize())
      {
        return true;
      }

      // try fuzzy matching
      Guesser::GetDefaultGuesser().WriteGuesses( PARAMETER, *m_AllowedParameters, ERROR_STREAM, PARAMETER_NAME);

      // return false because parameter was not found in list
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckAllowedNonConst::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write a short description
      OSTREAM << "allowed values: ";

      // write all entries in the list
      WriteList( OSTREAM);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckAllowedNonConst::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AllowedParameters, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckAllowedNonConst::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AllowedParameters, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param OSTREAM - the stream to which you should write
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckAllowedNonConst::WriteList( std::ostream &OSTREAM) const
    {
      //write opening parentheses
      OSTREAM << "{";

      // copy to output stream
      if( !m_AllowedParameters->IsEmpty())
      {
        std::copy
        (
          m_AllowedParameters->Begin(),
          m_AllowedParameters->End() - 1,
          std::ostream_iterator< std::string>( OSTREAM, ", ")
        );
        OSTREAM << m_AllowedParameters->LastElement();
      }

      // write closing parentheses
      OSTREAM << "}";

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
