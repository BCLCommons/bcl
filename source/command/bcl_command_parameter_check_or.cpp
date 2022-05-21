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
#include "command/bcl_command_parameter_check_or.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_validation_result.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckOr::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckOr())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckOr::ParameterCheckOr()
    {
    }

    //! @brief construct from two checks
    //! @param CHECK_A the first parameter check that will be considered
    //! @param CHECK_B the second parameter check that will be considered
    ParameterCheckOr::ParameterCheckOr
    (
      const ParameterCheckInterface &CHECK_A,
      const ParameterCheckInterface &CHECK_B
    ) :
      m_Choices
      (
        util::ShPtrVector< ParameterCheckInterface>::Create
        (
          util::CloneToShPtr( CHECK_A),
          util::CloneToShPtr( CHECK_B)
        )
      )
    {
    }

    //! @brief construct from three checks
    //! @param CHECK_A the first parameter check that will be considered
    //! @param CHECK_B the second parameter check that will be considered
    //! @param CHECK_C the third parameter check that will be considered
    ParameterCheckOr::ParameterCheckOr
    (
      const ParameterCheckInterface &CHECK_A,
      const ParameterCheckInterface &CHECK_B,
      const ParameterCheckInterface &CHECK_C
    ) :
      m_Choices
      (
        util::ShPtrVector< ParameterCheckInterface>::Create
        (
          util::CloneToShPtr( CHECK_A),
          util::CloneToShPtr( CHECK_B),
          util::CloneToShPtr( CHECK_C)
        )
      )
    {
    }

    //! @brief construct from four checks
    //! @param CHECK_A the first parameter check that will be considered
    //! @param CHECK_B the second parameter check that will be considered
    //! @param CHECK_C the third parameter check that will be considered
    //! @param CHECK_D the fourth parameter check that will be considered
    ParameterCheckOr::ParameterCheckOr
    (
      const ParameterCheckInterface &CHECK_A,
      const ParameterCheckInterface &CHECK_B,
      const ParameterCheckInterface &CHECK_C,
      const ParameterCheckInterface &CHECK_D
    ) :
      m_Choices
      (
        util::ShPtrVector< ParameterCheckInterface>::Create
        (
          util::CloneToShPtr( CHECK_A),
          util::CloneToShPtr( CHECK_B),
          util::CloneToShPtr( CHECK_C),
          util::CloneToShPtr( CHECK_D)
        )
      )
    {
    }

    //! @brief construct from vector of allowed checks
    //! @param CHECKS vector of parameter checks, any of which can be satisfied for IsAllowedParameter to return true
    ParameterCheckOr::ParameterCheckOr( const util::ShPtrVector< ParameterCheckInterface> &CHECKS) :
      m_Choices( CHECKS)
    {
    }

    //! @brief copy constructor
    //! @return pointer to new ParameterCheckOr object
    ParameterCheckOr *ParameterCheckOr::Clone() const
    {
      return new ParameterCheckOr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckOr::GetClassIdentifier() const
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
    bool ParameterCheckOr::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // create a temporary error stream; whose results will be written to ERROR_STREAM only if none of
      // the parameter checks is satisfied
      std::ostringstream temp_err_stream;

      // iterate over all possible choices
      for
      (
        util::ShPtrVector< ParameterCheckInterface>::const_iterator
          itr( m_Choices.Begin()), itr_end( m_Choices.End());
        itr != itr_end;
        ++itr
      )
      {
        std::ostringstream stream;

        // get a reference to the parameter check
        const ParameterCheckInterface &parameter_check( **itr);

        // if this parameter check was valid, return true; errors for other checks will be disregarded
        if( parameter_check.IsAllowedParameter( PARAMETER, PARAMETER_NAME, stream))
        {
          return true;
        }

        temp_err_stream << stream.str();
      }

      ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" is none of the following:\n";
      WriteList( ERROR_STREAM, 1);
      ERROR_STREAM << "Errors for each check: " << temp_err_stream.str() << '\n';

      // return false because parameter was not found in list
      return false;
    }

    //! @brief check the parameter, returning CheckResult
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns CheckResult
    io::ValidationResult ParameterCheckOr::Check
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // create a temporary error stream; whose results will be written to ERROR_STREAM only if none of
      // the parameter checks is satisfied
      std::ostringstream temp_err_stream;

      io::ValidationResult overall_result( io::e_Allowed);

      // iterate over all possible choices
      for
      (
        util::ShPtrVector< ParameterCheckInterface>::const_iterator
          itr( m_Choices.Begin()), itr_end( m_Choices.End());
        itr != itr_end;
        ++itr
      )
      {
        std::ostringstream stream;

        // get a reference to the parameter check
        const ParameterCheckInterface &parameter_check( **itr);

        // get the check result
        io::ValidationResult result( parameter_check.Check( PARAMETER, PARAMETER_NAME, stream));

        // handle trivial result
        if( result.IsAllowed())
        {
          return result;
        }

        // handle help request
        if( result.IsHelp())
        {
          if( overall_result.IsInvalid())
          {
            // clear the current temp error stream
            temp_err_stream.str( std::string());
          }
          overall_result = result;
        }
        else if( overall_result.IsAllowed()) // first invalid parameter
        {
          overall_result = result;
        }

        if( result == overall_result)
        {
          temp_err_stream << stream.str();
        }
      }

      if( overall_result.IsInvalid())
      {
        ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" is none of the following:\n";
        WriteList( ERROR_STREAM, 1);
        ERROR_STREAM << "Errors for each check: " << temp_err_stream.str() << '\n';
      }
      else // help
      {
        ERROR_STREAM << temp_err_stream.str() << '\n';
      }
      // return false because parameter was not found in list
      return overall_result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    //! @see WriteList
    std::ostream &ParameterCheckOr::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write a short description
      OSTREAM << "Any of the following:\n";

      // write all entries in the list
      WriteList( OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    //! @see ParameterCheckInterface::Read
    std::istream &ParameterCheckOr::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Choices, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckOr::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Choices, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT number of indentations
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckOr::WriteList( std::ostream &OSTREAM, const size_t &INDENT) const
    {
      // no choices, nothing to print
      if( m_Choices.IsEmpty())
      {
        return OSTREAM;
      }

      std::ostringstream temp_ostream;

      // write help for the first parameter check
      m_Choices.FirstElement()->WriteHelp( temp_ostream, INDENT);

      // get the string from the output stream
      std::string temp_string( temp_ostream.str());
      if( !temp_string.empty() && temp_string[ temp_string.length() - 1] != '\n')
      {
        temp_string += '\n';
      }
      OSTREAM << temp_string;

      // iterate over remaining options
      for
      (
        util::ShPtrVector< ParameterCheckInterface>::const_iterator
          itr( m_Choices.Begin() + 1), itr_end( m_Choices.End());
        itr != itr_end;
        ++itr
      )
      {
        // write Or to separate the next set of options
        io::Serialize::InsertIndent( OSTREAM, INDENT) << "Or ";

        // get a reference to the parameter check
        const ParameterCheckInterface &parameter_check( **itr);

        temp_ostream.str( std::string());
        // write help for this parameter check
        parameter_check.WriteHelp( temp_ostream, INDENT);
        temp_string = temp_ostream.str();
        if( !temp_string.empty() && temp_string[ temp_string.length() - 1] != '\n')
        {
          temp_string += '\n';
        }

        OSTREAM << temp_string;
      }

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
