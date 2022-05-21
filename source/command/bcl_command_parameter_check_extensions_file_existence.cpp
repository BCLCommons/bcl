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
#include "command/bcl_command_parameter_check_extensions_file_existence.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckExtensionsFileExistence::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckExtensionsFileExistence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckExtensionsFileExistence::ParameterCheckExtensionsFileExistence()
    {
    }

    //! @brief construct from allowed parameters
    //! @param EXTENSIONS vector of file extensions
    ParameterCheckExtensionsFileExistence::ParameterCheckExtensionsFileExistence( const storage::Vector< std::string> &EXTENSIONS) :
      m_Extensions( EXTENSIONS)
    {
    }

    //! @brief construct from allowed parameters
    //! @param EXTENSION single file extension
    ParameterCheckExtensionsFileExistence::ParameterCheckExtensionsFileExistence( const std::string &EXTENSION) :
      m_Extensions()
    {
      m_Extensions.PushBack( EXTENSION);
    }

    //! @brief virtual copy constructor
    //! @return pointer to ParameterCheckExtensionsFileExistence object
    ParameterCheckExtensionsFileExistence *ParameterCheckExtensionsFileExistence::Clone() const
    {
      return new ParameterCheckExtensionsFileExistence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckExtensionsFileExistence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is an allowed (i.e. it exists & has the correct extension
    //! @param PARAMETER - the parameter in question
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM - write errors here
    //! @return true if the parameter is allowed, false otherwise
    bool ParameterCheckExtensionsFileExistence::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // iterate over given extensions to check if those files exists
      for
      (
        std::vector< std::string>::const_iterator itr( m_Extensions.Begin()), itr_end( m_Extensions.End());
        itr != itr_end; ++itr
      )
      {
        const std::string name( PARAMETER + *itr);

        // if this file was not found
        if( !io::DirectoryEntry( name).DoesExist())
        {
          ERROR_STREAM << "File \"" << name << "\" does not exist!" << '\n';
          return false;
        } // if
      } // for

      return true;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM output stream to be written to
    //! @param INDENT the amount to indent each new line after the first
    //! @return ostream which was written to
    std::ostream &ParameterCheckExtensionsFileExistence::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write list of extensions
      OSTREAM << "file pattern <";
      WriteListOfExtensions( OSTREAM);
      OSTREAM << ">";

      // end
      return OSTREAM;
    }

    //! @brief writes the list of allowed values
    //! @param &OSTREAM - the stream to which you write
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckExtensionsFileExistence::WriteListOfExtensions( std::ostream &OSTREAM) const
    {
      // write opening parantheses
      OSTREAM << "\"";

      // write all entries in the list
      for
      (
        std::vector< std::string>::const_iterator itr( m_Extensions.Begin()), itr_end( m_Extensions.End());
        itr != itr_end;
        ++itr
      )
      {
        OSTREAM << "*" << *itr;
        if( itr != m_Extensions.End() - 1)
        {
          OSTREAM << ", ";
        }
      }
      // write closing parentheses
      OSTREAM << "\"";

      // end
      return OSTREAM;
    } // WriteListOfExtensions

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckExtensionsFileExistence::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Extensions, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckExtensionsFileExistence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Extensions, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
