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
#include "command/bcl_command_parameter_check_extension.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckExtension::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckExtension())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ParameterCheckExtension::ParameterCheckExtension()
    {
    }

    //! @brief construct from allowed parameters
    //! @param EXTENSION file extension
    ParameterCheckExtension::ParameterCheckExtension( const std::string &EXTENSION) :
      m_Extension( EXTENSION)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to ParameterCheckExtension object
    ParameterCheckExtension *ParameterCheckExtension::Clone() const
    {
      return new ParameterCheckExtension( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckExtension::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns if PARAMETER is allowed, i.e. if filename has allowed extension
    //! @param PARAMETER the parameter to check
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM the stream to which errors are written
    //! @return returns true if the parameter is allowed, false otherwise
    bool ParameterCheckExtension::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // determine the position of the Extension in the given PARAMETER
      const size_t pos( PARAMETER.rfind( m_Extension, PARAMETER.length()));

      // if extension was found and is at the end of PARAMETER return true
      if( util::IsDefined( pos) && ( pos + m_Extension.length() == PARAMETER.length()))
      {
        return true;
      }

      // otherwise it was not found or is in the wrong position, return false
      ERROR_STREAM << "Given parameter \"" << PARAMETER << "\" does not have the extension: " << m_Extension << '\n';
      return false;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM output stream to be written to
    //! @param INDENT the amount to indent each new line after the first
    //! @return ostream which was written to
    std::ostream &ParameterCheckExtension::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write extension
      OSTREAM << "file pattern <*" << m_Extension << "> ";

      //end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckExtension::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Extension, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckExtension::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Extension, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
