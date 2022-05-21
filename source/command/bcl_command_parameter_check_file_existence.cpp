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
#include "command/bcl_command_parameter_check_file_existence.h"

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
    const util::SiPtr< const util::ObjectInterface> ParameterCheckFileExistence::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckFileExistence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new ParameterCheckFileExistence
    ParameterCheckFileExistence *ParameterCheckFileExistence::Clone() const
    {
      return new ParameterCheckFileExistence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckFileExistence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns true if PARAMETER is an allowed parameters. i.e. true if the file exists
    //! @param PARAMETER - the parameter in question. This is the filename.
    //! @param PARAMETER_NAME the name of the parameter being checked
    //! @param ERROR_STREAM - the ostream to which errors are written
    //! @return bool - returns true if the filename given by PARAMETER exists, false otherwise
    bool ParameterCheckFileExistence::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {
      // if this file was not found
      if( !io::DirectoryEntry( PARAMETER).DoesExist())
      {
        ERROR_STREAM << "File \"" << PARAMETER << "\" does not exist!" << '\n';
        return false;
      }

      // file was found
      return true;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckFileExistence::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write list of allowed parameters. i.e. any filename
      OSTREAM << "any existent file";

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckFileExistence::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckFileExistence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace command
} // namespace bcl
