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

#ifndef BCL_COMMAND_PARAMETER_CHECK_FILE_EXISTENCE_H_
#define BCL_COMMAND_PARAMETER_CHECK_FILE_EXISTENCE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckFileExistence
    //! @brief This class is a command line helper class to check if files exist
    //!
    //! @see @link example_command_parameter_check_file_existence.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ParameterCheckFileExistence :
      public ParameterCheckInterface
    {

    public:

    //////////
    // data //
    //////////

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      //! @return pointer to new ParameterCheckFileExistence
      ParameterCheckFileExistence *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return const ref std::string - the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns true if PARAMETER is an allowed parameters. i.e. true if the file exists
      //! @param PARAMETER the parameter to check
      //! @param PARAMETER_NAME the name of the parameter being checked
      //! @param ERROR_STREAM the stream to which errors are written
      //! @return returns true if the parameter is allowed, false otherwise
      bool IsAllowedParameter
      (
        const std::string &PARAMETER,
        const std::string &PARAMETER_NAME,
        std::ostream &ERROR_STREAM
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the commandline
      //! @param OSTREAM - the stream to which you should write
      //! @param INDENT - amount to indent each new line after the first
      //! @return std::ostream &OSTREAM - return the stream after you write to it
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      //! @see ParameterCheckInterface::Write
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // ParameterCheckFileExistence

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_CHECK_FILE_EXISTENCE_H_
