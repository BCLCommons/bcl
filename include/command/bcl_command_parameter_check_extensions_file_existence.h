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

#ifndef BCL_COMMAND_PARAMETER_CHECK_EXTENSIONS_FILE_EXISTENCE_H_
#define BCL_COMMAND_PARAMETER_CHECK_EXTENSIONS_FILE_EXISTENCE_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckExtensionsFileExistence
    //! @brief a command line helper class to check for the existence of proper file extensions
    //! @details for example, you would set the m_Extensions to be ".pdb", ".fasta", and ".ascii"
    //! then when you ask IsAllowedParameter("1ubi"), it would return true if a 1ubi.pdb, a 1ubi.fasta,
    //! and a 1ubi.ascii exist.  If all of these files do not exist, it would return false.
    //!
    //! @see @link example_command_parameter_check_extensions_file_existence.cpp @endlink
    //! @author heinzes1, woetzen, karakam
    //! @date 10/12/2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ParameterCheckExtensionsFileExistence :
      public ParameterCheckInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::Vector< std::string> m_Extensions; //!< the extensions of the file

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ParameterCheckExtensionsFileExistence();

      //! @brief construct from allowed parameters
      //! @param EXTENSIONS vector of file extensions
      ParameterCheckExtensionsFileExistence( const storage::Vector< std::string> &EXTENSIONS);

      //! @brief construct from allowed parameters
      //! @param EXTENSION single file extension
      explicit ParameterCheckExtensionsFileExistence( const std::string &EXTENSION);

      //! @brief virtual copy constructor
      //! @return pointer to ParameterCheckExtensionsFileExistence object
      ParameterCheckExtensionsFileExistence *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return const ref std::string - the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns if PARAMETER is an allowed, i.e. if PARAMETER has an allowed extension and the file exists
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
      //! @param OSTREAM output stream to be written to
      //! @param INDENT the amount to indent each new line after the first
      //! @return ostream which was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const;

      //! @brief writes the list of allowed values
      //! @param &OSTREAM - the stream to which you write
      //! @return std::ostream &OSTREAM - return the stream after you write to it
      std::ostream &WriteListOfExtensions( std::ostream &OSTREAM) const;

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

    }; //class ParameterCheckExtensionsFilesExistence

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_CHECK_EXTENSIONS_FILE_EXISTENCE_H_
