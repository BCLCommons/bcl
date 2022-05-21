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

#ifndef BCL_COMMAND_PARAMETER_CHECK_FILE_IN_SEARCH_PATH_H_
#define BCL_COMMAND_PARAMETER_CHECK_FILE_IN_SEARCH_PATH_H_

// include the namespace header
#include "bcl_command.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_command_parameter_check_interface.h"
#include "io/bcl_io_directory.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ParameterCheckFileInSearchPath
    //! @brief Parameter check for whether files/directories exist with a search path that may be relative to the bcl exectuble
    //!        If the user provides a path of their own, the search path is ignored
    //!
    //! @see @link example_command_parameter_check_file_in_search_path.cpp @endlink
    //! @author mendenjl
    //! @date Jan 06, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ParameterCheckFileInSearchPath :
      public ParameterCheckInterface
    {

    private:

    //////////
    // data //
    //////////

      storage::Vector< std::string> m_Paths; //!< Paths to search;

    public:

    //////////
    // data //
    //////////

      //! @brief get a string variable that, when used with this class, will be replaced by the directory containing the bcl executable
      //! @return a string whose value will be replaced with the directory containing the bcl executable
      static const std::string &GetExecutableDirectoryVariable();

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param TARGET_NAME target file or directory name to search for in the default search path
      //! use the path the bcl was started from or up to 3 directories above it, before resorting to
      //! the directory on blue, or the current directory.  This allows both installed bcl executables to
      //! be able to find the models directory without relying on the often relative and often moved install path.
      //! It also allows modification and use of files in checkouts without needing to set this flag
      //! @param LAST_RESORT last file or directory to try before giving up; this is usually a meilerlab-specific path;
      //! The last resort filename allows MeilerLab users to copy locally-built executables around on the filesystem
      //! without needing to copy all associated files, unless they have in fact been modified
      //! @param TYPE the actual type of object to search for; be it a file or directory
      ParameterCheckFileInSearchPath
      (
        const std::string &TARGET_NAME = "",
        const std::string &LAST_RESORT = "",
        const io::Directory::EntryType &TYPE = io::Directory::e_File
      );

      //! @brief constructor
      //! @param TARGET_NAME target file or directory name to search for in the given paths
      //! @param PATHS ordered paths to search for the target file or directory
      //! @param TYPE the actual type of object to search for; be it a file or directory
      ParameterCheckFileInSearchPath
      (
        const std::string &TARGET_NAME,
        const storage::Vector< std::string> &PATHS,
        const io::Directory::EntryType &TYPE = io::Directory::e_File
      );

      //! @brief constructor
      //! @param PATHS ordered paths to search for a valid entry of the given type
      //! @param TYPE the actual type of object to search for; be it a file or directory
      ParameterCheckFileInSearchPath
      (
        const storage::Vector< std::string> &PATHS,
        const io::Directory::EntryType &TYPE = io::Directory::e_File
      );

      //! @brief virtual copy constructor
      //! @return pointer to new ParameterCheckFileInSearchPath
      ParameterCheckFileInSearchPath *Clone() const;

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

      //! @brief find the file in the search paths; if it exists, otherwise, return an empty string
      //! @param PARAMETER the first path to check; may be a default, blank, or may be something entered over the command line
      //! @return the first occurrance of the target file/dir in the search paths; empty if it is not in the search path
      std::string FindFile( const std::string &PARAMETER) const;

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

      //! @brief clean the parameter value; that is, make it usable by code that is independent of the parameter check
      //! Some checks may include or-type functionality, such as checking for a file among a list of paths
      //! The clean function updates the parameter value with the correct path
      //! @param PARAMETER the parameter to clean
      void Clean( std::string &PARAMETER) const;

      //! @brief return the type of directory entry being searched for
      io::Directory::EntryType GetSearchType() const;

      //! @brief get a string representation of all paths that will be searched
      std::string GetCandidatePathsString() const;

      //! @brief helper function; update all strings in m_Paths to be directories (end with /) or not
      void SetDirectories();

      //! @brief helper function; update all strings in m_Paths to be directories (end with /) or not
      void SetFiles();

      //! @brief function to simplify a given path by removing the pattern X/../
      static std::string SimplifyPath( const std::string &PATH);

    }; // ParameterCheckFileInSearchPath

  } // namespace command
} // namespace bcl

#endif //BCL_COMMAND_PARAMETER_CHECK_FILE_IN_SEARCH_PATH_H_
