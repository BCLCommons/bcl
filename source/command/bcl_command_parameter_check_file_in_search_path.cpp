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
#include "command/bcl_command_parameter_check_file_in_search_path.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "app/bcl_app_apps.h"
#include "io/bcl_io_directory_entry.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace command
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ParameterCheckFileInSearchPath::s_Instance
    (
      GetObjectInstances().AddInstance( new ParameterCheckFileInSearchPath( storage::Vector< std::string>()))
    );

    //! @brief get a string variable that, when used with this class, will be replaced by the directory containing the bcl executable
    //! @return a string whose value will be replaced with the directory containing the bcl executable
    const std::string &ParameterCheckFileInSearchPath::GetExecutableDirectoryVariable()
    {
      static const std::string s_bcl_var( "$BCL_EXE$");
      return s_bcl_var;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new ParameterCheckFileInSearchPath
    ParameterCheckFileInSearchPath *ParameterCheckFileInSearchPath::Clone() const
    {
      return new ParameterCheckFileInSearchPath( *this);
    }

  /////////////////
  // data access //
  /////////////////

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
    ParameterCheckFileInSearchPath::ParameterCheckFileInSearchPath
    (
      const std::string &TARGET_NAME,
      const std::string &LAST_RESORT,
      const io::Directory::EntryType &TYPE
    ) :
      m_Paths
      (
        GetVersion().IsLicense()
        // release/licensed copies are installed, so everything is usually together with the executable
        // but use the putative install directory as a backup -- its rarely correct because users can move the bcl
        // folder post-installation; moreover its just whatever the user types as the install directory, which is usually
        // a relative path.
        // The absolute last resort in all cases is the current directory
        ? storage::Vector< std::string>::Create
          (
            GetExecutableDirectoryVariable(),
            LAST_RESORT,
            ""
          )
        // for locally-built (unlicensed) versions, look in the super-directories of the bcl-exectuable 1st
        // because the build-folder often has the same folders as the working directory, so e.g. we don't want to get
        // build/linux64_release/histogram, rather, we want build/linux64_release/../../histogram == histogram
        // indeed, in this case, the folder the bcl is located in should be the last one searched due to the likelyhood
        // of finding similarly named folders in the build directory
        : storage::Vector< std::string>::Create
          (
            GetExecutableDirectoryVariable() + "/../..",
            GetExecutableDirectoryVariable() + "/../../..",
            GetExecutableDirectoryVariable() + "/..",
            GetExecutableDirectoryVariable() + "/bcl", // needed for ligand-design gui
            LAST_RESORT,
            "",
            GetExecutableDirectoryVariable()
          )
      )
    {
      SetDirectories();
      for( storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End()); itr != itr_end; ++itr)
      {
        if( *itr != LAST_RESORT)
        {
          *itr += TARGET_NAME;
        }
      }
      TYPE == io::Directory::e_File ? SetFiles() : SetDirectories();
    }

    //! @brief constructor
    //! @param TARGET_NAME target file or directory name to search for in the given paths
    //! @param PATHS ordered paths to search for the target file or directory
    //! @param TYPE the actual type of object to search for; be it a file or directory
    ParameterCheckFileInSearchPath::ParameterCheckFileInSearchPath
    (
      const std::string &TARGET_NAME,
      const storage::Vector< std::string> &PATHS,
      const io::Directory::EntryType &TYPE
    ) :
      m_Paths( PATHS)
    {
      if( !TARGET_NAME.empty())
      {
        // force all existing path components to be directories
        SetDirectories();

        // add or remove / depending on the desired entry type TYPE
        std::string to_append( TARGET_NAME);
        if( TARGET_NAME[ TARGET_NAME.size() - 1] == '/')
        {
          // remove / if type was supposed to be file
          if( TYPE == io::Directory::e_File)
          {
            to_append.erase( to_append.size() - 1);
          }
        }
        else if( TYPE == io::Directory::e_Dir)
        {
          // append / if it was supposed to be a directory
          to_append += '/';
        }

        // add the target name to each path
        for
        (
          storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
          itr != itr_end;
          ++itr
        )
        {
          *itr += to_append;
        }
      }
      else
      {
        GetSearchType() == io::Directory::e_File ? SetFiles() : SetDirectories();
      }
    }

    //! @brief constructor
    //! @param PATHS ordered paths to search for a valid entry of the given type
    //! @param TYPE the actual type of object to search for; be it a file or directory
    ParameterCheckFileInSearchPath::ParameterCheckFileInSearchPath
    (
      const storage::Vector< std::string> &PATHS,
      const io::Directory::EntryType &TYPE
    ) :
      m_Paths( PATHS)
    {
      TYPE == io::Directory::e_File ? SetFiles() : SetDirectories();
    }

    //! @brief returns class name
    //! @return const ref std::string - the class name
    const std::string &ParameterCheckFileInSearchPath::GetClassIdentifier() const
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
    bool ParameterCheckFileInSearchPath::IsAllowedParameter
    (
      const std::string &PARAMETER,
      const std::string &PARAMETER_NAME,
      std::ostream &ERROR_STREAM
    ) const
    {

      // this function is only called when the user has provided a path in which case the default search path is ignored
      io::DirectoryEntry entry( PARAMETER);
      if( entry.DoesExist() && entry.GetType() == GetSearchType())
      {
        return true;
      }
      ERROR_STREAM << PARAMETER << " is not a valid "
                   << ( GetSearchType() == io::Directory::e_Dir ? "directory" : "file")
                   << "for parameter " << PARAMETER_NAME;
      return false;
    } // IsAllowedParameter

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief writes the help for the commandline
    //! @param OSTREAM - the stream to which you should write
    //! @param INDENT - amount to indent each new line after the first
    //! @return std::ostream &OSTREAM - return the stream after you write to it
    std::ostream &ParameterCheckFileInSearchPath::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // determine default
      const std::string default_path( FindFile( ""));
      // write list of allowed parameters. i.e. any filename
      OSTREAM << "any " << ( GetSearchType() == io::Directory::e_Dir ? "directory" : "file")
              << ", if not provided, defaults to " << default_path
              << "; search path is: {" << GetCandidatePathsString() << "}\n";

      // end
      return OSTREAM;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ParameterCheckFileInSearchPath::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Paths, ISTREAM);
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    //! @see ParameterCheckInterface::Write
    std::ostream &ParameterCheckFileInSearchPath::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Paths, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief find the file in the search paths; if it exists, otherwise, return an empty string
    //! @param PARAMETER the first path to check; may be a default, blank, or may be something entered over the command line
    //! @return the first occurrance of the target file/dir in the search paths; empty if it is not in the search path
    std::string ParameterCheckFileInSearchPath::FindFile( const std::string &PARAMETER) const
    {
      if( PARAMETER.size())
      {
        io::DirectoryEntry entry( PARAMETER);
        if( entry.DoesExist() && entry.GetType() == GetSearchType())
        {
          return PARAMETER;
        }
      }

      // create a string replacer to replace the special string with the bcl executable path
      const util::StringReplacement replacer
      (
        util::StringReplacement::e_Any,
        GetExecutableDirectoryVariable(),
        io::DirectoryEntry( app::Apps::GetExecutablePath()).GetDirectory().GetPath()
      );

      // get the type of entry being searched for
      io::Directory::EntryType entry_type( GetSearchType());
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string candidate_path( *itr);
        replacer.ReplaceEachIn( candidate_path);
        io::DirectoryEntry candidate( candidate_path);
        if( candidate.DoesExist() && candidate.GetType() == entry_type)
        {
          return SimplifyPath( candidate.GetFullName() + ( entry_type == io::Directory::e_File ? "" : "/"));
        }
      }

      return std::string();
    }

    //! @brief clean the parameter value; that is, make it usable by code that is independent of the parameter check
    //! Some checks may include or-type functionality, such as checking for a file among a list of paths
    //! The clean function updates the parameter value with the correct path
    //! @param PARAMETER the parameter to clean
    void ParameterCheckFileInSearchPath::Clean( std::string &PARAMETER) const
    {
      std::string file( FindFile( PARAMETER));
      if( !file.empty())
      {
        PARAMETER = file;
      }
    }

    //! @brief return the type of directory entry being searched for
    io::Directory::EntryType ParameterCheckFileInSearchPath::GetSearchType() const
    {
      return
        !m_Paths.IsEmpty() && m_Paths.FirstElement()[ m_Paths.FirstElement().size() - 1] == '/'
        ? io::Directory::e_Dir
        : io::Directory::e_File;
    }

    //! @brief get a string representation of all paths that will be searched
    std::string ParameterCheckFileInSearchPath::GetCandidatePathsString() const
    {
      // create a string replacer to replace the special string with the bcl executable path
      const util::StringReplacement replacer
      (
        util::StringReplacement::e_Any,
        GetExecutableDirectoryVariable(),
        io::DirectoryEntry( app::Apps::GetExecutablePath()).GetDirectory().GetPath()
      );

      std::string candidates;
      for
      (
        storage::Vector< std::string>::const_iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string candidate_path( *itr);
        replacer.ReplaceEachIn( candidate_path);
        candidates += SimplifyPath( candidate_path);
        if( itr + 1 != itr_end)
        {
          // add colon for next possibility
          candidates += ":";
        }
      }
      return candidates;
    }

    //! @brief helper function; update all strings in m_Paths to be directories (end with /) or not
    void ParameterCheckFileInSearchPath::SetDirectories()
    {
      // force all paths in m_Paths to have a '/' at the end
      for
      (
        storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string &str( *itr);
        if( !str.empty() && str[ str.size() - 1] != '/')
        {
          str += '/';
        }
      }
    }

    //! @brief helper function; update all strings in m_Paths to be directories (end with /) or not
    void ParameterCheckFileInSearchPath::SetFiles()
    {
      // force all paths in m_Paths to not have a '/' at the end
      for
      (
        storage::Vector< std::string>::iterator itr( m_Paths.Begin()), itr_end( m_Paths.End());
        itr != itr_end;
        ++itr
      )
      {
        std::string &str( *itr);
        while( !str.empty() && str[ str.size() - 1] == '/')
        {
          str.erase( str.size() - 1);
        }
      }
    }

    //! @brief function to simplify a given path by removing the pattern X/../
    std::string ParameterCheckFileInSearchPath::SimplifyPath( const std::string &PATH)
    {
      // remove duplicate slashes, also ./ as a directory or in the middle of the path
      util::StringReplacement replace_double_slash( util::StringReplacement::e_Any, "//", "/");
      util::StringReplacement replace_slash_dot_slash( util::StringReplacement::e_Any, "/./", "/");
      std::string temp_string( PATH);
      replace_double_slash.ReplaceAllIn( temp_string);
      replace_slash_dot_slash.ReplaceAllIn( temp_string);
      while( util::StartsWith( temp_string, "./"))
      {
        temp_string.erase( 0, 2);
      }

      if( temp_string.size() < size_t( 4))
      {
        return temp_string;
      }

      // walk through the string, find anytime /../ appears
      size_t last_found_pos( 1);
      size_t prev_directory_pos( temp_string.find( "/../", last_found_pos));
      while( prev_directory_pos < temp_string.size())
      {
        size_t start_pos( temp_string.rfind( '/', prev_directory_pos - 1));
        if( start_pos == std::string::npos)
        {
          start_pos = 0;
        }
        else
        {
          ++start_pos;
        }
        // remove the unnecessary part of the string
        temp_string.erase( start_pos, prev_directory_pos + 4 - start_pos);
        last_found_pos = start_pos ? start_pos - 1 : start_pos;
        prev_directory_pos = temp_string.find( "/../", last_found_pos);
      }
      return temp_string;
    }

  } // namespace command
} // namespace bcl
