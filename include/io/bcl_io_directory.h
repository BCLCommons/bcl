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

#ifndef BCL_IO_DIRECTORY_H_
#define BCL_IO_DIRECTORY_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Directory
    //! @brief class that handles a directory
    //!
    //! @see @link example_io_directory.cpp @endlink
    //! @author woetzen
    //! @date Aug 2, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Directory :
      public util::ObjectInterface
    {

    public:

    ///////////
    // types //
    ///////////

      //! @enum EntryType
      //! @brief type of directory entries
      enum EntryType
      {
        e_File,    //!< file type
        e_Dir,     //!< directory type
        e_Unknown, //!< if type is unknown
        s_MaxTypes //!< number of types
      };

    private:

    //////////
    // data //
    //////////

      //! path of directory
      std::string m_Path;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Directory();

      //! @brief construct from path
      //! @param PATH path for directory
      Directory( const std::string &PATH);

      //! @brief Clone function
      //! @return pointer to new Directory
      Directory *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get path
      //! @return the path for this directory
      const std::string &GetPath() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add directory to filename
      //! @param NAME the name of the file
      //! @return full name including this path
      std::string AppendFilename( const std::string &NAME) const;

      //! @brief check if directory physically exist
      //! @return true, if is does exist, false if it does not exist as a directory (e.g. as file or something else)
      bool DoesExist() const;

      //! @brief list content of directory
      //! @return list of Directory entries
      storage::List< DirectoryEntry> ListEntries() const;

      //! @brief list entries of the directory that match the given type and suffix
      //! @param TYPE type of entry (file/directory, etc., s_MaxTypes for any)
      //! @param PREFIX if given, only list entries that beging with the given string
      //! @param SUFFIX if given, only list entries with the given suffix
      //! @return list of Directory entries
      storage::List< DirectoryEntry> ListEntries
      (
        const EntryType &TYPE,
        const std::string &PREFIX = std::string(),
        const std::string &SUFFIX = std::string()
      ) const;

      //! @brief list files of the directory that match the given type and suffix
      //! @param TYPE type of entry (file/directory, etc., s_MaxTypes for any)
      //! @param PREFIX if given, only list entries that beging with the given string
      //! @param SUFFIX if given, only list entries with the given suffix
      //! @param RECURSIVE true to auto-descend into sub-directories
      //! @note when using RECURSIVE, entries will be given in breadth-first order
      //! @return list of Directory entries
      storage::List< DirectoryEntry> ListFiles
      (
        const std::string &PREFIX = std::string(),
        const std::string &SUFFIX = std::string(),
        const bool &RECURSIVE = false
      ) const;

      //! @brief create the directory
      //! @return true if it was successful, false if it did exist or could not be created
      bool Make();

      //! @brief remove the directory
      //! @param RECURSIVE default false
      //! @return true, if removal was successful, false if directory did not exist or removal failed (if it was not empty and recursive was not passed)
      bool Remove( const bool RECURSIVE = false);

      //! @brief clear the directory, without deleting it
      //! @param RECURSIVE default false
      //! @return true, if removal of content was successful, false if directory did not exist or removal failed (if sub directories were not empty and recursive was not passed)
      bool Clear( const bool RECURSIVE = false);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief create a directory
      //! @param DIRECTORY name of directory to be created
      //! @return the Directory associated with the created directory, check if is exists to confirm successful creation
      static Directory MkDir( const std::string &DIRECTORY);

      //! @brief clean directory
      //! changes empty to '.' and removes additional '/' or '\\'
      void CleanPath();

    }; // class Directory

  } // namespace io
} // namespace bcl

#endif // BCL_IO_DIRECTORY_H_
