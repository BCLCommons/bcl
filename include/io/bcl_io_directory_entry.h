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

#ifndef BCL_IO_DIRECTORY_ENTRY_H_
#define BCL_IO_DIRECTORY_ENTRY_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_directory.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DirectoryEntry
    //! @brief this class mimics the POSIX struct dirent
    //! @details The POSIX dirent struct not consistently implemented (and does not even exist for MS VS)
    //! beside the name of the directory entry, it will contain the type (FILE, DIRECTORY, LINK ...) and other things
    //! possibly later, like permissions, owner, time of modification ...
    //!
    //! @see @link example_io_directory_entry.cpp @endlink
    //! @author woetzen
    //! @date Aug 2, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DirectoryEntry :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      std::string             m_Name;      //!< name of entry
      Directory::EntryType    m_Type;      //!< type of entry
      util::ShPtr< Directory> m_Directory; //!< directory this entry is contained in

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
      DirectoryEntry();

      //! @brief construct from name
      //! @param NAME name of file, path or link
      DirectoryEntry( const std::string &NAME);

      //! @brief construct from directory and name
      //! @param SP_DIR directory the file is contained
      //! @param NAME name of file, path or link
      DirectoryEntry( const util::ShPtr< Directory> &SP_DIR, const std::string &NAME);

      //! @brief construct from directory, name, and type
      //! @param SP_DIR directory the file is contained
      //! @param NAME name of file, path or link
      //! @param TYPE type of the directory entry
      DirectoryEntry( const util::ShPtr< Directory> &SP_DIR, const std::string &NAME, const Directory::EntryType &TYPE);

      //! @brief Clone function
      //! @return pointer to new DirectoryEntry
      DirectoryEntry *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to name
      //! @return name of entry
      const std::string &GetName() const;

      //! @brief get Directory
      const Directory &GetDirectory() const;

      //! @brief fullname including the path the entry is in
      //! @return fullpath and name
      std::string GetFullName() const;

      //! @brief access to type
      //! @return type of directory entry
      Directory::EntryType GetType() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief check if entry is of given type
      //! @param TYPE the desired entry type
      //! @return true if the type matches
      bool IsType( const Directory::EntryType &TYPE) const;

      //! @brief check if the entry does exist
      //! @return true, if the entry exists, false otherwise
      bool DoesExist() const;

      //! @brief rename the directory entry
      //! @param NEW the new entry
      //! @param ALLOW_OVERWRITE whether to allow overwriting NEW if it already exists
      //! @return true, if successful
      bool Rename( const DirectoryEntry &NEW, const bool &ALLOW_OVERWRITE = false);

      //! @brief remove the directory entry
      //! @return true, if successful
      bool Remove();

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

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief set attributes
      //! sets all the entry attributes - can be necessary to be called, if the file in the file system has changed
      //! @return true if successful, false if there was an error (e.g file did not exist or other problems)
      bool SetAttributes();

    }; // class DirectoryEntry

  } // namespace io
} // namespace bcl

#endif // BCL_IO_DIRECTORY_ENTRY_H_
