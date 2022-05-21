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
#include "io/bcl_io_directory_entry.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically
#include <cstdio>
#include <sys/stat.h>

namespace bcl
{
  namespace io
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DirectoryEntry::s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DirectoryEntry::DirectoryEntry() :
      m_Name( ""),
      m_Type( Directory::e_Unknown),
      m_Directory( new Directory( "."))
    {
    }

    //! @brief construct from name
    //! @param NAME name of file, path or link
    DirectoryEntry::DirectoryEntry( const std::string &NAME) :
      m_Name( ""),
      m_Type( Directory::e_Unknown),
      m_Directory()
    {
      // split name into path and filename
      storage::VectorND< 2, std::string> path_and_filename( File::SplitToPathAndFileName( NAME));
      m_Name = path_and_filename.Second();
      if( path_and_filename.First().empty())
      {
        path_and_filename.First() = ".";
      }
      m_Directory = util::ShPtr< Directory>( new Directory( path_and_filename.First()));

      // set attributes
      SetAttributes();
    }

    //! @brief construct from directory and name
    //! @param SP_DIR directory the file is contained
    //! @param NAME name of file, path or link
    DirectoryEntry::DirectoryEntry( const util::ShPtr< Directory> &SP_DIR, const std::string &NAME) :
      m_Name( NAME),
      m_Type( Directory::e_Unknown),
      m_Directory( SP_DIR)
    {
      // check that directory is defined
      BCL_Assert( SP_DIR.IsDefined(), "entry cannot be created from undefined directory pointer!");

      // set the attributes
      SetAttributes();
    }

    //! @brief construct from directory, name, and type
    //! @param SP_DIR directory the file is contained
    //! @param NAME name of file, path or link
    //! @param TYPE type of the directory entry
    DirectoryEntry::DirectoryEntry
    (
      const util::ShPtr< Directory> &SP_DIR,
      const std::string &NAME,
      const Directory::EntryType &TYPE
    ) :
      m_Name( NAME),
      m_Type( TYPE),
      m_Directory( SP_DIR)
    {
      // check that directory is defined
      BCL_Assert( SP_DIR.IsDefined(), "entry cannot be created from undefined directory pointer!");
    }

    //! @brief Clone function
    //! @return pointer to new DirectoryEntry
    DirectoryEntry *DirectoryEntry::Clone() const
    {
      return new DirectoryEntry( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DirectoryEntry::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to name
    //! @return name of entry
    const std::string &DirectoryEntry::GetName() const
    {
      return m_Name;
    }

    //! @brief get Directory
    const Directory &DirectoryEntry::GetDirectory() const
    {
      return *m_Directory;
    }

    //! @brief fullname including the path the entry is in
    //! @return fullpath and name
    std::string DirectoryEntry::GetFullName() const
    {
      return m_Directory->AppendFilename( m_Name);
    }

    //! @brief access to type
    //! @return type of directory entry
    Directory::EntryType DirectoryEntry::GetType() const
    {
      return m_Type;
    }

    //! @brief check if entry is of given type
    //! @param TYPE the desired entry type
    //! @return true if the type matches
    bool DirectoryEntry::IsType( const Directory::EntryType &TYPE) const
    {
      return m_Type == TYPE;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief check if the entry does exist
    //! @return true, if the entry exists, false otherwise
    bool DirectoryEntry::DoesExist() const
    {
      // check if directory is defined
      if( !m_Directory.IsDefined())
      {
        return false;
      }
      // create stat object "file_attributes" to hold the information about "FILENAME"
      struct stat file_attributes;

      const std::string path( GetFullName());
      // create int "status" and initialize with results of stat function
      // stat function returns 0 if the status information is retrieved. Otherwise returns -1 and errno is set
      // see http://www.cplusplus.com/reference/clibrary/cerrno/errno.html for more information about errno
      const int status( stat( path.c_str(), &file_attributes));

      // true if "status" is zero, meaning the information about "FILENAME" could be retrieved
      return status == 0;
    }

    //! @brief remove the directory entry
    //! @return true, if successful
    bool DirectoryEntry::Remove()
    {
      // check that file exists
      if( !DoesExist())
      {
        return false;
      }

      int error = ::remove( GetFullName().c_str());

      // update attributes
      SetAttributes();

      // end
      return error == 0;
    }

    //! @brief rename the directory entry
    //! @param NEW the new entry
    //! @param ALLOW_OVERWRITE whether to allow overwriting NEW if it already exists
    //! @return true, if successful
    bool DirectoryEntry::Rename( const DirectoryEntry &NEW, const bool &ALLOW_OVERWRITE)
    {
      if( !DoesExist())
      {
        return false;
      }

      if( !ALLOW_OVERWRITE && NEW.DoesExist())
      {
        return false;
      }

      int error = ::rename( GetFullName().c_str(), NEW.GetFullName().c_str());

      if( error == 0)
      {
        *this = NEW;
      }

      // update attributes
      SetAttributes();

      // return true if there were no errors
      return !error;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DirectoryEntry::Read( std::istream &ISTREAM)
    {
      // read member
      Serialize::Read( m_Name     , ISTREAM);
      Serialize::Read( m_Directory, ISTREAM);
      m_Type = Directory::e_Unknown;

      // set the attributes
      if( !m_Directory.IsDefined() || ( DoesExist() && !SetAttributes()))
      {
        BCL_MessageCrt( "unable to set attributes for file: " + m_Name);
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DirectoryEntry::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      Serialize::Write( m_Name     , OSTREAM, INDENT) << '\n';
      Serialize::Write( m_Directory, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set attributes
    //! sets all the entry attributes
    //! @return true if successful, false if there was an error (e.g file did not exist or other problems)
    bool DirectoryEntry::SetAttributes()
    {
      // create stat object "file_attributes" to hold the information about "FILENAME"
      struct stat file_attributes;

      // full path of entry
      const std::string path( GetFullName());

      // create int "status" and initialize with results of stat function
      // stat function returns 0 if the status information is retrieved. Otherwise returns -1 and errno is set
      // see http://www.cplusplus.com/reference/clibrary/cerrno/errno.html for more information about errno
      const int status( stat( path.c_str(), &file_attributes));

      // true if "status" is zero, meaning the information about path could be retrieved
      if( status != 0)
      {
        m_Type = Directory::e_Unknown;
        return false;
      }

      // find the type
      switch( file_attributes.st_mode & S_IFMT)
      {
        case S_IFDIR:
          m_Type = Directory::e_Dir;
          break;
        case S_IFREG:
          m_Type = Directory::e_File;
          break;
        case S_IFCHR:
        case S_IFMT:
        default:
          m_Type = Directory::e_Unknown;
          break;
      }

      // success
      return true;
    }

  } // namespace io
} // namespace bcl
