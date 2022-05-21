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
#include "io/bcl_io_directory.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <cstdio>
#include <sys/stat.h>
#if defined (_MSC_VER) or defined(__MINGW32__)
#  include <direct.h> // _mkdir
#  include <windows.h>
#  include <io.h>

// mimic unix opendir, closedir, readdir, rewinddir
struct dirent
{
  char *d_name;
  unsigned int attrib;
};

// local namespace to avoid exporting symbols
namespace
{
  //! @brief helper function to create an io::DirectoryEntry with info from a dirent, if available
  //! @param DIRECTORY the directory of interest
  //! @param ENTRY pointer to entry in that directory
  //! @return a new directory entry
  //! @note using information given on most systems about the filetype,
  //!       returned by readdir, expensive calls to stat (esp. on NFS) are avoided
  bcl::io::DirectoryEntry CreateDirectoryEntry
  (
    const bcl::util::ShPtr< bcl::io::Directory> &DIRECTORY,
    const dirent *ENTRY
  )
  {
    // null entry
    if( !ENTRY)
    {
      return bcl::io::DirectoryEntry();
    }
    else if( ENTRY->attrib == _A_SUBDIR)
    {
      // subdirectory
      return bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name, bcl::io::Directory::e_Dir);
    }
    // regular files, readonly files hidden files, system files
    else if
    (
      ENTRY->attrib == _A_ARCH
      || ENTRY->attrib == _A_NORMAL
      || ENTRY->attrib == _A_HIDDEN
      || ENTRY->attrib == _A_RDONLY
      || ENTRY->attrib == _A_SYSTEM
    )
    {
      return bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name, bcl::io::Directory::e_File);
    }
    // unknown file types, have to perform the usual stat command to determine type
    return bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name);
  }
}

struct DIR
{
    long                handle; // -1 for failed rewind
    struct _finddata_t  info;
    struct dirent       result; // d_name null iff first time
    char                *name;  // null-terminated char string
};

DIR *opendir( const char *name)
{
  DIR *dir = 0;

  if( name && name[ 0])
  {
    size_t base_length = strlen( name);
    const char *all = // search pattern must end with suitable wildcard
        strchr( "/\\", name[ base_length - 1]) ? "*" : "/*";

    if
    (
         ( dir = ( DIR *) malloc( sizeof *dir)) != 0
      && ( dir->name = ( char *) malloc( base_length + strlen( all) + 1)) != 0
    )
    {
      strcat( strcpy( dir->name, name), all);

      if( ( dir->handle = ( long) _findfirst( dir->name, &dir->info)) != -1)
      {
        dir->result.d_name = 0;
      }
      else // rollback
      {
        free( dir->name);
        free( dir);
        dir = 0;
      }
    }
    else // rollback
    {
      free( dir);
      dir   = 0;
      errno = ENOMEM;
    }
  }
  else
  {
    errno = EINVAL;
  }

  return dir;
}

int closedir( DIR *dir)
{
  int result = -1;

  if( dir)
  {
    if( dir->handle != -1)
    {
      result = _findclose( dir->handle);
    }

    free( dir->name);
    free( dir);
  }

  if( result == -1) // map all errors to EBADF
  {
    errno = EBADF;
  }

  return result;
}

struct dirent *readdir( DIR *dir)
{
  struct dirent *result = 0;

  if( dir && dir->handle != -1)
  {
    if( !dir->result.d_name || _findnext( dir->handle, &dir->info) != -1)
    {
      result         = &dir->result;
      result->d_name = dir->info.name;
      result->attrib = dir->info.attrib;
    }
  }
  else
  {
    errno = EBADF;
  }

  return result;
}

void rewinddir( DIR *dir)
{
  if( dir && dir->handle != -1)
  {
    _findclose( dir->handle);
    dir->handle = ( long) _findfirst( dir->name, &dir->info);
    dir->result.d_name = 0;
  }
  else
  {
    errno = EBADF;
  }
}
#elif defined(__GNUC__)
#  include <dirent.h>

// local namespace to avoid exporting symbols
namespace
{
#ifdef _DIRENT_HAVE_D_TYPE
  //! @brief helper function to get the type of file, if known
  //! @param DIRECTORY the directory of interest
  //! @param ENTRY pointer to entry in that directory
  //! @return a new directory entry
  //! @note using information given on most systems about the filetype,
  //!       returned by readdir, expensive calls to stat (esp. on NFS) are avoided
  bcl::io::DirectoryEntry CreateDirectoryEntry
  (
    const bcl::util::ShPtr< bcl::io::Directory> &DIRECTORY,
    const dirent *ENTRY
  )
  {
    // null entry
    if( !ENTRY)
    {
      return bcl::io::DirectoryEntry();
    }
    else if( ENTRY->d_type == DT_DIR)
    {
      // subdirectory
      return bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name, bcl::io::Directory::e_Dir);
    }
    else if( ENTRY->d_type == DT_REG)
    {
      // subdirectory
      return bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name, bcl::io::Directory::e_File);
    }

    // other files, links, sockets, etc.; have to determine the underlying file type (especially for links) using a
    // normal stat call.  This conditional is also reached if the file system always returns DT_UNKNOWN
    return bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name);
  }
#else // no d_type parameter
  //! @brief helper function to get the type of file, if known
  //! @param DIRECTORY the directory of interest
  //! @param ENTRY pointer to entry in that directory
  //! @return a new directory entry
  //! @note using information given on most systems about the filetype,
  //!       returned by readdir, expensive calls to stat (esp. on NFS) are avoided
  bcl::io::DirectoryEntry CreateDirectoryEntry
  (
    const bcl::util::ShPtr< bcl::io::Directory> &DIRECTORY,
    const dirent *ENTRY
  )
  {
    return ENTRY ? bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name) : bcl::io::DirectoryEntry();
  }
#endif

}
#else
// local namespace to avoid exporting symbols
namespace
{
  //! @brief helper function to get the type of file, if known
  //! @param DIRECTORY the directory of interest
  //! @param ENTRY pointer to entry in that directory
  //! @return a new directory entry
  //! @note using information given on most systems about the filetype,
  //!       returned by readdir, expensive calls to stat (esp. on NFS) are avoided
  bcl::io::DirectoryEntry CreateDirectoryEntry
  (
    const bcl::util::ShPtr< bcl::io::Directory> &DIRECTORY,
    const dirent *ENTRY
  )
  {
    // other file system; have to use the normal stat command
    return ENTRY ? bcl::io::DirectoryEntry( DIRECTORY, ENTRY->d_name) : bcl::io::DirectoryEntry();
  }
}
#endif

namespace bcl
{
  namespace io
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> Directory::s_Instance
    (
      GetObjectInstances().AddInstance( new Directory())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Directory::Directory() :
      m_Path( ".")
    {
    }

    //! @brief construct from path
    //! @param PATH path for directory
    Directory::Directory( const std::string &PATH) :
      m_Path( PATH)
    {
      CleanPath();
    }

    //! @brief Clone function
    //! @return pointer to new Directory
    Directory *Directory::Clone() const
    {
      return new Directory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Directory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get path
    //! @return the path for this directory
    const std::string &Directory::GetPath() const
    {
      return m_Path;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add directory to filename
    //! @param NAME the name of the file
    //! @return full name including this path
    std::string Directory::AppendFilename( const std::string &NAME) const
    {
      // creat path
      std::string path( m_Path);

      // only add path, if it actually contains a path, otherwise it points to current working directory
      // NOTE: It is CRITICAL that PATH_SEPARATOR not be added to path, otherwise calls to stat fail on _native_
      // windows systems (but not wine)
      if( !m_Path.empty() && !NAME.empty())
      {
        path += PATH_SEPARATOR;
      }
      path += NAME;

      // end
      return path;
    }

    //! @brief check if directory physically exist
    //! @return true, if is does exist, false if it does not exist as a directory (e.g. as file or something else)
    bool Directory::DoesExist() const
    {
      // create stat object "file_attributes" to hold the information about "FILENAME"
      struct stat file_attributes;

      // create int "status" and initialize with results of stat function
      // stat function returns 0 if the status information is retrieved. Otherwise returns -1 and errno is set
      // see http://www.cplusplus.com/reference/clibrary/cerrno/errno.html for more information about errno
      const int status( stat( m_Path.c_str(), &file_attributes));

      // true if "status" is zero, meaning the information about "FILENAME" could be retrieved
      if( status != 0)
      {
        // attributes of "FILENAME" could not be gotten therefore "FILENAME" does not exist so return false
        return false;
      }

      // check if the attributes mode indicate, that this is a directory
      const bool is_directory( ( file_attributes.st_mode & S_IFMT) == S_IFDIR);
      return is_directory;
    }

    //! @brief list content of directory
    //! @param DIRECTORY director object
    //! @return list of Directory entries
    storage::List< DirectoryEntry> Directory::ListEntries() const
    {
      // create directory that entrie will be associated with
      util::ShPtr< Directory> sp_dir( Clone());

      // entries
      storage::List< DirectoryEntry> entries;

      // open the directory
      DIR *directory( opendir( m_Path.c_str()));

      // check if directory could be opend
      if( directory == NULL)
      {
        BCL_MessageCrt( "unable to open path! " + m_Path);
        return entries;
      }

      // get all netried in the directory
      struct dirent *directory_entry( NULL);

      // loop as long as there is a directory
      while( ( directory_entry = readdir( directory)) != NULL)
      {
        // check that the directory entry is not '.' or '..'
        if( std::string( directory_entry->d_name) == "." || std::string( directory_entry->d_name) == "..")
        {
          continue;
        }
        // insert the entry
        entries.PushBack( CreateDirectoryEntry( sp_dir, directory_entry));
      }

      // close directory
      closedir( directory);

      // end
      return entries;
    }

    //! @brief list entries of the directory that match the given type and suffix
    //! @param TYPE type of entry (file/directory, etc., s_MaxTypes for any)
    //! @param PREFIX if given, only list entries that beging with the given string
    //! @param SUFFIX if given, only list entries with the given suffix
    //! @return list of Directory entries
    storage::List< DirectoryEntry> Directory::ListEntries
    (
      const EntryType &TYPE,
      const std::string &PREFIX,
      const std::string &SUFFIX
    ) const
    {
      // create directory that entrie will be associated with
      util::ShPtr< Directory> sp_dir( Clone());

      BCL_MessageVrb( "Listing entries in " + m_Path);

      // entries
      storage::List< DirectoryEntry> entries;

      // open the directory
      DIR *directory( opendir( m_Path.c_str()));

      // check if directory could be opend
      if( directory == NULL)
      {
        BCL_MessageCrt( "unable to open path! " + m_Path);
        return entries;
      }

      // get all netried in the directory
      struct dirent *directory_entry( NULL);

      // loop as long as there is a directory
      while( ( directory_entry = readdir( directory)) != NULL)
      {
        // create a new entry
        DirectoryEntry entry( CreateDirectoryEntry( sp_dir, directory_entry));

        // test that the type matches or that the type does not matter
        if( entry.GetType() != TYPE && TYPE != s_MaxTypes)
        {
          // type mismatch, continue
          continue;
        }

        // check that the directory entry is not '.' or '..'
        if( std::string( directory_entry->d_name) == "." || std::string( directory_entry->d_name) == "..")
        {
          continue;
        }

        // get a reference to the string
        const std::string &name( entry.GetName());

        // test for suffix or prefix mismatch
        if( !util::EndsWith( name, SUFFIX) || !util::StartsWith( name, PREFIX))
        {
          // suffix mismatch, continue
          continue;
        }

        // add the entry
        entries.PushBack( entry);
      }

      // close directory
      closedir( directory);

      // end
      return entries;
    }

    //! @brief list files of the directory that match the given type and suffix
    //! @param TYPE type of entry (file/directory, etc., s_MaxTypes for any)
    //! @param PREFIX if given, only list entries that beging with the given string
    //! @param SUFFIX if given, only list entries with the given suffix
    //! @param RECURSIVE true to auto-descend into sub-directories
    //! @note when using RECURSIVE, entries will be given in breadth-first order
    //! @return list of Directory entries
    storage::List< DirectoryEntry> Directory::ListFiles
    (
      const std::string &PREFIX,
      const std::string &SUFFIX,
      const bool &RECURSIVE
    ) const
    {
      // create directory that entrie will be associated with
      util::ShPtr< Directory> sp_dir( Clone());

      BCL_MessageVrb( "Listing entries in " + m_Path);

      // entries
      storage::List< DirectoryEntry> entries;

      // get all netried in the directory
      struct dirent *directory_entry( NULL);

      // a list of all sub-directories found
      util::ShPtrList< Directory> sub_directories( 1, util::ShPtr< Directory>( Clone()));
      while( !sub_directories.IsEmpty())
      {
        util::ShPtrList< Directory> new_subdirectories;
        for
        (
          util::ShPtrList< Directory>::const_iterator itr( sub_directories.Begin()), itr_end( sub_directories.End());
          itr != itr_end;
          ++itr
        )
        {
          // open the directory
          DIR *directory( opendir( ( *itr)->m_Path.c_str()));
          BCL_MessageVrb( "Opened directory or file at " + ( *itr)->m_Path);

          // check if directory could be opend
          if( directory == NULL)
          {
            BCL_MessageCrt( "unable to open path! " + ( *itr)->m_Path);
            return entries;
          }

          // loop as long as there is a directory
          while( ( directory_entry = readdir( directory)) != NULL)
          {
            // create a new entry
            DirectoryEntry entry( CreateDirectoryEntry( *itr, directory_entry));

            // get a reference to the string
            const std::string &name( entry.GetName());

            // check that the directory entry is not '.' or '..'
            if( name == "." || name == "..")
            {
              continue;
            }

            // test for directories
            if( entry.GetType() == e_Dir)
            {
              if( RECURSIVE)
              {
                new_subdirectories.PushBack( util::ShPtr< Directory>( new Directory( entry.GetFullName())));
              }
              continue;
            }

            // test for suffix or prefix mismatch
            if( !util::EndsWith( name, SUFFIX) || !util::StartsWith( name, PREFIX))
            {
              // suffix mismatch, continue
              continue;
            }

            // add the entry
            entries.PushBack( entry);
          }

          // close directory
          closedir( directory);
        }
        sub_directories = new_subdirectories;
      }

      // end
      return entries;
    }

    //! @brief create the directory
    //! @return true if it was successful, false if it did exist or could not be created
    bool Directory::Make()
    {
      // check if directory exists
      if( DoesExist())
      {
        return false;
      }

      // call the correct mkdir function
      int error( -1);
#if defined(__MINGW32__) || defined(__MINGW64__)
      error = mkdir( m_Path.c_str());
#elif defined (__GNUC__)
      umask( 0);
      error = mkdir( m_Path.c_str(), 0777);
#elif defined (_MSC_VER)
      error = _mkdir( m_Path.c_str());
#endif
      return error == 0;
    }

    //! @brief remove the directory
    //! @param RECURSIVE, default false
    //! @return true, if removal was successful, false if directory did not exist or removal failed (if it was not empty and recursive was not passed)
    bool Directory::Remove( const bool RECURSIVE)
    {
      bool clear_success( true);

      // clear directory, if recursive
      if( RECURSIVE)
      {
        clear_success &= Clear( RECURSIVE);
      }

      // check if recursive delete was successful, if requested
      if( clear_success)
      {
#if defined(__MINGW32__)
        int error = rmdir( m_Path.c_str());
#else
        int error = ::remove( m_Path.c_str());
#endif
        return error == 0;
      }

      // error recursively removing
      return false;
    }

    //! @brief clear the directory, without deleting it
    //! @param RECURSIVE, default false
    //! @return true, if removal of content was successful, false if directory did not exist or removal failed (if sub directories were not empty and recursive was not passed)
    bool Directory::Clear( const bool RECURSIVE)
    {
      // check if directory exists, to not falsely remove a file
      if( !DoesExist())
      {
        return false;
      }

      // list all directory entries
      storage::List< DirectoryEntry> entries( ListEntries());

      bool success( true);

      // iterate over entries and delete them
      for
      (
        storage::List< DirectoryEntry>::iterator itr( entries.Begin()), itr_end( entries.End());
        itr != itr_end && success;
        ++itr
      )
      {
        DirectoryEntry &current_entry( *itr);
        // skip unknown types
        if( !current_entry.DoesExist() || current_entry.IsType( Directory::e_Unknown))
        {
          continue;
        }
        // delete directory entries
        else if( current_entry.IsType( Directory::e_Dir))
        {
          // skip .  and .. and directories, if recursive is not requested
          if( current_entry.GetName() == "." || current_entry.GetName() == "..")
          {
            continue;
          }

          // create Directory object for that entry and remove recursively
          Directory current_dir( AppendFilename( current_entry.GetName()));
          success &= current_dir.Remove( RECURSIVE);
        }
        // regular files and links
        else
        {
          success &= current_entry.Remove();
        }
      }

      // end
      return success;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Directory::Read( std::istream &ISTREAM)
    {
      // read member
      Serialize::Read( m_Path, ISTREAM);

      // clean the path
      CleanPath();

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Directory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      Serialize::Write( m_Path, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a directory
    //! @param DIRECTORY name of directory to be created
    //! @return the Directory associated with the created directory, check if is exists to confirm successful creation
    Directory Directory::MkDir( const std::string &DIRECTORY)
    {
      // create the object
      Directory dir( DIRECTORY);

      // return dir is it exists
      if( dir.DoesExist())
      {
        return DIRECTORY;
      }

      // create directory
      dir.Make();

      // return the directory
      return dir;
    }

    //! @brief clean directory
    //! changes empty to '.' and removes additional '/' or '\\'
    void Directory::CleanPath()
    {
      if( m_Path.empty())
      {
        m_Path = ".";
      }
      // last position
      size_t last( m_Path.length() - 1);

      // remove additional slashes in the end
      while( last > 0 && ( m_Path[ last] == '/' || m_Path[ last] == '\\'))
      {
        m_Path.erase( last, 1);
        last = m_Path.length() - 1;
      }
    }

  } // namespace io
} // namespace bcl
