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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "io/bcl_io_directory.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_directory.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoDirectory :
    public ExampleInterface
  {
  public:

    ExampleIoDirectory *Clone() const
    { return new ExampleIoDirectory( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constuctor
      io::Directory default_dir;

      // construct from path
      io::Directory dir_from_path( AddExampleOutputPathToFilename( default_dir, ""));

      // construct from input file path.
      io::Directory dir_from_input_path
      (
        GetExamples().GetExamplePath() + s_ExampleInputFolderName + PATH_SEPARATOR
      );

      // construct from path
      io::Directory dir_test( AddExampleOutputPathToFilename( default_dir, "test"));

      // clone
      util::ShPtr< io::Directory> sp_dir( dir_test.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( sp_dir->GetClassIdentifier(), GetStaticClassName( default_dir));

      // path
      BCL_ExampleCheck( dir_test.GetPath(), AddExampleOutputPathToFilename( default_dir, "test"));

    ////////////////
    // operations //
    ////////////////

      // append filename
      BCL_ExampleCheck( dir_from_path.AppendFilename( "test.pdb"), AddExampleOutputPathToFilename( dir_from_path, "test.pdb"));

      // check if directory exists
      BCL_ExampleCheck( dir_test.DoesExist(), false);

      // create directory
      const bool make_success( dir_test.Make());
      BCL_ExampleIndirectCheck( make_success, true, "Make directory");

      // remove directory
      const bool remove_success( dir_test.Remove());
      BCL_ExampleIndirectCheck( remove_success, true, "Remove directory");

      // mkdir
      io::Directory mk_test( io::Directory::MkDir( AddExampleOutputPathToFilename( default_dir, "test2")));
      BCL_ExampleIndirectCheck( mk_test.DoesExist(), true, "MkDir");

      // create subdirectory and file
      io::Directory sub_dir( mk_test.AppendFilename( "test3"));
      const bool make_subdir_success( sub_dir.Make());
      BCL_ExampleIndirectCheck( make_subdir_success, true, "make sub directory");

      // determine the number of items in the newly-created directory; this may vary depending on system and environment
      // so don't check it against a particular number
      const size_t original_files_sub_dir( sub_dir.ListEntries().GetSize());
      BCL_ExampleIndirectCheck
      (
        sub_dir.ListEntries( io::Directory::e_Dir).GetSize(),
        0,
        "make sub directory should not have created any sub-sub directories"
      );

      // file
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, sub_dir.AppendFilename( "testfile.txt"));
      write << "hello";
      io::File::CloseClearFStream( write);

      BCL_ExampleCheck( sub_dir.ListEntries( io::Directory::e_File, "", ".txt").GetSize(), 1);
      BCL_ExampleCheck( sub_dir.ListEntries( io::Directory::e_File, "test").GetSize(), 1);

      // check that the new file is found
      if( BCL_ExampleIndirectCheck( sub_dir.ListEntries().GetSize(), original_files_sub_dir + 1, "new files are seen"))
      {
        // check that the GetDirectory().GetPath() on a directory element is identical to the path of the directory
        BCL_ExampleCheck
        (
          sub_dir.ListEntries().FirstElement().GetDirectory().GetPath(),
          sub_dir.GetPath()
        );
      }

      // try non recursive clear
      const bool clear_mk_test( mk_test.Clear());
      BCL_ExampleIndirectCheck( clear_mk_test, false, "MkDir dir clear");

      // try non recursive removal
      const bool remove_mk_test( mk_test.Remove());
      BCL_ExampleIndirectCheck( remove_mk_test, false, "MkDir dir remove");

      // recursive removal (test indirectly recursive clear)
      const bool remove_recursive_mk_test( mk_test.Remove( true));
      BCL_ExampleIndirectCheck( remove_recursive_mk_test, true, "MkDir dir recursive remove");

      // list all directory entries
      const storage::List< io::DirectoryEntry> entries( dir_from_input_path.ListEntries());
      size_t number_dirs( 0);
      for( storage::List< io::DirectoryEntry>::const_iterator itr( entries.Begin()), itr_end( entries.End()); itr != itr_end; ++itr)
      {
        if( itr->IsType( io::Directory::e_Dir) && itr->GetName()[ 0] != '.') // only count non-special directories
        {
          ++number_dirs;
        }
      }
      BCL_ExampleIndirectCheck( number_dirs, s_NumberExampleInputTypes, "ListEntries number of directories");

    //////////////////////
    // input and output //
    //////////////////////

      // read and write bcl object
      WriteBCLObject( dir_from_path);
      io::Directory dir_read;
      ReadBCLObject( dir_read);

      // check if the path is correct
      BCL_ExampleIndirectCheck( dir_read.GetPath(), dir_from_path.GetPath(), "write and read");

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoDirectory

  const ExampleClass::EnumType ExampleIoDirectory::s_Instance
  (
    GetExamples().AddEnum( ExampleIoDirectory())
  );

} // namespace bcl
