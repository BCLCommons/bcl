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
#include "io/bcl_io_directory_entry.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io_directory_entry.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIoDirectoryEntry :
    public ExampleInterface
  {
  public:

    ExampleIoDirectoryEntry *Clone() const
    { return new ExampleIoDirectoryEntry( *this);}

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

      // default constructor
      io::DirectoryEntry dir_ent_default;

      // other names
      std::string path_name( AddExampleOutputPathToFilename( dir_ent_default, ""));
      path_name.erase( path_name.length() - 1, 1);
      const std::string filename( "testfile.txt");
      const std::string fullfilename( AddExampleOutputPathToFilename( dir_ent_default, filename));

      // directory
      util::ShPtr< io::Directory> sp_dir( new io::Directory( path_name));

      // construct from name
      io::DirectoryEntry dir_ent_from_name( fullfilename);

      // construct from dir and name
      io::DirectoryEntry dir_ent_dir_name( sp_dir, filename);

      // clone
      util::ShPtr< io::DirectoryEntry> sp_dir_ent( dir_ent_dir_name.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( sp_dir_ent->GetClassIdentifier(), GetStaticClassName( dir_ent_default));

      // name
      BCL_ExampleCheck( dir_ent_dir_name.GetName(), filename);
      BCL_ExampleCheck( dir_ent_from_name.GetName(), filename);

      // directory
      BCL_ExampleCheck( dir_ent_dir_name.GetDirectory().GetPath(), path_name);
      BCL_ExampleCheck( dir_ent_from_name.GetDirectory().GetPath(), path_name);

      // full name
      BCL_ExampleCheck( dir_ent_dir_name.GetFullName(), fullfilename);
      BCL_ExampleCheck( dir_ent_from_name.GetFullName(), fullfilename);

      // type
      BCL_ExampleCheck( dir_ent_dir_name.GetType(), io::Directory::e_Unknown);
      BCL_ExampleCheck( dir_ent_from_name.GetType(), io::Directory::e_Unknown);

    ////////////////
    // operations //
    ////////////////

      // is type
      BCL_ExampleCheck( dir_ent_dir_name.IsType( io::Directory::e_Unknown), true);
      BCL_ExampleCheck( dir_ent_from_name.IsType( io::Directory::e_File), false);

      // does exist
      BCL_ExampleCheck( dir_ent_dir_name.DoesExist(), false);
      BCL_ExampleCheck( dir_ent_from_name.DoesExist(), false);

      // create the file
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, fullfilename);
      write << "hello";
      io::File::CloseClearFStream( write);

      // does exist
      BCL_ExampleCheck( dir_ent_dir_name.DoesExist(), true);
      BCL_ExampleCheck( dir_ent_from_name.DoesExist(), true);

      // update the attributes
      BCL_ExampleCheck( dir_ent_dir_name.SetAttributes(), true);
      BCL_ExampleCheck( dir_ent_from_name.SetAttributes(), true);

      BCL_ExampleIndirectCheck( dir_ent_dir_name.IsType( io::Directory::e_File), true, "setattributes");
      BCL_ExampleIndirectCheck( dir_ent_from_name.IsType( io::Directory::e_File), true, "setattributes");

      // remove
      BCL_ExampleCheck( dir_ent_dir_name.Remove(), true);
      BCL_ExampleCheck( dir_ent_from_name.Remove(), false);

      // update attributes
      BCL_ExampleCheck( dir_ent_dir_name.SetAttributes(), false);
      BCL_ExampleCheck( dir_ent_from_name.SetAttributes(), false);
      BCL_ExampleIndirectCheck( dir_ent_dir_name.IsType( io::Directory::e_Unknown), true, "setattributes");
      BCL_ExampleIndirectCheck( dir_ent_from_name.IsType( io::Directory::e_Unknown), true, "setattributes");

    //////////////////////
    // input and output //
    //////////////////////

      // write and read
      WriteBCLObject( dir_ent_from_name);
      io::DirectoryEntry dir_ent_read;
      ReadBCLObject( dir_ent_read);

      // check
      BCL_ExampleIndirectCheck( dir_ent_read.GetDirectory().GetPath(), dir_ent_from_name.GetDirectory().GetPath(), "read and write");
      BCL_ExampleIndirectCheck( dir_ent_read.GetName(), dir_ent_from_name.GetName(), "read and write");

    //////////////////////
    // helper functions //
    //////////////////////

      // set attributes was tested above

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIoDirectoryEntry

  const ExampleClass::EnumType ExampleIoDirectoryEntry::s_Instance
  (
    GetExamples().AddEnum( ExampleIoDirectoryEntry())
  );

} // namespace bcl
