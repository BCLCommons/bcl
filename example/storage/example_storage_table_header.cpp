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
#include "storage/bcl_storage_table_header.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_table_header.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageTableHeader :
    public ExampleInterface
  {
  public:

    ExampleStorageTableHeader *Clone() const
    {
      return new ExampleStorageTableHeader( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // use default constructor
      BCL_MessageStd( "Calling default constructor");
      storage::TableHeader table_header_default;
      BCL_Example_Check
      (
        table_header_default.IsEmpty(),
        "The default table header should be empty but it has a size of " + util::Format()( table_header_default.GetSize())
      );

      // construct vector of string
      storage::Vector< std::string> column_names_a( storage::Vector< std::string>::Create( "score", "rmsd", "rmsd100"));
      storage::Vector< std::string> column_names_b( storage::Vector< std::string>::Create( "rmsd", "score", "rmsd100"));

      // construct table header from the columns
      BCL_MessageStd( "Calling constructors with vector of strings");
      storage::TableHeader table_header_a( column_names_a);
      BCL_Example_Check
      (
        table_header_a == column_names_a,
        "The table_header_a differs from column_names_a\n" + util::Format()( table_header_a) +
          "\nvs\n" + util::Format()( column_names_a)
      );
      storage::TableHeader table_header_b( column_names_b);
      BCL_Example_Check
      (
        table_header_b == column_names_b,
        "The table_header_b differs from column_names_b\n" + util::Format()( table_header_b) +
          "\nvs\n" + util::Format()( column_names_b)
      );

      // construct using copy constructor
      BCL_MessageStd( "Calling copy constructor");
      storage::TableHeader table_header_a_copy( table_header_a);
      BCL_Example_Check
      (
        table_header_a_copy == table_header_a,
        "The table_header_a_copy differs from table_header_a_copy\n" + util::Format()( table_header_a_copy) +
          "\nvs\n" + util::Format()( table_header_a_copy)
      );

      // contruct using clone operator
      BCL_MessageStd( "Calling clone");
      util::ShPtr< storage::TableHeader> sp_table_header( table_header_a.Clone());
      BCL_Example_Check
      (
        *sp_table_header == table_header_a,
        "The sp_table_header differs from table_header_a\n" + util::Format()( *sp_table_header) +
          "\nvs\n" + util::Format()( table_header_a)
      );

    /////////////////
    // data access //
    /////////////////

      // test class identifier
      BCL_MessageStd( "The class identifier is: " + table_header_a.GetClassIdentifier());
      BCL_Example_Check
      (
        table_header_a.GetClassIdentifier() == "bcl::storage::TableHeader",
        "The class identifier should be  bcl::storage::TableHeader"
      );

      // test static class name
      BCL_MessageStd( "The static class name is: " + GetStaticClassName( table_header_a));
      BCL_Example_Check
      (
        GetStaticClassName( table_header_a) == "bcl::storage::TableHeader",
        "The static class name should be  bcl::storage::TableHeader"
      );

    ///////////////
    // operators //
    ///////////////

      // test operator [] for table_header_a
      BCL_MessageStd( "testing operator[] for table_header_a");
      BCL_Example_Check
      (
        table_header_a[ "score"] == 0, "score should have the index 0 not " + util::Format()( table_header_a[ "score"])
      );
      BCL_Example_Check
      (
        table_header_a[ "rmsd"] == 1, "rmsd should have the index 1 not " + util::Format()( table_header_a[ "rmsd"])
      );

      // test operator [] for table_header_b
      BCL_MessageStd( "testing operator[] for table_header_b");
      BCL_Example_Check
      (
        table_header_b[ "rmsd"] == 0, "rmsd should have the index 0 not " + util::Format()( table_header_b[ "rmsd"])
      );
      BCL_Example_Check
      (
        table_header_b[ "score"] == 1, "score should have the index 1 not " + util::Format()( table_header_b[ "score"])
      );

      BCL_MessageStd( "Calling operator =");
      table_header_b = table_header_a;
      BCL_Example_Check
      (
        table_header_a == table_header_b,
        "The table_header_b is not correctly assigned to table_header_a!!\ntable_header_a\n" +
        util::Format()( table_header_a) + "table_header_b\n" +
        util::Format()( table_header_b)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // initialize read write stream
      io::IFStream read;
      io::OFStream write;

      // test formatted output
      BCL_MessageStd( "testing formatted output to file");
      // open file to output the table in formatted way
      const std::string formatted_filename( AddExampleOutputPathToFilename( table_header_default, "table_header_formatted.txt"));
      BCL_ExampleMustOpenOutputFile( write, formatted_filename);

      // write formatted with default format
      BCL_MessageStd( "Outputting with default width( 10)");
      table_header_a.WriteFormatted( write);
      write << '\n';

      // write formatted with width of 10
      BCL_MessageStd( "Outputting with width of 15");
      table_header_a.WriteFormatted( write, util::Format().W( 15));
      write << '\n';

      // write formatted with width of 15
      BCL_MessageStd( "Outputting with width of 20");
      table_header_a.WriteFormatted( write, util::Format().W( 20));
      write << '\n';

      // reset stream
      io::File::CloseClearFStream( write);

      // test class input output to file
      BCL_MessageStd( "testing read and write functionalities");
      WriteBCLObject( table_header_a);
      storage::TableHeader table_header_read;
      ReadBCLObject( table_header_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_ExampleIndirectCheck( table_header_a, table_header_read, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageTableHeader

  const ExampleClass::EnumType ExampleStorageTableHeader::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageTableHeader())
  );

} // namespace bcl

