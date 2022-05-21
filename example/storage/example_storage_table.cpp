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
#include "storage/bcl_storage_table.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_table.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageTable :
    public ExampleInterface
  {
  public:

    ExampleStorageTable *Clone() const
    {
      return new ExampleStorageTable( *this);
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
      // initialize table header
      storage::TableHeader table_header( storage::Vector< std::string>::Create( "score", "rmsd", "rmsd100"));
      util::ShPtr< storage::TableHeader> sp_table_header( table_header.Clone());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      storage::Table< double> table_default;
      BCL_ExampleCheck( storage::Table< double>().IsEmpty(), true);

      // constructor from a table_header
      storage::Table< double> table_a( table_header);
      BCL_ExampleCheck( storage::Table< double>( table_header).GetHeader(), table_header);

      // construct from a ShPtr< TableHeader>
      storage::Table< double> table_b( sp_table_header);
      BCL_ExampleCheck( storage::Table< double>( sp_table_header).GetHeader(), *sp_table_header);

      // copy constructor
      BCL_ExampleCheck( storage::Table< double>( table_a), table_a);

      // contruct using clone operator
      BCL_ExampleCheck( *util::ShPtr< storage::Table< double> >( table_a.Clone()), table_a);

    /////////////////
    // data access //
    /////////////////

      // test class identifier
      BCL_ExampleCheck( table_a.GetClassIdentifier(), GetStaticClassName( table_a));

      // check GetSize
      // test static class name
      BCL_MessageStd( "Checking GetSize function");
      BCL_ExampleCheck( storage::Table< double>().GetSize(), 0);

      // inserting elements to table_a
      BCL_MessageStd( "Inserting 5 elements into table_a");
      table_a.InsertRow( "0", storage::Vector< double>::Create( -500, 5.0, 6.0));
      table_a.InsertRow( "1", storage::Vector< double>::Create( -600, 3.2, 5.0));
      table_a.InsertRow( "2", storage::Vector< double>::Create( -700, 2.8, 4.0));
      table_a.InsertRow( "3", storage::Vector< double>::Create( -400, 5.6, 3.0));
      table_a.InsertRow( "4", storage::Vector< double>::Create( -300, 5.4, 2.0));

      // check that all insertions were correct
      BCL_ExampleIndirectCheck( table_a.GetSize(), 5, "GetSize() after 5 calls to InsertRow");

      // check construct from iterator range
      storage::Table< double> table_c( table_a.Begin(), table_a.End());
      BCL_ExampleCheck( storage::Table< double>( table_a.Begin(), table_a.End()), table_a);

      // create the row names list
      storage::List< std::string> row_names_list;
      row_names_list.PushBack( "0");
      row_names_list.PushBack( "1");
      row_names_list.PushBack( "2");
      row_names_list.PushBack( "3");
      row_names_list.PushBack( "4");

      // check the GetRowNames function
      BCL_ExampleCheck( table_a.GetRowNames(), row_names_list);

      // check HasRow function
      BCL_MessageStd( "Checking HasRow with an existing row name");
      BCL_ExampleCheck( table_a.HasRow( "1"), true);

      // check HasRow function
      BCL_MessageStd( "Checking HasRow with a non-existing row name");
      BCL_ExampleCheck( !table_a.HasRow( "5"), true);

      // check the RowIndex function
      BCL_ExampleCheck( table_a.RowIndex( "1"), 1);

    ///////////////
    // operators //
    ///////////////

      // check the tag 0 gives the correct
      BCL_ExampleCheck( table_a[ "0"].GetData(), storage::Vector< double>::Create( -500, 5.0, 6.0));

    ////////////////
    // operations //
    ////////////////

      // remove the element with tag 4
      table_a.RemoveRow( "4");

      // assert the removal was successful
      BCL_ExampleIndirectCheck( table_a.HasRow( "4"), false, "table_a.RemoveRow( \"4\")");

      // inserting elements to table_a
      table_b.InsertRow( "5", storage::Vector< double>::Create( -300,  3.2,  5.0));
      table_b.InsertRow( "6", storage::Vector< double>::Create( -200,  2.7,  4.0));
      table_b.InsertRow( "7", storage::Vector< double>::Create( -100,  3.6, 11.0));
      table_b.InsertRow( "8", storage::Vector< double>::Create(    0, 13.4, 12.0));

      // check that all insertions were correct
      BCL_ExampleIndirectCheck( table_b.GetSize(), 4, "InsertRow called 4 times");

      // check the InsertRows function
      storage::Vector< std::string> row_names_vector
      (
        storage::Vector< std::string>::Create
        (
          "row_a",
          "row_b",
          "row_c",
          "this_the_name_of_the_last_row"
        )
      );
      table_c.InsertRows( row_names_vector);
      BCL_ExampleIndirectCheck( table_c.GetSize(), 9, "InsertRows");

      // test the append function
      BCL_MessageStd( "Testing the append function by appending table_b to table_a");
      table_a.Append( table_b);
      BCL_ExampleIndirectCheck( table_a.GetSize(), 8, "Append");

      // initialize format object
      util::Format format1;
      format1.W( 10).FFP( 3).R();

      // initialize write stream
      io::OFStream write;
      // test formatted output
      BCL_MessageStd( "testing formatted output to file for different sorts");
      // open file to output the table in formatted way
      const std::string formatted_filename( AddExampleOutputPathToFilename( table_default, "table_sorted.txt"));
      BCL_ExampleMustOpenOutputFile( write, formatted_filename);

      // sort table_a by score column
      BCL_MessageStd( "Sorting table_a by score column");
      table_a.SortByColumn( "score");
      write << "Table sorted by score column" << '\n';
      table_a.WriteFormatted( write, format1);

      storage::Table< double>::const_iterator row_itr_a( table_a.Begin());
      storage::Table< double>::const_iterator row_itr_b( table_a.Begin());
      storage::Table< double>::const_iterator row_itr_end( table_a.End());
      ++row_itr_b;
      for
      ( ;
        row_itr_a != row_itr_end && row_itr_b != row_itr_end;
        ++row_itr_a, ++row_itr_b
      )
      {
        BCL_Example_Check
        (
          row_itr_a->Second()[ "score"] <= row_itr_b->Second()[ "score"],
          "The scores are not in correct order " + util::Format()( row_itr_a->Second()[ "score"]) +
          " vs " + util::Format()( row_itr_b->Second()[ "score"])
        );
      }

      // sort table_a by rmsd column
      BCL_MessageStd( "Sorting table_a by rmsd column");
      table_a.SortByColumn( "rmsd");
      write << "Table sorted by rmsd column" << '\n';
      table_a.WriteFormatted( write, format1);
      row_itr_a = table_a.Begin(); row_itr_b = table_a.Begin(); row_itr_end = table_a.End();
      ++row_itr_b;
      for
      (
        ;
        row_itr_a != row_itr_end && row_itr_b != row_itr_end;
        ++row_itr_a, ++row_itr_b
      )
      {
        BCL_Example_Check
        (
          row_itr_a->Second()[ "rmsd"] <= row_itr_b->Second()[ "rmsd"],
          "The rmsd are not in correct order " + util::Format()( row_itr_a->Second()[ "rmsd"]) +
          " vs " + util::Format()( row_itr_b->Second()[ "rmsd"])
        );
      }

      // sort table_a by row names
      BCL_MessageStd( "Sorting table_a by row names");
      table_a.SortByRowName( std::less< std::string>());
      write << "Table sorted by row name" << '\n';
      table_a.WriteFormatted( write, format1);
      row_itr_a = table_a.Begin(); row_itr_b = table_a.Begin(); row_itr_end = table_a.End();
      ++row_itr_b;
      for
      ( ;
        row_itr_a != row_itr_end && row_itr_b != row_itr_end;
        ++row_itr_a, ++row_itr_b
      )
      {
        BCL_Example_Check
        (
          row_itr_a->First() <= row_itr_b->First(),
          "The row names are not in correct order " + row_itr_a->First() + " vs " + row_itr_b->First()
        );
      }
      write << "Table transposed" << '\n';
      storage::Table< double> transposed_table( table_a.GetTransposedTable());
      transposed_table.WriteFormatted( write, format1);

      // reset stream
      io::File::CloseClearFStream( write);

    //////////////////////
    // input and output //
    //////////////////////

      // initialize read stream
      io::IFStream read;

      // test class input output to file
      BCL_MessageStd( "testing read and write functionalities");

      // write object
      WriteBCLObject( table_a);
      // read the object back in
      storage::Table< double> table_read;
      ReadBCLObject( table_read);

      // compare the objects
      BCL_MessageStd( "compare written and read objects");
      BCL_ExampleIndirectCheck( table_a, table_read, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      // test get widest row name
      BCL_MessageStd( "Testing GetWidestRowName function");
      BCL_ExampleCheck( table_c.GetWidestRowName(), "this_the_name_of_the_last_row");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageTable

  const ExampleClass::EnumType ExampleStorageTable::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageTable())
  );

} // namespace bcl

