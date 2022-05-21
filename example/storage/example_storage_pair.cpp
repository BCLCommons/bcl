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
#include "storage/bcl_storage_pair.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_pair.cpp
  //!
  //! @author heinzes1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStoragePair :
    public ExampleInterface
  {
  public:

    ExampleStoragePair *Clone() const
    {
      return new ExampleStoragePair( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      const util::Format format_pair;

      // test default constructor
      storage::Pair< int, std::string> def_constructor;

      // test construction from two values and accessing the vaules by First and Second functions
      BCL_MessageStd( "construct from two values  " + util::Format()( 123) + " and LYS");
      storage::Pair< int, std::string> residue( 123, "LYS");
      BCL_MessageStd( "test access to elements  " + util::Format()( 123) + " and LYS");
      BCL_MessageStd( "residue: first element:  " + util::Format()( residue.First()));
      BCL_MessageStd( "residue: second element: " + residue.Second());
      BCL_Example_Check
      (
        residue.First() == 123 && residue.Second() == "LYS",
        "construction with two values and accessing with First and Second functions did not work"
      );

      // test construction from a std::pair
      BCL_MessageStd( "construct from std::pair  " + util::Format()( 456) + " and GLY");
      storage::Pair< int, std::string> from_std( std::pair< int, std::string>( 456, "GLY"));
      BCL_MessageStd( "from_std:  first element:  " + util::Format()( from_std.First()));
      BCL_MessageStd( "from_std: second element: " + from_std.Second());
      BCL_Example_Check
      (
        from_std.First() == 456 && from_std.Second() == "GLY",
        "construction from std::pair and accessing with First and Second functions did not work"
      );

      // test copy constructor
      BCL_MessageStd( "copy constructor: copy from_std  ");
      storage::Pair< int, std::string> copy( from_std);
      BCL_MessageStd( "copy: first element:  " + util::Format()( copy.First()));
      BCL_MessageStd( "copy: second element: " + copy.Second());
      BCL_Example_Check
      (
        copy.First() == 456 && copy.Second() == "GLY",
        "copy constructor and accessing with First and Second functions did not work"
      );

      // test virtual copy constructor (Clone)
      BCL_MessageStd( "Clone constructor: copy from_std  ");
      util::ShPtr< storage::Pair< int, std::string> > sp_to_base( from_std.Clone());
      BCL_MessageStd( "copy: first element:  " + util::Format()( sp_to_base->First()));
      BCL_MessageStd( "copy: second element: " + sp_to_base->Second());
      BCL_Example_Check
      (
        sp_to_base->First() == 456 && sp_to_base->Second() == "GLY",
        "Clone constructor and accessing with First and Second functions did not work"
      );

      // test GetClassIdentifier function
      BCL_MessageStd( "GetClassIdentifier: " + residue.GetClassIdentifier());
      BCL_Example_Check
      (
        residue.GetClassIdentifier() == "bcl::storage::Pair<int,std::string>",
        "GetClassIdentifier does not return the correct name bcl::storage::Pair<int,std::string>"
      );

      // test GetStaticClassName function
      BCL_MessageStd( "GetStaticClassName: " + GetStaticClassName( residue));
      BCL_Example_Check
      (
        GetStaticClassName( residue) == std::string( "bcl::storage::Pair<int,std::string>"),
        "GetStaticClassName does not return bcl::storage::Pair<int,std::string>"
      );

      // test data manipulation of First
      BCL_MessageStd( "test data manipulation of First: ");
      BCL_MessageStd( "residue: first element: before " + util::Format()( residue.First()));
      residue.First() = 789;
      BCL_MessageStd( "residue: first element: after " + util::Format()( residue.First()));
      BCL_Example_Check
      (
        residue.First() == 789, "data manipulation of First did not work correctly"
      );

      // test data manipulation of Second
      BCL_MessageStd( "test data manipulation of Second: ");
      BCL_MessageStd( "residue: Second element: before " + util::Format()( residue.Second()));
      residue.Second() = "VAL";
      BCL_MessageStd( "residue: Second element: after " + util::Format()( residue.Second()));
      BCL_Example_Check
      (
        residue.Second() == "VAL", "data manipulation of Second did not work correctly"
      );

      // test equal operator with two bcl::storage::Pair
      BCL_MessageStd( "test operator = with two bcl::storage::Pair: residue = from_std");
      BCL_MessageStd( "from_std:  first element:  " + util::Format()( from_std.First()));
      BCL_MessageStd( "from_std: second element: " + from_std.Second());
      BCL_MessageStd( "residue: first element:  " + util::Format()( residue.First()));
      BCL_MessageStd( "residue: second element: " + residue.Second());
      residue = from_std;
      BCL_MessageStd( "after: residue: first element:  " + util::Format()( residue.First()));
      BCL_MessageStd( "after: residue: second element: " + residue.Second());
      BCL_Example_Check
      (
        residue.Second() == from_std.Second() && residue.First() == from_std.First(), "operator = did not work properly"
      );

      // test equal operator with  bcl::storage::Pair and std::pair
      BCL_MessageStd( "test operator = with  bcl::storage::Pair and std::pair: ");
      BCL_MessageStd( "residue: first element:  " + util::Format()( residue.First()));
      BCL_MessageStd( "residue: second element: " + residue.Second());
      residue = std::pair< int, std::string>( 123, "PHE");
      BCL_MessageStd( "after: residue: first element:  " + util::Format()( residue.First()));
      BCL_MessageStd( "after: residue: second element: " + residue.Second());
      BCL_Example_Check
      (
        residue.Second() == "PHE" && residue.First() == 123, "operator = did not work properly"
      );

      // test output
      BCL_MessageStd( "this is the complete Pair: " + util::Format()( residue));

      // test operator == in true scenario
      BCL_MessageStd
      (
        "test operator == with true scenario: " + util::Format()( copy == from_std)
      );
      BCL_Example_Check
      (
        copy == from_std, "operator == returns false when it should return true"
      );

      // test operator == in false scenario
      BCL_MessageStd
      (
        "test operator == with false scenario: " + util::Format()( copy == residue)
      );
      BCL_Example_Check
      (
        !( copy == residue), "operator == returns true when it should return false"
      );

      // test operator != in true scenario
      BCL_MessageStd
      (
        "test operator != with true scenario: " + util::Format()( copy != residue)
      );
      BCL_Example_Check
      (
        copy != residue, "operator != returns false when it should return true"
      );

      // test operator != in false scenario
      BCL_MessageStd
      (
        "test operator != with false scenario: " + util::Format()( copy != from_std)
      );
      BCL_Example_Check
      (
        !( copy != from_std), "operator != returns true when it should return false"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStoragePair

  const ExampleClass::EnumType ExampleStoragePair::s_Instance
  (
    GetExamples().AddEnum( ExampleStoragePair())
  );

} // namespace bcl
